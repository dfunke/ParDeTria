import numpy as np

from scipy.stats import norm

import matplotlib.pylab as plt
import matplotlib.patches as mpatches
import matplotlib.lines as mlines

from matplotlib.ticker import FormatStrFormatter
from mpl_toolkits.mplot3d import Axes3D

import os

def ensure_dir(f):
    d = os.path.dirname(f)
    if not os.path.exists(d):
        os.makedirs(d)

def getFig(title, xlabel, ylabel, zlabel = None, grid = 'both'):
    if not zlabel:
        fig, ax = plt.subplots()
        if not grid is 'none':
            ax.grid(b=True, which='both', axis=grid)
    else:
        fig, ax = plt.subplots(subplot_kw={'projection' :'3d'})
        ax.set_zlabel(zlabel)

    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)

    return fig, ax


def closeFig(file):

    if not file.startswith('plot/output/'):
        file = 'plot/output/' + file

    if not '.' in file:
        file = file + '.png'

    ensure_dir(file)

    plt.savefig(file)

    plt.close()

dir = 'data'
impl = ['cgal', 'tria', 'root']
labels =  {'cgal' : 'CGAL', 'tria' : 'Triangle', 'root' : 'ROOT'}
colors  = {'cgal' : 'blue', 'tria' : 'green', 'root' : 'red'}
markers = {'cgal' : 'x', 'tria' : 'D', 'root' : 'o'}

funcNames = [  'sin'
              ,'sinsin'
              ,'cubic'
              ,'sqrtsincos']

funcLabels = [ r'$\sin(4 xy)$'
              ,r'$ \frac{\sin{x}}{x} \cdot \frac{\sin{y}}{y} + 0.2$'
              ,r'$x^3 - 3x + y^3 -3y$'
              ,r'$\sqrt{\sin{x}+\cos{y}}$']


timings = {}

print("Reading in timing data")
for tag in impl:
    print("\t%s" % tag)
    timings[tag] = np.genfromtxt(os.path.join(dir, 'timings_%s.csv' % tag),
                                 names=True)

#setup unique Ns and Fs
N = np.unique(timings[impl[0]]['n']).astype(int)
F = np.unique(timings[impl[0]]['f']).astype(int)

gen = {}
ip = {}
tt = {}

print("Processing timing data")
for tag in impl:
    print("\t%s" % tag)
    gen[tag] = []
    ip[tag] = []
    tt[tag] = []
    for n in N:
        dat = timings[tag][timings[tag]['n'] == n]
        gen[tag].append(dat['gen'])
        ip[tag].append(dat['ip'])
        tt[tag].append(dat['gen'] + 1e4 * dat['ip'])

print("Generating timing plots")

fig, ax = getFig("Triangulation Runtime", r"$n$ sampling points", r"$t$ [ns]", None, "x")
start = .75
legItems = []
for tag in impl:
    ax.boxplot(gen[tag], positions = np.arange(start, len(N)+start, 1),
               widths=.2,
               boxprops={'color' : colors[tag]},
               whiskerprops={'color' : colors[tag]},
               capprops={'color' : colors[tag]},
               flierprops={'color' : colors[tag]},
               medianprops={'color' : colors[tag]})
    legItems.append(mpatches.Patch(color=colors[tag], label=labels[tag]))
    start += .25

ax.set_yscale('log')

ax.axhline(60000000000, color='k')
ax.set_xlim(0, len(N)+1)
ax.set_xticklabels(N[1::2])
ax.set_xticks(np.arange(2, len(N)+1, 2))

ax.legend(handles=legItems, bbox_to_anchor=(1.1, 1.1))

closeFig("runtime/gen")

fig, ax = getFig("Interpolation Time per Point", r"$n$ sampling points", r"$t$ [ns]", None, "x")
start = .75
legItems = []
for tag in impl:
    ax.boxplot(ip[tag], positions = np.arange(start, len(N)+start, 1),
               widths=.2,
               boxprops={'color' : colors[tag]},
               whiskerprops={'color' : colors[tag]},
               capprops={'color' : colors[tag]},
               flierprops={'color' : colors[tag]},
               medianprops={'color' : colors[tag]})
    legItems.append(mpatches.Patch(color=colors[tag], label=labels[tag]))
    start += .25

ax.set_yscale('log')

ax.set_xlim(0, len(N)+1)
ax.set_xticklabels(N[1::2])
ax.set_xticks(np.arange(2, len(N)+1, 2))

ax.legend(handles=legItems, bbox_to_anchor=(1.1, 1.1))

closeFig("runtime/ip")

fig, ax = getFig(r"Runtime Triangulation + $10^4$ Interpolations", r"$n$ sampling points", r"$t$ [ns]", None, "x")
start = .75
legItems = []
for tag in impl:
    ax.boxplot(tt[tag], positions = np.arange(start, len(N)+start, 1),
               widths=.2,
               boxprops={'color' : colors[tag]},
               whiskerprops={'color' : colors[tag]},
               capprops={'color' : colors[tag]},
               flierprops={'color' : colors[tag]},
               medianprops={'color' : colors[tag]})
    legItems.append(mpatches.Patch(color=colors[tag], label=labels[tag]))
    start += .25

ax.set_yscale('log')

ax.set_xlim(0, len(N)+1)
ax.set_xticklabels(N[1::2])
ax.set_xticks(np.arange(2, len(N)+1, 2))

ax.legend(handles=legItems, bbox_to_anchor=(1.1, 1.1))

closeFig("runtime/total")

################################################################################
#precision figure

precision = {}

print("Reading in precision data")
for tag in impl:
    print("\t%s" % tag)
    precision[tag] = np.genfromtxt(os.path.join(dir, 'precision_%s.csv' % tag),
                                 names=True)

reference = precision['cgal']

print("Processing precision plots")
for f in F:

    print("\t%s" % funcNames[f])

    ipError = {}
    for tag in impl:
        ipError[tag] = []

        for n in N:
            err = []
            nums = precision[tag]
            nums = nums[(nums['n'] == n) & (nums['f'] == f)]
            for x in nums:
                err.append((x['ip'] - x['z']) / x['z'])
            ipError[tag].append(err)

    #boxplots of errors over N

    fig, ax = getFig(r"Interpolation Error - %s" % funcLabels[f], r"$n$ sampling points", r"Error - $\frac{\hat{z} - z}{z}$")

    start = .75
    legItems = []
    for tag in impl:
        ax.boxplot(ipError[tag], positions = np.arange(start, len(N)+start, 1),
                   widths=.2, whis=[5,95],
                   boxprops={'color' : colors[tag]},
                   whiskerprops={'color' : colors[tag]},
                   capprops={'color' : colors[tag]},
                   flierprops={'color' : colors[tag], 'marker' : '.'},
                   medianprops={'color' : colors[tag]})
        legItems.append(mpatches.Patch(color=colors[tag], label=labels[tag]))
        start += .25

    ax.axhline(0, color='k')

    ax.set_yscale('symlog', linthreshy=1e-3)
    ax.set_ylim(-100, 100)

    ax.set_xlim(0, len(N)+1)
    ax.set_xticklabels(N[1::2])
    ax.set_xticks(np.arange(2, len(N)+1, 2))

    ax.legend(handles=legItems, bbox_to_anchor=(1.1, 1.1))

    closeFig("precision/bp_error_%s" % funcNames[f])

    #histograms of individual error distributions
    for i,n in enumerate(N):
        fig, ax = getFig(r"Interpolation Error n = %i - %s" % (n, funcLabels[f]),
                r"Error - $\frac{\hat{z} - z}{z}$", "Frequency")

        legItems = []
        for tag in impl:
            data = ipError[tag][i]
            #select data in range [-100,100]
            data = [x for x in data if -100 < x and x < 100]

            hist, bin_edges = np.histogram(data,
                                           bins=100, range=(-100, 100), density=True)
            bin_centres = (bin_edges[:-1] + bin_edges[1:])/2

            (mu, sigma) = norm.fit(data)
            gauss = norm.pdf(bin_centres, loc=mu, scale=sigma)

            plt.plot(bin_centres, hist,
                     color=colors[tag], linestyle='None', marker=markers[tag])

#            plt.plot(bin_centres, gauss,
#                     color=colors[tag])

            legItems.append(mlines.Line2D([], [], color=colors[tag], marker=markers[tag],
                                          label=labels[tag]))

        ax.set_yscale('symlog', linthreshy=1e-4)
        ax.legend(handles=legItems)
        closeFig("precision/dist/%s_%i"%(funcNames[f] ,n))

    #3d plots of functions
#    for n in N:
#        fig, ax = getFig(r"Interpolation Quality n = %i - %s" % (n, funcLabels[f]), "x", "y", "z")
#
#        dNorm = reference[(reference['n'] == n) & (reference['f'] == f)]
#        ax.plot_wireframe(dNorm['x'], dNorm['y'], dNorm['z'], rstride=10, cstride=10,
#                          label='real', color='grey', alpha=0.4)
#
#        for tag in impl:
#            nums = precision[tag]
#            nums = nums[(nums['n'] == n) & (nums['f'] == f)]
#            ax.plot_wireframe(nums['x'],
#                              nums['y'],
#                              nums['ip'],
#                              rstride=10, cstride=10,
#                              label=labels[tag], color=colors[tag], alpha=0.4)
#        ax.legend(loc=2)
#        closeFig("precision/3d/%s_%i"%(funcNames[f] ,n))

print("done")
