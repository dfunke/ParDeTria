import numpy as np
import matplotlib.pylab as plt
import matplotlib.patches as mpatches
from matplotlib.ticker import FormatStrFormatter
from mpl_toolkits.mplot3d import Axes3D

import os

def ensure_dir(f):
    d = os.path.dirname(f)
    if not os.path.exists(d):
        os.makedirs(d)
        
def getFig(title, xlabel, ylabel, zlabel = None):
    if not zlabel:
        fig, ax = plt.subplots()
        ax.grid(b=True, which='both')
    else:
        fig, ax = plt.subplots(subplot_kw={'projection' :'3d'})
        ax.set_zlabel(zlabel)

    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    
    return fig, ax


def closeFig(file):
    
    if not '/' in file:
        file = 'plot/output/' + file
        
    if not '.' in file:
        file = file + '.png'
    
    ensure_dir(file)
    
    plt.savefig(file)
    
    plt.close()
    
dir = 'data'
impl = ['cgal', 'tria']

funcNames = [  'sin'
              ,'sinsin'
              ,'cubic'
              ,'sqrtsincos']

funcLabels = [ r'$\sin(4 xy)$'
              ,r'$ \frac{\sin{x}}{x} \cdot \frac{\sin{y}}{y} + 0.2$'
              ,r'$x^3 - 3x + y^3 -3y$'
              ,r'$\sqrt{\sin{x}+\cos{y}}$']


timings = {}

for tag in impl:
    timings[tag] = np.genfromtxt(os.path.join(dir, 'timings_%s.csv' % tag),
                                 names=True)

#gen time figure
fig, ax = getFig("Generation Time over Number of Points", "n", "t [ns]")

for tag in impl:
    ax.errorbar(timings[tag]['n'],
                timings[tag]['gen'],
                yerr=timings[tag]['gen_std'],
                label=tag)

ax.set_xscale('log')
ax.set_yscale('log')

ax.legend(loc=2)

closeFig("runtime_root_gen")

#gen time per point figure
fig, ax = getFig("Generation Time per Point", "n", "t/n [ns]")

for tag in impl:
    ax.errorbar(timings[tag]['n'],
            timings[tag]['gen'] / timings[tag]['n'],
            yerr=timings[tag]['gen_std'] / timings[tag]['n'],
            label=tag)

ax.set_xscale('log')
ax.set_yscale('log')

ax.legend(loc=2)

closeFig("runtime_root_gen_pp")

#interpolation over points figure
fig, ax = getFig("Interpolation Time over Number of Points", "n", "t [ns]")

for tag in impl:
    ax.errorbar(timings[tag]['n'],
            timings[tag]['ip'],
            yerr=timings[tag]['ip_std'],
            label=tag)

ax.set_xscale('log')
ax.set_yscale('log')

ax.legend(loc=2)

closeFig("runtime_root_ip")

#precision figure

precision = {}
colors = {'cgal' : 'blue', 'tria' : 'green'}

for tag in impl:
    precision[tag] = np.genfromtxt(os.path.join(dir, 'precision_%s.csv' % tag),
                                 names=True)

norm = precision['cgal']
N = np.unique(norm['n']).astype(int)
F = np.unique(norm['f']).astype(int)

for f in F:
    ipError = {}
    for tag in impl:
        ipError[tag] = []

        for n in N:
            chi2 = []
            nums = precision[tag]
            nums = nums[(nums['n'] == n) & (nums['f'] == f)]
            for x in nums:
                chi2.append((x['ip'] - x['z'])**2)
            ipError[tag].append(chi2)

    #boxplots

    fig, ax = getFig(r"Interpolation Error - %s" % funcLabels[f], "n", "Error")

    start = .75
    legItems = []
    for tag in impl:
        ax.boxplot(ipError[tag], positions = np.arange(start, len(N)+start, 1),
                   widths=.25,
                   boxprops={'color' : colors[tag]},
                   whiskerprops={'color' : colors[tag]},
                   capprops={'color' : colors[tag]},
                   flierprops={'color' : colors[tag]})
        legItems.append(mpatches.Patch(color=colors[tag], label=tag))
        start += .5

    ax.set_yscale('log')

    ax.set_xlim(0, len(N)+1)
    #ax.set_xticklabels(["%.0f" % x for x in N[0::2]])
    ax.set_xticklabels(N[0::2])
    ax.set_xticks(np.arange(1, len(N), 2))

    ax.legend(handles=legItems)

    closeFig("interpolation_error_boxplot_%s" % funcNames[f])


#    #mean interpolation error
#
#    fig, ax = getFig("Mean Interpolation Error", "n", "Error")
#
#    for tag in impl:
#        ax.errorbar(N, np.mean(ipError[tag], axis=1), yerr=np.std(ipError[tag], axis=1), label=tag)
#
#    ax.set_xscale('log')
#    ax.set_yscale('log')
#
#    ax.legend(loc=2)
#    closeFig("interpolation_error")
#
    for n in N:
        fig, ax = getFig(r"Interpolation Quality n = %i - %s" % (n, funcLabels[f]), "x", "y", "z")

        dNorm = norm[(norm['n'] == n) & (norm['f'] == f)]
        ax.plot_wireframe(dNorm['x'], dNorm['y'], dNorm['z'], rstride=10, cstride=10,
                          label='real', color='grey', alpha=0.4)

        for tag in impl:
            nums = precision[tag]
            nums = nums[(nums['n'] == n) & (nums['f'] == f)]
            ax.plot_wireframe(nums['x'],
                              nums['y'],
                              nums['ip'],
                              rstride=10, cstride=10, 
                              label=tag, color=colors[tag], alpha=0.4)
        ax.legend(loc=2)
        closeFig("precision_%s_%i"%(funcNames[f] ,n))

print("done")
