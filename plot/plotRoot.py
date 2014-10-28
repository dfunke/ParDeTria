import numpy as np
import matplotlib.pylab as plt
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

for tag in impl:
    precision[tag] = np.genfromtxt(os.path.join(dir, 'precision_%s.csv' % tag),
                                 names=True)

norm = precision['cgal']
N = np.unique(norm['n'])

ipError = {}
for tag in impl:
    ipError[tag] = []

    for n in N:
        chi2 = []
        nums = precision[tag]
        nums = nums[nums['n'] == n]
        for x in nums:
            chi2.append((x['ip'] - x['z'])**2)
        ipError[tag].append(chi2)

fig, ax = getFig("Mean Interpolation Error", "n", "Error")

for tag in impl:
    ax.errorbar(N, np.mean(ipError[tag], axis=1), yerr=np.std(ipError[tag], axis=1), label=tag)

ax.set_xscale('log')
ax.set_yscale('log')

ax.legend(loc=2)
closeFig("interpolation_error")
    #fig, ax = getFig("Interpolation Quality n = %i" % n, "x", "y", "z")

    #ax.plot_wireframe(norm['x'], norm['y'], norm['z'], label="real")

    #for tag in impl:
    #    ax.plot_wireframe(precision[tag]['x'],
    #                      precision[tag]['y'],
    #                      precision[tag]['ip'],
    #                      label=tag)
    #ax.legend(loc=2)
    #closeFig("precision_%i"%n)

print("done")
