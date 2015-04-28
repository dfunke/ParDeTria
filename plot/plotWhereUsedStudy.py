#! /usr/bin/env python3

import numpy as np
import scipy as sp

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

from matplotlib import rcParams
rcParams['xtick.direction'] = 'out'
rcParams['ytick.direction'] = 'out'

dir = 'data'
files = [ f for f in os.listdir(dir)
        if f.startswith('whereUsed_collisions_3_c_') and f.endswith('.csv')]

impl = ['rol', 'xor', 'xorSeed']
labels =  {'rol' : 'rotate', 'xor' : 'xor', 'xorSeed' : 'xor with seed'}
colors  = {'rol' : 'blue', 'xor' : 'green', 'xorSeed' : 'red'}
markers = {'rol' : 'x', 'xor' : 'D', 'xorSeed' : 'o'}

colData = {}
for f in files:
    print("Reading in data from %s" % f)
    alg = f.replace('whereUsed_collisions_3_c_', '').replace('.csv','')

    rawData = np.genfromtxt(os.path.join(dir, f), names=True, dtype=int)
    rawData.sort(order='n')

    N = np.unique(rawData['n'])

    print("Processing data")
    colData[alg] = []
    for n in N:
        colData[alg].append(rawData[rawData['n'] == n]['collisions'])

    print("Generating plots")

    fig, ax = getFig("Collisions over Points", r"$n$", r"$c$")

    ax.boxplot(colData[alg], whis=[5,95],
               boxprops={'color' : colors[alg]},
               whiskerprops={'color' : colors[alg]},
               capprops={'color' : colors[alg]},
               flierprops={'color' : colors[alg], 'marker' : markers[alg], 'markersize' : .9},
               medianprops={'color' : 'magenta'})

    ax.set_xlim(0, len(N)+1)
    ax.set_xticklabels((N-8)[0::9])
    ax.set_xticks(np.arange(1, len(N)+1, 9))

    _, ymax = ax.get_ylim()
    ymax += 1
    ax.set_ylim(0, ymax)

    #ax.legend()
    fig.text(0.12, 0.9, labels[alg], va='bottom', ha='left')

    closeFig("whereUsed_study/00_collisions_%s" % alg)

print("Generating combined plots")

fig, ax = getFig("Collisions over Points", r"$n$", r"$c$")

start = .75
legItems = []
for alg in colData:
    ax.boxplot(colData[alg], positions = np.arange(start, len(N)+start, 1),
               widths=.2, whis=[5,95],
               boxprops={'color' : colors[alg]},
               whiskerprops={'color' : colors[alg]},
               capprops={'color' : colors[alg]},
               flierprops={'color' : colors[alg], 'marker' : markers[alg], 'markersize' : .5},
               medianprops={'color' : 'magenta'})
    legItems.append(mpatches.Patch(color=colors[alg], label=labels[alg]))
    start += .25

ax.set_xlim(0, len(N)+1)
ax.set_xticklabels((N-8)[0::9])
ax.set_xticks(np.arange(1, len(N)+1, 9))

_, ymax = ax.get_ylim()
ymax += 1
ax.set_ylim(0, ymax)

ax.legend(handles=legItems, bbox_to_anchor=(1.1, 1.1))

closeFig("whereUsed_study/00_collisions_combined")


print("Done")
