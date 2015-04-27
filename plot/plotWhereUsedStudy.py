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

print("Reading in data")

rawData = np.genfromtxt(os.path.join(dir, 'whereUsed_collisions_3_c.csv'), names=True, dtype=int)
rawData.sort(order='n')

N = np.unique(rawData['n'])

print("Processing data")
colData = []
for n in N:
    colData.append(rawData[rawData['n'] == n]['collisions'])

print("Generating plots")

fig, ax = getFig("Collisions over Points", r"$n$", r"$c$")

ax.boxplot(colData, whis=[5,95])

ax.set_xlim(0, len(N)+1)
ax.set_xticklabels((N-8)[0::9])
ax.set_xticks(np.arange(1, len(N)+1, 9))

#ax.legend()

closeFig("whereUsed_study/04_collisions")


print("Done")
