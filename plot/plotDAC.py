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

# data format
# n splitter provenance base_case edge_tria valid nPoints nSimplices nEdgePoints nEdgeSimplices time mem max_mem


print("Reading in data")
data = np.genfromtxt(os.path.join(dir, 'triangulation_report.csv'), dtype=None, names=True)

# validity

fig, ax = getFig("Validity", r"valid", r"")

ax.hist(data['valid'])

closeFig("report/validity")

#setup unique data
splitters = np.unique(data['splitter'])

# plot edge size
fig, ax = getFig("Edge Size", r"$n$", r"edge points")

for s in splitters:
    ldata = data[data['splitter'] == s]
    ax.plot(ldata['nPoints'], ldata['nEdgePoints'], 's', label=s)

ax.legend()

closeFig("report/edgePoints")

# plot edge size
fig, ax = getFig("Edge Size", r"$n$", r"edge simplices")

for s in splitters:
    ldata = data[data['splitter'] == s]
    ax.plot(ldata['nPoints'], ldata['nEdgeSimplices'], 's', label=s)

ax.legend()

closeFig("report/edgeSimplices")
