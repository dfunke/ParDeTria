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

# define some aux data for labeling

class LabeledItem:

    def __init__(self, tag, label, color, marker):
        self.tag = tag
        self.label = label
        self.color = color
        self.marker = marker

splitters = [  LabeledItem(b'0', '0', 'b', 's')
             , LabeledItem(b'1', '1', 'g', '^')
             , LabeledItem(b'c', 'cycle', 'y', 'D')
             , LabeledItem(b'd', 'two-way', 'r', 'o')
            ]

print("Reading in data")
data = np.genfromtxt(os.path.join(dir, 'triangulation_report.csv'), dtype=None, names=True)

# validity
# check for invalid triangulations, only provenance '0' is of concern, as this is the entire triangulation
fig, ax = getFig("Validity", r"valid", r"")

ldata = data[data['provenance'] == b'0']
bins, _, _ = ax.hist(ldata['valid'], bins=2)

ax.text(0.1, 0.1, 'invalid = %i' %  bins[0], color='white')
ax.text(0.6, 0.1, 'valid = %i' % bins[1], color='white')

closeFig("report/validity")

# plot edge size
# only non base case triangulations are relevant
fig, ax = getFig("Edge Size", r"$n$", r"edge points")

for s in splitters:
    ldata = data[(data['splitter'] == s.tag) & (data['base_case'] == 0)]
    ax.plot(ldata['nPoints'], ldata['nEdgePoints'], label=s.label, color=s.color, marker=s.marker, ls='None')

ax.legend()

closeFig("report/edgePoints")

# plot edge size
fig, ax = getFig("Edge Size", r"$n$", r"edge simplices")

for s in splitters:
    ldata = data[(data['splitter'] == s.tag) & (data['base_case'] == 0)]
    ax.plot(ldata['nPoints'], ldata['nEdgeSimplices'], label=s.label, color=s.color, marker=s.marker, ls='None')

ax.legend()

closeFig("report/edgeSimplices")

# plot runtime
fig, ax = getFig("Runtime", r"$n$", r"$t$")

for s in splitters:
    ldata = data[data['splitter'] == s.tag]
    ax.plot(ldata['n'], ldata['time'], label=s.label, color=s.color, marker=s.marker, ls='None')

ax.legend()

closeFig("report/runtime")

# plot runtime scaled
fig, ax = getFig("Runtime", r"$n$", r"$\frac{t}{n \log n}$")

for s in splitters:
    ldata = data[data['splitter'] == s.tag]
    ax.plot(ldata['n'], ldata['time'] / (ldata['n'] * np.log2(ldata['n'])), label=s.label, color=s.color, marker=s.marker, ls='None')

ax.legend()

closeFig("report/runtime_scaled")

# plot mem
fig, ax = getFig("Memory", r"$n$", r"$mem$")

for s in splitters:
    ldata = data[data['splitter'] == s.tag]
    ax.plot(ldata['n'], ldata['mem'], label=s.label, color=s.color, marker=s.marker, ls='None')

ax.legend()

closeFig("report/memory")

print("Done")
