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

print("Reading in data")
edge2 = np.genfromtxt(os.path.join(dir, 'edge_study_2.csv'), names=True)
edge2.sort(order='nP')

edge3 = np.genfromtxt(os.path.join(dir, 'edge_study_3.csv'), names=True)
edge3.sort(order='nP')

partition2d= np.genfromtxt(os.path.join(dir, 'partition_study_2_d.csv'), names=True)
partition2d.sort(order='nP')

partition2c= np.genfromtxt(os.path.join(dir, 'partition_study_2_c.csv'), names=True)
partition2c.sort(order='nP')

partition3d= np.genfromtxt(os.path.join(dir, 'partition_study_3_d.csv'), names=True)
partition3d.sort(order='nP')

partition3c= np.genfromtxt(os.path.join(dir, 'partition_study_3_c.csv'), names=True)
partition3c.sort(order='nP')

print("Generating plots")

fig, ax = getFig("Number of Simplices Over Number of Points", r"$p$", r"s")

ax.plot(edge2['nP'], edge2['eS'], label="2D")
ax.plot(edge3['nP'], edge3['eS'], label="3D")
ax.plot(edge2['nP'], edge2['nP'], ls='--', color='k')

ax.legend()

closeFig("edge_study/01_simplices")

fig, ax = getFig("Edge Points Over Number of Points", r"$p$", r"edge points")

ax.plot(edge2['nP'], edge2['eP'], label="2D")
ax.plot(edge3['nP'], edge3['eP'], label="3D")
ax.plot(edge2['nP'], edge2['nP'], ls='--', color='k')

ax.legend()

closeFig("edge_study/02_edgePoints")

fig, ax = getFig("Percentage of Edge Points over Points", r"$p$", r"%")

ax.plot(edge2['nP'], edge2['eP'] / edge2['nP'], label="2D")
ax.plot(edge3['nP'], edge3['eP'] / edge3['nP'], label="3D")
ax.axhline(.5, ls='--', color='black')

ax.legend()

closeFig("edge_study/03_edgePoints_perc")


fig, ax = getFig("Edge Simplices of Number of Simplices", r"$s$", r"edge simplices")

ax.plot(edge2['nS'], edge2['eS'], label="2D")
ax.plot(edge3['nS'], edge3['eS'], label="3D")
ax.plot(edge3['nS'], edge3['nS'], ls='--', color='k')

ax.legend()

closeFig("edge_study/04_edgeSimplices")

fig, ax = getFig("Percentage of Edge Simplices of Simplices", r"$s$", r"%")

ax.plot(edge2['nS'], edge2['eS'] / edge2['nS'], label="2D")
ax.plot(edge3['nS'], edge3['eS'] / edge3['nS'], label="3D")
ax.axhline(.5, ls='--', color='black')

ax.legend()

closeFig("edge_study/05_edgeSimplices_perc")

fig, ax = getFig("Edge Points Over Number of Points", r"$p$", r"edge points")

ax.plot(partition2d['nP'], partition2d['eP'], label="2D - d part")
ax.plot(partition2c['nP'], partition2c['eP'], label="2D - c part")
ax.plot(partition3d['nP'], partition3d['eP'], label="3D - d part")
ax.plot(partition3c['nP'], partition3c['eP'], label="3D - c part")
ax.plot(partition2d['nP'], partition2d['nP'], ls='--', color='k')

ax.legend()

closeFig("edge_study/06_partition_edgePoints")

fig, ax = getFig("Percentage of Edge Points Over Points", r"$p$", r"%")

ax.plot(partition2d['nP'], partition2d['eP'] / partition2d['nP'], label="2D - d part")
ax.plot(partition2c['nP'], partition2c['eP'] / partition2c['nP'], label="2D - c part")
ax.plot(partition3d['nP'], partition3d['eP'] / partition3d['nP'], label="3D - d part")
ax.plot(partition3c['nP'], partition3c['eP'] / partition3c['nP'], label="3D - c part")
ax.axhline(.5, ls='--', color='black')

ax.legend()

closeFig("edge_study/07_partition_edgePoints_perc")

print("Done")
