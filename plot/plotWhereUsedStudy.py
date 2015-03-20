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

dir = 'data'

print("Reading in data")
whereUsedSet2 = np.genfromtxt(os.path.join(dir, 'whereUsed_study_set_2_c.csv'), names=True)
whereUsedSet2.sort(order='nP')

whereUsedArray2 = np.genfromtxt(os.path.join(dir, 'whereUsed_study_array_2_c.csv'), names=True)
whereUsedArray2.sort(order='nP')

whereUsedSet3 = np.genfromtxt(os.path.join(dir, 'whereUsed_study_set_3_c.csv'), names=True)
whereUsedSet3.sort(order='nP')

whereUsedArray3 = np.genfromtxt(os.path.join(dir, 'whereUsed_study_array_3_c.csv'), names=True)
whereUsedArray3.sort(order='nP')

N2D = np.unique(whereUsedSet2['nP'])
N3D = np.unique(whereUsedSet3['nP'])

print("Generating plots")

for n in N2D:

    fig, ax = getFig("Length of WhereUsed List Over Number of Points", r"$p$", r"$l$")

    plt.hist(whereUsedSet2[whereUsedSet2['nP'] == n]['nWU'], label='Set')
    plt.hist(whereUsedArray2[whereUsedArray2['nP'] == n]['nWU'], label='Array - active')
    plt.hist(whereUsedArray2[whereUsedArray2['nP'] == n]['nIWU'], label='Array - deleted')

    ax.legend()

    closeFig("whereUsed_study/01_length2D_%i" % n)

for n in N3D:

#Lenght of array/set
    fig, ax = getFig("Length of WhereUsed List Over Number of Points", r"$p$", r"$l$")

    plt.hist([whereUsedArray3[whereUsedArray3['nP'] == n]['nWU'] - whereUsedArray3[whereUsedArray3['nP'] == n]['nIWU'],
              whereUsedArray3[whereUsedArray3['nP'] == n]['nIWU']],
             label=['Array - active', 'Array - deleted'], align='right', histtype='barstacked', log=True, bins=100)
    plt.hist(whereUsedSet3[whereUsedSet3['nP'] == n]['nWU'], label = 'Set', align='left', log=True, bins=20)


    ax.legend()

    closeFig("whereUsed_study/01_length3D_%i" % n)

#Percent of deleted items

    fig, ax = getFig("Percentage of Deleted in WhereUsed List Over Number of Points", r"$p$", r"%")

    plt.hist(whereUsedArray3[whereUsedArray3['nP'] == n]['nIWU'] / whereUsedArray3[whereUsedArray3['nP'] == n]['nWU'],
             label='Array - % deleted')


    ax.legend()

    closeFig("whereUsed_study/02_deleted3D_%i" % n)

print("Done")
