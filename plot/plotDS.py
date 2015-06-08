import numpy as np
import matplotlib.pylab as plt

import os

def ensure_dir(f):
    d = os.path.dirname(f)
    if not os.path.exists(d):
        os.makedirs(d)
        
def getFig(title, xlabel, ylabel):
    fig, ax = plt.subplots()
    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    
    ax.grid(b=True, which='both')
    
    return fig, ax


def closeFig(file):
    
    if not '/' in file:
        file = 'plot/output/ds/' + file
        
    if not '.' in file:
        file = file + '.png'
    
    ensure_dir(file)
    
    plt.savefig(file)
    
    plt.close()

data = np.genfromtxt("data/benchmark_set.csv", names=True)

#runtime figure
fig, ax = getFig("Runtime", "threads", "speedup")

ax.plot(data['threads'], data['tbb'], label="TBB")
ax.plot(data['threads'], data['own'], label="own")

ax.legend(loc=1)

closeFig("time")

#runtime figure
fig, ax = getFig("Relative Speedup", "threads", "speedup")

ax.plot(data['threads'], data['tbb'][0] / data['tbb'], label="TBB")
ax.plot(data['threads'], data['own'][0] / data['own'], label="own")

ax.legend(loc=2)

closeFig("rel")

#runtime figure
fig, ax = getFig("Absolute Speedup", "threads", "speedup")

ax.plot(data['threads'], data['tbb'][0] / data['tbb'], label="TBB")
ax.plot(data['threads'], data['tbb'][0] / data['own'], label="own")

ax.legend(loc=2)

closeFig("abs")

print("done")
