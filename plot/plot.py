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
        file = 'plot/output/' + file
        
    if not '.' in file:
        file = file + '.png'
    
    ensure_dir(file)
    
    plt.savefig(file)
    
    plt.close()
    
dir = 'data'
files = [ f for f in os.listdir(dir) if '.csv' in f ]

print("processing files: ", files)

data = {}
for f in files:
    data[f.replace('.csv','')] = np.genfromtxt(os.path.join(dir, f), names=True)
    
#determine order at last data point
tags = sorted(data.keys(), key = lambda k: np.max(data[k]['time']), reverse=True)

#runtime figure
fig, ax = getFig("Runtime over Number of Points", "n", "t [ns]")

for tag in tags:
    ax.errorbar(data[tag]['n'], data[tag]['time'], yerr=data[tag]['time_std'], label=tag)

ax.legend(loc=2)

closeFig("runtime")

#runtime figure
fig, ax = getFig("Speedup over ROOT with Number of Points", "n", "t [ns]")

baseline = data['root']['time']
for tag in reversed(tags): #reverse tags, as speedup is opposite to runtime
    ax.plot(data[tag]['n'], baseline / data[tag]['time'], label=tag)

ax.legend(loc=2)

closeFig("speedup")

print("done")