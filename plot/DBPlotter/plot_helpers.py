import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import copy

class Series:

    def __init__(self, o = None):
        if o is None:
            self.xvalues = []
            self.yvalues = []
            self.color = None
            self.marker = None
            self.linestyle = None
            self.label = None
        else:
            self.xvalues = copy.copy(o.xvalues)
            self.yvalues = copy.copy(o.yvalues)

            self.color = o.color
            self.marker = o.marker
            self.linestyle = o.linestyle
            self.label = o.label
            
    def plot(self, ax : plt.Axes):
        kargs = {}
        if not self.label is None: kargs['label'] = self.label
        if not self.marker is None: kargs['marker'] = self.marker
        if not self.color is None: kargs['color'] = self.color
        if not self.linestyle is None: kargs['linestyle'] = self.linestyle

        ax.plot(self.xvalues, self.yvalues, **kargs)

class HLine(Series):

    def __init__(self, o = None):
        if o is None:
            self.xvalues = None
            self.yvalue = None
            self.color = None
            self.marker = None
            self.linestyle = None
            self.label = None
        else:
            if(not o.xvalues is None):
                self.xvalues = copy.copy(o.xvalues)
            else:
                self.xvalues = None
            self.yvalue = o.yvalue

            self.color = o.color
            self.marker = o.marker
            self.linestyle = o.linestyle
            self.label = o.label

    def plot(self, ax : plt.Axes):
        kargs = {}
        if not self.label is None: kargs['label'] = self.label
        if not self.marker is None: kargs['marker'] = self.marker
        if not self.color is None: kargs['color'] = self.color
        if not self.linestyle is None: kargs['linestyle'] = self.linestyle

        if self. xvalues is None:
            ax.axhline(self.yvalue, **kargs)
        else:
            ax.plot(self.xvalues, np.full_like(self.xvalues, self.yvalue), **kargs)

class Plot:

    def __init__(self):
        self.series = []
        self.title = None
        self.xlabel = None
        self.ylabel = None
        self.desc = None

    def _ensure_dir(f : str):
        import os
        d = os.path.dirname(f)
        if not os.path.exists(d):
            os.makedirs(d)

    def _sanitize_file(f: str) -> str:
        if not f.startswith('plot/output/dbdump/'):
            f = 'plot/output/dbdump/' + f

        if not '.' in f:
            f = f + '.png'

        return f

    def addSeries(self, series : Series):
        self.series.append(series)

    def plot(self, filename : str, legend_cols : int = 1, legend_loc : int = 0):

        fig, ax = plt.subplots()

        ax.set_title(self.title)
        ax.set_xlabel(self.xlabel)
        ax.set_ylabel(self.ylabel)

        for s in self.series:
            s.plot(ax)

        if not self.desc is None:
            fig.text(0.12, 0.9, self.desc, va='bottom', ha='left')

        ax.legend(loc=legend_loc, ncol=legend_cols)

        filename = Plot._sanitize_file(filename)
        Plot._ensure_dir(filename)

        plt.savefig(filename)
        plt.close()