__author__ = 'dfunke'

import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

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
        kwargs = {}
        if not self.label is None: kwargs['label'] = self.label
        if not self.marker is None: kwargs['marker'] = self.marker
        if not self.color is None: kwargs['color'] = self.color
        if not self.linestyle is None: kwargs['linestyle'] = self.linestyle

        if self. xvalues is None:
            ax.axhline(self.yvalue, **kwargs)
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

    def _plot(self, ax : plt.Axes):
        for s in self.series:
            s.plot(ax)

    def _legend(self, obj, **kwargs):
        if not 'loc' in kwargs:
            kwargs['loc'] = 0

        obj.legend(**kwargs)

    def plot(self, filename : str,
             figureArgs : dict = {}, axesArgs : dict = {},
             legendArgs : dict = {}, figureLegend : bool = False):

        fig = plt.figure(**figureArgs)
        ax = fig.add_subplot(111, **axesArgs)

        ax.set_title(self.title)
        ax.set_xlabel(self.xlabel)
        ax.set_ylabel(self.ylabel)

        self._plot(ax)

        if not self.desc is None:
            fig.text(0.12, 0.9, self.desc, va='bottom', ha='left')

        if figureLegend:
            self._legend(fig, **legendArgs)
        else:
            self._legend(ax, **legendArgs)

        filename = Plot._sanitize_file(filename)
        Plot._ensure_dir(filename)

        plt.savefig(filename)
        plt.close()

class StackPlot(Plot):

    def _legend(self, obj, **kwargs):

        handles = []
        labels = []

        for s, c in zip(self.series, self.colors):
            handles.append(plt.Rectangle((0, 0), 1, 1, fc=c))
            labels.append(s.label)

        obj.legend(handles=handles[::-1], labels=labels[::-1], **kwargs)

    def _plot(self, ax : plt.Axes):

        if len(self.series) < 1:
            # empty plot
            return

        x = self.series[0].xvalues
        y = []

        for s in self.series:
            if (x != s.xvalues).any():
                # for the stack plot to make sense,
                # the x values need to be equal for all series
                print("Unequal x values for stack plot")
                return

            y.append(s.yvalues)

        self.colors = plt.get_cmap('jet')(np.linspace(0, 1.0, len(self.series)))
        ax.stackplot(x, y, colors=self.colors)