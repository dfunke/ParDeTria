import numpy as np
import pandas as pd
import matplotlib as plt

class Series:

    def __init__(self):
        self.xvalues = []
        self.yvalues = []
        self.color = None
        self.marker = None
        self.linestyle = None
        self.label = None

class Plot:

    def __init__(self):
        self.series = []
        self.title = None
        self.xlabel = None
        self.ylabel = None

    def _ensure_dir(f : str):
        import os
        d = os.path.dirname(f)
        if not os.path.exists(d):
            os.makedirs(d)

    def _sanitize_file(f: str) -> str:
        if not '/' in f:
            f = 'plot/output/' + f

        if not '.' in f:
            f = f + '.png'

        return f

    def plot(self, filename : str):

        fig, ax = plt.subplots()

        fig.set_title(self.title)
        ax.set_xlabel(self.xlabel)
        ax.set_ylabel(self.ylabel)

        for s in self.series:
            ax.plot(s.xvalues, s.yvalues, label=s.label, marker=s.marker, color=s.color, linestyle=s.linestyle)

        ax.legend()

        filename = Plot._sanitize_file(filename)
        Plot._ensure_dir(filename)

        plt.savefig(filename)
        plt.close()