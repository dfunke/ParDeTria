#! /usr/bin/env python3

__author__ = 'dfunke'

import sys
import pyejdb
import numpy as np
import pandas as pd
import math
import copy

import plot_helpers as ph
import db_helpers as dh

def plotCGAL():
    lastRunCGAL = dh.getLastRun(database, 'pureCGAL')
    maxPoints = dh.selectMax(database, 'pureCGAL', 'nP', {'run-number' : lastRunCGAL})

    pureCGAL = dh.load(database, 'pureCGAL', {'run-number' : lastRunCGAL, 'nP' : maxPoints})
    charsCGAL = dh.getCharacteristics(pureCGAL)

    pRuntime = ph.Plot()
    pRuntime.title = "CGAL MT: Runtime over Threads"
    pRuntime.xlabel = "threads"
    pRuntime.ylabel = "time [ms]"
    pRuntime.desc = r"$10^%i$ points" % math.log10(charsCGAL['nP'][0])

    pSpeedup = ph.Plot()
    pSpeedup.title = "CGAL MT: Speedup over Threads"
    pSpeedup.xlabel = "threads"
    pSpeedup.ylabel = "speedup"
    pSpeedup.desc = r"$10^%i$ points" % math.log10(charsCGAL['nP'][0])

    for occ in charsCGAL['occupancy']:
        sRuntime = ph.Series()
        sRuntime.xvalues = charsCGAL['threads']
        sRuntime.yvalues = [np.mean(x) / 1e6 for x in dh.select(pureCGAL, {'occupancy' : occ}).sort(columns='threads')['times']]

        sRuntime.label = "CGAL - Lock Granularity %s" % occ

        pRuntime.addSeries(sRuntime)

        sSpeedup = ph.Series(sRuntime)
        sSpeedup.yvalues = sSpeedup.yvalues[0] / sSpeedup.yvalues

        pSpeedup.addSeries(sSpeedup)

    pRuntime.plot("pureCGAL_time_over_threads.png")
    pSpeedup.plot("pureCGAL_speedup_over_threads.png")

def plotBaseCase():
    lastRunBenchmarks = dh.getLastRun(database, 'benchmarks')
    maxPoints = dh.selectMax(database, 'benchmarks', 'nP', {'run-number' : lastRunBenchmarks})

    benchmarks = dh.load(database, 'benchmarks', {'run-number' : lastRunBenchmarks, 'nP' : maxPoints})
    charsBenchmarks = dh.getCharacteristics(benchmarks)

    pRuntime = ph.Plot()
    pRuntime.title = "Runtime over Threads"
    pRuntime.xlabel = "threads"
    pRuntime.ylabel = "time [ms]"
    pRuntime.desc = r"$10^%i$ points" % math.log10(charsBenchmarks['nP'][0])

    pSpeedup = ph.Plot()
    pSpeedup.title = "Speedup over Threads"
    pSpeedup.xlabel = "threads"
    pSpeedup.ylabel = "speedup"
    pSpeedup.desc = r"$10^%i$ points" % math.log10(charsBenchmarks['nP'][0])

    for bs in charsBenchmarks['basecase']:
        # sequential base solver
        sRSeq = ph.Series()
        sRSeq.xvalues = charsBenchmarks['threads']
        sRSeq.yvalues = [np.mean(x) / 1e6 for x in dh.select(benchmarks, {'basecase' : bs, 'parallel-base': False}).sort(columns='threads')['times']]

        sRSeq.label = "seq. Base %i" % bs

        pRuntime.addSeries(sRSeq)

        sSSeq = ph.Series(sRSeq)
        sSSeq.yvalues = sSSeq.yvalues[0] / sSSeq.yvalues

        pSpeedup.addSeries(sSSeq)
        
    for bs in charsBenchmarks['basecase']:
        # parallel base solver
        
        sRPar = ph.Series()
        sRPar.xvalues = charsBenchmarks['threads']
        sRPar.yvalues = [np.mean(x) / 1e6 for x in dh.select(benchmarks, {'basecase' : bs, 'parallel-base': True}).sort(columns='threads')['times']]

        sRPar.label = "par. Base %i" % bs

        pRuntime.addSeries(sRPar)

        sSPar = ph.Series(sRPar)
        sSPar.yvalues = sSPar.yvalues[0] / sSPar.yvalues

        pSpeedup.addSeries(sSPar)

    pRuntime.plot("basecase_time_over_threads.png", {'ncol' : 2})
    pSpeedup.plot("basecase_speedup_over_threads.png", {'ncol' : 2, 'loc' : 2})

def plotComparison():
    ########################################################################################################################

    lastRunBenchmarks = dh.getLastRun(database, 'benchmarks')
    maxPoints = dh.selectMax(database, 'benchmarks', 'nP', {'run-number' : lastRunBenchmarks})

    benchmarks = dh.load(database, 'benchmarks', {'run-number' : lastRunBenchmarks, 'nP' : maxPoints})

    charsBenchmarks = dh.getCharacteristics(benchmarks)

    pRuntime = ph.Plot()
    pRuntime.title = "Runtime over Threads"
    pRuntime.xlabel = "threads"
    pRuntime.ylabel = "time [ms]"
    pRuntime.desc = r"$10^%i$ points" % math.log10(charsBenchmarks['nP'][0])

    pSpeedupRel = ph.Plot()
    pSpeedupRel.title = "Relative Speedup over Threads"
    pSpeedupRel.xlabel = "threads"
    pSpeedupRel.ylabel = "speedup"
    pSpeedupRel.desc = r"$10^%i$ points" % math.log10(charsBenchmarks['nP'][0])

    pSpeedupAbs = ph.Plot()
    pSpeedupAbs.title = "Absolute Speedup over Threads"
    pSpeedupAbs.xlabel = "threads"
    pSpeedupAbs.ylabel = "speedup"
    pSpeedupAbs.desc = r"$10^%i$ points" % math.log10(charsBenchmarks['nP'][0])

    #################################################
    # CGAL sequential

    sCGAL = ph.HLine()
    sCGAL.xvalues = charsBenchmarks['threads']
    sCGAL.yvalue = np.mean(list(dh.select(benchmarks, {'alg' : 'c'}).sort(columns='threads')['times'])) / 1e6
    sCGAL.label = "CGAL sequential"
    pRuntime.addSeries(sCGAL)

    sSCGAL = ph.HLine(sCGAL)
    sSCGAL.yvalue = 1
    pSpeedupRel.addSeries(sSCGAL)
    pSpeedupAbs.addSeries(sSCGAL)

    #################################################
    # CGAL MT

    sCGALMT = ph.Series()
    sCGALMT.xvalues = charsBenchmarks['threads']
    sCGALMT.yvalues = [np.mean(x) / 1e6 for x in dh.select(benchmarks, {'alg' : 'm'}).sort(columns='threads')['times']]
    sCGALMT.label = "CGAL MT"
    pRuntime.addSeries(sCGALMT)

    sSCGALMT = ph.Series(sCGALMT)
    sSCGALMT.yvalues = sCGALMT.yvalues[0] / sCGALMT.yvalues
    pSpeedupRel.addSeries(sSCGALMT)

    sASCGALMT = ph.Series(sCGALMT)
    sASCGALMT.yvalues = sCGAL.yvalue / sCGALMT.yvalues
    pSpeedupAbs.addSeries(sASCGALMT)

    #################################################
    # small base case, seq. base

    sSeqBase = ph.Series()
    sSeqBase.xvalues = charsBenchmarks['threads']
    sSeqBase.yvalues = [np.mean(x) / 1e6 for x in dh.select(benchmarks, {'alg': 'd', 'parallel-base' : False, 'basecase' : np.min(charsBenchmarks['basecase'])}).sort(columns='threads')['times']]
    sSeqBase.label = "Seq. Base %i" % np.min(charsBenchmarks['basecase'])
    pRuntime.addSeries(sSeqBase)

    sSSeqBase = ph.Series(sSeqBase)
    sSSeqBase.yvalues = sSeqBase.yvalues[0] / sSeqBase.yvalues
    pSpeedupRel.addSeries(sSSeqBase)

    sASSeqBase = ph.Series(sSeqBase)
    sASSeqBase.yvalues = sCGAL.yvalue / sSeqBase.yvalues
    pSpeedupAbs.addSeries(sASSeqBase)

    #################################################
    # large base case, par. base

    sParBase = ph.Series()
    sParBase.xvalues = charsBenchmarks['threads']
    sParBase.yvalues = [np.mean(x) / 1e6 for x in dh.select(benchmarks, {'alg': 'd', 'parallel-base' : True, 'basecase': np.max(charsBenchmarks['basecase'])}).sort(columns='threads')['times']]
    sParBase.label = "Par. Base %i" % np.max(charsBenchmarks['basecase'])
    pRuntime.addSeries(sParBase)

    ssParBase = ph.Series(sParBase)
    ssParBase.yvalues = sParBase.yvalues[0] / sParBase.yvalues
    pSpeedupRel.addSeries(ssParBase)

    sASParBase = ph.Series(sParBase)
    sASParBase.yvalues = sCGAL.yvalue / sParBase.yvalues
    pSpeedupAbs.addSeries(sASParBase)

    pRuntime.plot("comparison_time_over_threads.png")
    pSpeedupRel.plot("comparison_rel_speedup_over_threads.png")
    pSpeedupAbs.plot("comparison_abs_speedup_over_threads.png")

def plotImprovement():
    ########################################################################################################################

    lastRunBenchmarks = dh.getLastRun(database, 'benchmarks')

    if lastRunBenchmarks < 2:
        # no improvement to plot
        return

    maxPointsNew = dh.selectMax(database, 'benchmarks', 'nP', {'run-number' : lastRunBenchmarks})
    maxPointsOld = dh.selectMax(database, 'benchmarks', 'nP', {'run-number' : lastRunBenchmarks-1})
    maxPoints = min(maxPointsOld, maxPointsNew)

    new = dh.load(database, 'benchmarks', {'run-number' : lastRunBenchmarks, 'nP' : maxPoints})
    old = dh.load(database, 'benchmarks', {'run-number' : lastRunBenchmarks-1,'nP' : maxPoints})

    charsBenchmarks = dh.getCharacteristics(new)

    pRuntime = ph.Plot()
    pRuntime.title = "Runtime over Threads"
    pRuntime.xlabel = "threads"
    pRuntime.ylabel = "time [ms]"
    pRuntime.desc = r"$10^%i$ points" % math.log10(charsBenchmarks['nP'][0])

    pImprovement = ph.Plot()
    pImprovement.title = "Improvement over Threads"
    pImprovement.xlabel = "threads"
    pImprovement.ylabel = "ratio"
    pImprovement.desc = r"$10^%i$ points" % math.log10(charsBenchmarks['nP'][0])

    #################################################
    # old small base case, seq. base

    sOldSeqBase = ph.Series()
    sOldSeqBase.xvalues = charsBenchmarks['threads']
    sOldSeqBase.yvalues = [np.mean(x) / 1e6 for x in dh.select(old, {'alg': 'd', 'parallel-base' : False, 'basecase' : np.min(charsBenchmarks['basecase'])}).sort(columns='threads')['times']]
    sOldSeqBase.label = "Old - Seq. Base %i" % np.min(charsBenchmarks['basecase'])
    pRuntime.addSeries(sOldSeqBase)


    #################################################
    # old large base case, par. base

    sOldParBase = ph.Series()
    sOldParBase.xvalues = charsBenchmarks['threads']
    sOldParBase.yvalues = [np.mean(x) / 1e6 for x in dh.select(old, {'alg': 'd', 'parallel-base' : True, 'basecase': np.max(charsBenchmarks['basecase'])}).sort(columns='threads')['times']]
    sOldParBase.label = "Old - Par. Base %i" % np.max(charsBenchmarks['basecase'])
    pRuntime.addSeries(sOldParBase)

    #################################################
    # new small base case, seq. base

    sNewSeqBase = ph.Series()
    sNewSeqBase.xvalues = charsBenchmarks['threads']
    sNewSeqBase.yvalues = [np.mean(x) / 1e6 for x in dh.select(new, {'alg': 'd', 'parallel-base' : False, 'basecase' : np.min(charsBenchmarks['basecase'])}).sort(columns='threads')['times']]
    sNewSeqBase.label = "New - Seq. Base %i" % np.min(charsBenchmarks['basecase'])
    pRuntime.addSeries(sNewSeqBase)


    #################################################
    # new large base case, par. base

    sNewParBase = ph.Series()
    sNewParBase.xvalues = charsBenchmarks['threads']
    sNewParBase.yvalues = [np.mean(x) / 1e6 for x in dh.select(new, {'alg': 'd', 'parallel-base' : True, 'basecase': np.max(charsBenchmarks['basecase'])}).sort(columns='threads')['times']]
    sNewParBase.label = "New - Par. Base %i" % np.max(charsBenchmarks['basecase'])
    pRuntime.addSeries(sNewParBase)

    sImpSeqBase = ph.Series(sNewSeqBase)
    sImpSeqBase.yvalues = [old / new for old, new in zip(sOldSeqBase.yvalues, sImpSeqBase.yvalues)]
    pImprovement.addSeries(sImpSeqBase)

    sImpParBase = ph.Series(sNewParBase)
    sImpParBase.yvalues = [old / new for old, new in zip(sOldParBase.yvalues, sImpParBase.yvalues)]
    pImprovement.addSeries(sImpParBase)

    pRuntime.plot("improvement_time_over_threads.png", {'ncol' : 2})
    pImprovement.plot("improvement_ratio_over_threads.png")

def plotProfiling():

    lastRun = dh.getLastRun(database, 'profiling')
    prof = dh.load(database, 'profiling', {'run-number' : lastRun}).sort(columns='nP')

    chars = dh.getCharacteristics(prof)
    counters = dh.getCounters(prof)

    stopCrit = dh.getStoppingCriterion(prof)

    #plot each counter over N
    for counter in counters:
        label = counter.replace('counter_', '')
        print("Plotting metric %s over n" % label)

        plt = ph.Plot()
        plt.title = r"%s over $n$" % label
        plt.xlabel = r"$n$"
        plt.ylabel = "%s" % label

        for sc in chars[stopCrit.field]:
            data = dh.select(prof, {stopCrit.field : sc})

            series = ph.Series()
            series.xvalues = data['nP']
            series.yvalues = data[counter]
            series.label = "%s: %i" % (stopCrit.label, sc)

            plt.addSeries(series)

        plt.plot("profiling/01_%s_over_n.png" % label)

    #plot stackplots for each stop crit over n
    # accumulate totals
    totalOps = {}

    for sc in chars[stopCrit.field]:
        print("Plotting stackplot for %s %i" % (stopCrit.label, sc))

        pCt = ph.StackPlot()
        pCt.title = r"Ops over $n$"
        pCt.xlabel = r"$n$"
        pCt.ylabel = "a.u."
        pCt.desc = "%s %i" % (stopCrit.label, sc)

        data = dh.select(prof, {stopCrit.field : sc})

        for counter in counters:
            label = counter.replace('counter_', '')

            series = ph.Series()
            series.xvalues = data['nP']
            series.yvalues = data[counter]
            series.label = label

            if not sc in totalOps:
                totalOps[sc] = copy.copy(series.yvalues)
            else:
                totalOps[sc] += series.yvalues

            pCt.addSeries(series)

        pCt.plot("profiling/02_ops_over_n_%s%i.png" % (stopCrit.shortLabel, sc),
                 figureArgs={'figsize': (11,6)},
                 axesArgs={'position' : (0.075, 0.1, 0.6, 0.8)},
                 legendArgs={'ncol' : 1, 'loc' : 5,
                             'fontsize' : 'small'}, figureLegend=True)

        pSh = ph.StackPlot()
        pSh.title = r"Ops over $n$"
        pSh.xlabel = r"$n$"
        pSh.ylabel = "%"
        pCt.desc = "%s %i" % (stopCrit.label, sc)

        for counter in counters:
            label = counter.replace('counter_', '')

            series = ph.Series()
            series.xvalues = data['nP']
            series.yvalues = data[counter] / totalOps[sc]
            series.label = label

            pSh.addSeries(series)

        pSh.plot("profiling/03_ops_perc_over_n_%s%i.png" % (stopCrit.shortLabel, sc),
                 figureArgs={'figsize': (11,6)},
                 axesArgs={'position' : (0.075, 0.1, 0.6, 0.8)},
                 legendArgs={'ncol' : 1, 'loc' : 5,
                             'fontsize' : 'small'}, figureLegend=True)

    # plot total operations
    plt = ph.Plot()
    plt.title = r"Total Operations over $n$"
    plt.xlabel = r"$n$"
    plt.ylabel = "ops"

    for sc in chars[stopCrit.field]:
        data = dh.select(prof, {stopCrit.field : sc})

        series = ph.Series()
        series.xvalues = data['nP']
        series.yvalues = totalOps[sc]
        series.label = "%s: %i" % (stopCrit.label, sc)

        plt.addSeries(series)

    plt.plot("profiling/01_totalOps_over_n.png")

    # plot box plots of basecase size over base case
    for sc in chars[stopCrit.field]:
        print("Plotting boxplot for %s %i" % (stopCrit.label, sc))

        pBS = ph.BoxPlot()
        pBS.title = r"Base Size over $n$"
        pBS.xlabel = r"$n$"
        pBS.ylabel = "$bc$"
        pBS.desc = "%s %i" % (stopCrit.label, sc)

        data = dh.select(prof, {stopCrit.field : sc})
        lN = np.sort(np.unique(data['nP']))

        series = ph.BoxPlotSeries()
        series.color = 'blue'
        series.marker = 'D'
        series.label = '%s %i' % (stopCrit.label, sc)

        series.xvalues = lN
        series.yvalues = []

        for n in lN:
            series.yvalues.append(list(dh.select(data, {'nP' : n})['basecaseSize']))

        pBS.addSeries(series)

        labelStep = 1
        if len(lN) > 10:
            labelStep = 4
        if len(lN) > 20:
            labelStep = 9

        pBS.plot("profiling/04_baseSize_over_n_%s%i" % (stopCrit.shortLabel, sc),
                 plotArgs={'labelStep' : labelStep})

if len(sys.argv) != 2:
    print("Specify database")
    sys.exit(4)

database = sys.argv[1]

# print("Plotting pure CGAL")
# plotCGAL()
#
# print("Plot Base Case")
# plotBaseCase()
#
# print("Plotting Comparision")
# plotComparison()
#
# print("Plotting Improvement")
# plotImprovement()

print("Plotting Profiling")
plotProfiling()