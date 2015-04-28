#! /usr/bin/env python3

__author__ = 'dfunke'

import sys
import pyejdb
import numpy as np
import pandas as pd
import math

import plot_helpers as ph

# util functions
def _getLastEntry(database : str, collection : str, query : dict = None) -> pyejdb.bson.BSON_LazyDict:
    # Open database
    ejdb = pyejdb.EJDB(database, pyejdb.JBOREADER)

    res = ejdb.findOne(collection, query, hints= {
                   "$orderby" : [ ("run-number", -1) ],
                   "$fields" :  { "run-number" : 1, "git-rev" : 1 }})
    ejdb.close()

    return res

def getLastCommit(database : str, collection : str, query : dict = None) -> str:

    res = _getLastEntry(database, collection, query)

    if not res is None:
        print("%s: %s" % (res['run-number'], res['git-rev']))
        return res['git-rev']
    else:
        return None

def getLastRun(database : str, collection : str, query : dict = None) -> str:

    res = _getLastEntry(database, collection, query)

    if not res is None:
        print("Selecting run: %s" % res['run-number'])
        return res['run-number']
    else:
        return None

def load(database: str, collection: str, query : dict = None) -> pd.DataFrame:
    # Open database
    ejdb = pyejdb.EJDB(database, pyejdb.JBOREADER)

    lData = list()
    with ejdb.find(collection, query) as cur:

        for p in cur:
            lData.append(dict(p.items()))

    ejdb.close()

    if len(lData) == 0:
        print("%s-%s: no data selected by %s" % (database, collection, str(query)))
        sys.exit(4)

    return pd.DataFrame(lData).convert_objects(convert_dates=True, convert_numeric=True)


def select(dataset : pd.DataFrame, **filters) -> pd.DataFrame:
    mask = np.ones(len(dataset), dtype=np.bool)

    for field, value in filters.items():
        mask &= (dataset[field] == value)

    return dataset[mask]


def getCharacteristics(dataset : pd.DataFrame) -> dict:
    characteristics = dict()

    for col in dataset.columns:
        # skip over id and datetime column
        if col == '_id' or 'time' in col:
            continue

        # check whether column contains a plain characteristic
        if isinstance(dataset[col][0], list):
            continue

        nonNull = pd.notnull(dataset[col])

        if np.any(nonNull):
            characteristics[col] = np.sort(np.unique(dataset[col][nonNull]))

    # print(characteristics)
    return characteristics

def getSensitivityMap(dataset : pd.DataFrame, chars : dict, seriesChars : list) -> dict:
    map = dict()

    for series in seriesChars:
        map[series] = dict()
        for sValue in chars[series]:
            map[series][sValue] = dict()
            lData = select(dataset, **{series : sValue})

            for c in [x for x in chars if x != series]:
                map[series][sValue][c] = np.any(pd.notnull(lData[c]))

    return map

def plot(dataset : pd.DataFrame, chars : dict, xname : str, yname : str, seriesName : str, seriesValue : str, sensitiveChars : list, **filters):
    xvalues = np.sort(chars[xname])
    print(seriesName, "-", seriesValue, ": plot", yname, " over", xname, ":", xvalues ,"with filters", filters)

    query = filters.copy()
    query[seriesName] = seriesValue
    lData = select(dataset, **query).sort(xname)

    if len(lData) <= len(xvalues):
        # basecase - plot data

        if len(lData) == len(xvalues):
            yvalues = [np.mean(x) for x in lData[yname]]
        else:
            # we have less y values, create a mapping
            yvalues = list()
            for x in xvalues:
                if x in lData[xname]:
                    yvalues.append(np.mean(list(lData[lData[xname] == x][yname])))
                else:
                    yvalues.append(0)

        p = ph.Plot()
        p.title = "%s - %s: %s over %s" % (seriesName, seriesValue, yname, xname)

        p.xlabel = xname
        p.ylabel = yname
        p.desc = "Filters: %s" % str(filters)

        d = ph.Series()
        d.xvalues = xvalues
        d.yvalues = yvalues
        d.label = "%s: %s" % (seriesName, seriesValue)

        p.addSeries(d)

        # filename: seriesName-seriesValue_yname_xname_filters

        filterStr = "filters"
        for f in filters:
            filterStr += "_%s-%s" % (f, filters[f])

        p.plot("%s/%s-%s_%s_%s_%s.png" % (seriesName, seriesName, seriesValue, yname, xname, filterStr))

    else:
        # we need to reduce the dataset by further filters

        if not sensitiveChars:
            print("Dataset can't be reduced further - aborting")
            print(lData)
            sys.exit(4)

        # if len(sensitiveChars) == 1 and len(chars[sensitiveChars[0]]) <= 5 and len(lData) <= len(xvalues) * len(chars[sensitiveChars[0]]):
            # we have one sensitive characteristic left, with a sensible amount of catagories,
            # try to plot them as series


        cSensitiveChars = sensitiveChars.copy()
        filterName = cSensitiveChars.pop()

        for filterValue in chars[filterName]:
            nFilter = filters.copy()
            nFilter[filterName] = filterValue

            plot(dataset, chars, xname, yname, seriesName, seriesValue, cSensitiveChars, **nFilter)

def plotCGAL():
    lastRunCGAL = getLastRun(database, 'pureCGAL')

    pureCGAL = load(database, 'pureCGAL', {'run-number' : lastRunCGAL})
    charsCGAL = getCharacteristics(pureCGAL)

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
        sRuntime.yvalues = [np.mean(x) / 1e6 for x in select(pureCGAL, occupancy = occ).sort(columns='threads')['times']]

        sRuntime.label = "CGAL - Lock Granularity %s" % occ

        pRuntime.addSeries(sRuntime)

        sSpeedup = ph.Series(sRuntime)
        sSpeedup.yvalues = sSpeedup.yvalues[0] / sSpeedup.yvalues

        pSpeedup.addSeries(sSpeedup)

    pRuntime.plot("pureCGAL_time_over_threads.png")
    pSpeedup.plot("pureCGAL_speedup_over_threads.png")

def plotBaseCase():
    lastRunBenchmarks = getLastRun(database, 'benchmarks')

    benchmarks = load(database, 'benchmarks', {'run-number' : lastRunBenchmarks})
    charsBenchmarks = getCharacteristics(benchmarks)

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
        sRSeq.yvalues = [np.mean(x) / 1e6 for x in select(benchmarks, **{'basecase' : bs, 'parallel-base': False}).sort(columns='threads')['times']]

        sRSeq.label = "seq. Base %i" % bs

        pRuntime.addSeries(sRSeq)

        sSSeq = ph.Series(sRSeq)
        sSSeq.yvalues = sSSeq.yvalues[0] / sSSeq.yvalues

        pSpeedup.addSeries(sSSeq)
        
    for bs in charsBenchmarks['basecase']:
        # parallel base solver
        
        sRPar = ph.Series()
        sRPar.xvalues = charsBenchmarks['threads']
        sRPar.yvalues = [np.mean(x) / 1e6 for x in select(benchmarks, **{'basecase' : bs, 'parallel-base': True}).sort(columns='threads')['times']]

        sRPar.label = "par. Base %i" % bs

        pRuntime.addSeries(sRPar)

        sSPar = ph.Series(sRPar)
        sSPar.yvalues = sSPar.yvalues[0] / sSPar.yvalues

        pSpeedup.addSeries(sSPar)

    pRuntime.plot("basecase_time_over_threads.png", legend_cols=2)
    pSpeedup.plot("basecase_speedup_over_threads.png", legend_cols=2, legend_loc=2)

def plotComparison():
    ########################################################################################################################

    lastRunBenchmarks = getLastRun(database, 'benchmarks')

    benchmarks = load(database, 'benchmarks', {'run-number' : lastRunBenchmarks})
    charsBenchmarks = getCharacteristics(benchmarks)

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
    sCGAL.yvalue = np.mean(list(select(benchmarks, alg='c').sort(columns='threads')['times'])) / 1e6
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
    sCGALMT.yvalues = [np.mean(x) / 1e6 for x in select(benchmarks, alg='m').sort(columns='threads')['times']]
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
    sSeqBase.yvalues = [np.mean(x) / 1e6 for x in select(benchmarks, **{'alg': 'd', 'parallel-base' : False, 'basecase' : np.min(charsBenchmarks['basecase'])}).sort(columns='threads')['times']]
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
    sParBase.yvalues = [np.mean(x) / 1e6 for x in select(benchmarks, **{'alg': 'd', 'parallel-base' : True, 'basecase': np.max(charsBenchmarks['basecase'])}).sort(columns='threads')['times']]
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

    lastRunBenchmarks = int(getLastRun(database, 'benchmarks'))

    if lastRunBenchmarks < 2:
        # no improvement to plot
        return

    new = load(database, 'benchmarks', {'run-number' : str(lastRunBenchmarks)})
    old = load(database, 'benchmarks', {'run-number' : str(lastRunBenchmarks-1)})

    charsBenchmarks = getCharacteristics(new)

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
    sOldSeqBase.yvalues = [np.mean(x) / 1e6 for x in select(old, **{'alg': 'd', 'parallel-base' : False, 'basecase' : np.min(charsBenchmarks['basecase'])}).sort(columns='threads')['times']]
    sOldSeqBase.label = "Old - Seq. Base %i" % np.min(charsBenchmarks['basecase'])
    pRuntime.addSeries(sOldSeqBase)

    
    #################################################
    # old large base case, par. base

    sOldParBase = ph.Series()
    sOldParBase.xvalues = charsBenchmarks['threads']
    sOldParBase.yvalues = [np.mean(x) / 1e6 for x in select(old, **{'alg': 'd', 'parallel-base' : True, 'basecase': np.max(charsBenchmarks['basecase'])}).sort(columns='threads')['times']]
    sOldParBase.label = "Old - Par. Base %i" % np.max(charsBenchmarks['basecase'])
    pRuntime.addSeries(sOldParBase)
    
    #################################################
    # new small base case, seq. base

    sNewSeqBase = ph.Series()
    sNewSeqBase.xvalues = charsBenchmarks['threads']
    sNewSeqBase.yvalues = [np.mean(x) / 1e6 for x in select(new, **{'alg': 'd', 'parallel-base' : False, 'basecase' : np.min(charsBenchmarks['basecase'])}).sort(columns='threads')['times']]
    sNewSeqBase.label = "New - Seq. Base %i" % np.min(charsBenchmarks['basecase'])
    pRuntime.addSeries(sNewSeqBase)

    
    #################################################
    # new large base case, par. base

    sNewParBase = ph.Series()
    sNewParBase.xvalues = charsBenchmarks['threads']
    sNewParBase.yvalues = [np.mean(x) / 1e6 for x in select(new, **{'alg': 'd', 'parallel-base' : True, 'basecase': np.max(charsBenchmarks['basecase'])}).sort(columns='threads')['times']]
    sNewParBase.label = "New - Par. Base %i" % np.max(charsBenchmarks['basecase'])
    pRuntime.addSeries(sNewParBase)

    sImpSeqBase = ph.Series(sNewSeqBase)
    sImpSeqBase.yvalues = [old / new for old, new in zip(sOldSeqBase.yvalues, sImpSeqBase.yvalues)]
    pImprovement.addSeries(sImpSeqBase)
    
    sImpParBase = ph.Series(sNewParBase)
    sImpParBase.yvalues = [old / new for old, new in zip(sOldParBase.yvalues, sImpParBase.yvalues)]
    pImprovement.addSeries(sImpParBase)

    pRuntime.plot("improvement_time_over_threads.png", legend_cols=2)
    pImprovement.plot("improvement_ratio_over_threads.png")


if len(sys.argv) != 2:
    print("Specify database")
database = sys.argv[1]

print("Plotting pure CGAL")
plotCGAL()

print("Plot Base Case")
plotBaseCase()

print("Plotting Comparision")
plotComparison()

print("Plotting Improvement")
plotImprovement()