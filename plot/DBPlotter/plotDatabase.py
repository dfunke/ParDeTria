from dbm.dumb import _Database

__author__ = 'dfunke'

import sys
import pyejdb
import numpy as np
import pandas as pd
import matplotlib.pylab as plt

import plot_helpers

# util functions
def _getLastEntry(database : str, collection : str, query : dict = None) -> pyejdb.bson.BSON_LazyDict:
    # Open database
    ejdb = pyejdb.EJDB(database, pyejdb.JBOREADER)

    res = ejdb.findOne(collection, query, hints= {
                   "$orderby" : [ ("datetime", -1) ],
                   "$fields" :  { "datetime" : 1, "git-rev" : 1 }})
    ejdb.close()

    return res

def getLastCommit(database : str, collection : str, query : dict = None) -> str:

    res = _getLastEntry(database, collection, query)

    if not res is None:
        print("%s: %s" % (res['datetime'], res['git-rev']))
        return res['git-rev']
    else:
        return None

def getLastRun(database : str, collection : str, query : dict = None) -> str:

    res = _getLastEntry(database, collection, query)

    if not res is None:
        print("Selecting run from: %s" % res['datetime'])
        return res['datetime']
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
        if col == '_id' or col == 'datetime':
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

        fig, ax = plt.subplots()
        ax.plot(xvalues, yvalues, label="%s: %s" % (seriesName, seriesValue))
        ax.text(0,0, "Filters: %s" % str(filters), transform=ax.transAxes)

        ax.set_xlabel(xname)
        ax.set_ylabel(yname)
        ax.legend()

        # filename: seriesName-seriesValue_yname_xname_filters

        filterStr = "filters"
        for f in filters:
            filterStr += "_%s-%s" % (f, filters[f])

        plt.savefig("output/dbdump/%s-%s_%s_%s_%s.png" % (seriesName, seriesValue, yname, xname, filterStr))
        plt.close()

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



if len(sys.argv) != 2:
    print("Specify database")

database = sys.argv[1]
collection = 'benchmarks'
lastRun = getLastRun(database, collection)

df = load(database, collection, {'datetime' : lastRun})
chars = getCharacteristics(df)

seriesChars = ['alg', 'parallel-base']

# delete characteristics with only one entry
keys = list(chars.keys())
for chr in keys:
    if len(chars[chr]) == 1:
        del chars[chr]

sensitivity = getSensitivityMap(df, chars, seriesChars)
print("Sensitive parameters:", sensitivity)

# plot(df, chars, 'occupancy', 'times', 'alg', 'm', [], **{'nP': 10, 'threads': 4})

for series in seriesChars:
    print("Series:", series)

    for sValue in chars[series]:
        print("Value:", sValue)
        sensitiveChars = [c for c in chars if c not in seriesChars and sensitivity[series][sValue][c] == True]
        print("Sensitive Chars:", sensitiveChars)

        for xname in sensitiveChars:
            plot(df, chars, xname, 'times', series, sValue, [x for x in sensitiveChars if x != xname])