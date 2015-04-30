__author__ = 'dfunke'

import sys
import pyejdb
import numpy as np
import pandas as pd

# util functions
def selectMax(database : str, collection : str, field : str, query : dict = None, addFields : list = []) -> pyejdb.bson.BSON_LazyDict:
    # Open database
    ejdb = pyejdb.EJDB(database, pyejdb.JBOREADER)

    hints = {}
    hints["$orderby"] = [(field, -1)]

    hints["$fields"] = {}
    hints["$fields"][field] = 1
    for f in addFields:
        hints["$fields"][f] = 1

    res = ejdb.findOne(collection, query, hints= hints)
    ejdb.close()

    if not res is None:
        print("max %s: %s" % (field, res.items()))
        if addFields:
            return res
        else:
            return res[field]
    else:
        return None

def getLastCommit(database : str, collection : str, query : dict = None) -> str:

    res = selectMax(database, collection, 'run-number', query, ['git-rev'])

    if not res is None:
        print("%s: %s" % (res['run-number'], res['git-rev']))
        return res['git-rev']
    else:
        return None

def getLastRun(database : str, collection : str, query : dict = None) -> str:

    res = selectMax(database, collection, 'run-number', query)

    if not res is None:
        print("Selecting run: %s" % res)
        return res
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
        return None

    return pd.DataFrame(lData).convert_objects(convert_dates=True, convert_numeric=True)


def select(dataset : pd.DataFrame, filters : dict) -> pd.DataFrame:
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

        # ignore counter fields
        if col.startswith('counter_'):
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