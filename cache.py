import re
from functools import partial
import itertools
import numpy as np
import pandas as pd
from scipy import stats

from .helpers import *
from .predict import *

class hlaPredCache(dict):
    """Load 8,9,10,11-mer binding affinities into a big dictionary
    ba[(hla,peptide)]=9.1

    TODO:
     (1) Improve handling of class I and class II epitopes (and core info)
     (2) Integrate better with the prediction requester above"""
    def __init__(self, baseFn=None, kmers=[8, 9, 10, 11], warn=True, oldFile=False, useRand=False, newFile=False):
        dict.__init__(self)
        self.repAsteriskPattern = re.compile(r'\*')

        if oldFile:
            columnNames = ['method', 'hla', 'peptide', 'ic50']
            readCSVFunc = partial(pd.read_csv, names=columnNames, header=None)
        elif newFile:
            columnNames = ['method', 'hla', 'peptide', 'ic50']
            def _readfunc(fn):
                tmp = pd.read_csv(fn, usecols=['peptide', 'allele', 'prediction_method_name', 'pred'])
                tmp = tmp.rename({'prediction_method_name':'method', 'allele':'hla', 'pred':'ic50'}, axis=1)
                tmp.loc[:, 'hla'] = tmp.hla.str.slice(4, None)
                tmp.loc[:, 'hla'] = tmp.hla.map(partial(re.compile(r':').sub, ''))
                return tmp
            readCSVFunc = _readfunc
        else:
            columnNames = ['method', 'hla', 'peptide', 'core', 'ic50']
            readCSVFunc = partial(pd.read_csv, names=columnNames, skiprows=1)
        

        if not baseFn is None:
            self.name = baseFn.split('.')[0]
            self.predictionMethod = baseFn.split('.')[-1]
            fileList = ['%s.%d.out' % (baseFn, k) for k in kmers]
            self.fileList = fileList

            for fn in fileList:
                predDf = readCSVFunc(fn)
                if kmers == [15] and oldFile:
                    predDf['core'] = predDf.ic50.map(lambda s: s.split("'")[1])
                    predDf['ic50'] = predDf.ic50.map(lambda s: s[1:].split(",")[0])
                    predDf['ic50'] = predDf.ic50.map(float64)

                # predDf['hla'] = predDf.hla.map(partial(re.sub, self.repAsteriskPattern, '_'))
                predDf['hla'] = predDf.hla.map(partial(self.repAsteriskPattern.sub, '_'))
                self.update({(h, p):v for h, p, v in zip(predDf['hla'], predDf['peptide'], predDf['ic50'])})
        else:
            self.predictionMethod = ''
            self.name = ''

        self.warn = warn
        self.useRand = useRand

    def __str__(self):
        return '%s (%s)' % (self.name, self.predictionMethod)

    def __getitem__(self, key):
        """Calls the getItem() method with default options"""
        return self.getItem(key)
    def getItem(self, key, useRand = False):
        """Returns the requested prediction.
        First tries to return pred quickly, if not found then does minimal
        error checking for * in HLA or BADAA in mer.
        
        Warns for missing (hla,mer) keys before returning nan
            (can be suppressed with self.warn = False)

        Does not warn for invalid peptides, returns nan"""
        
        if self.useRand or useRand:
            hla, pep = key
            key = (hla, self.uPep[self.transMat[self.uPep.index(pep), self.uHLA.index(hla)]])
        try:
            """First, try to get the prediction the fastest way possible"""
            val = dict.__getitem__(self, key)
        except KeyError:
            hla, peptide = key
            hla = re.sub(self.repAsteriskPattern, '_', hla)

            if not isvalidmer(peptide):
                if self.warn:
                    pass
                    """Only want to warn about missing predictions, not invalid peptides"""
                    #print 'Invalid peptide %s, returning nan' % (peptide)
                val = np.nan
            else:
                try:
                    val = dict.__getitem__(self, (hla, peptide))
                except KeyError:
                    if self.warn:
                        print('HLA prediction not found (%s,%s), returning nan' % (hla, peptide))
                    val = np.nan
        return val
    def getRand(self, key):
        return self.getItem(key, useRand = True)
    def permutePeptides(self, seed = None):
        if not seed is None:
            np.random.seed(seed)
        uPep = list({k[1] for k in list(self.keys())})
        uHLA = list({k[0] for k in list(self.keys())})
        self.transMat = np.zeros((len(uPep), len(uHLA)), dtype = int)
        for i in np.arange(len(uHLA)):
            self.transMat[:, i] = np.random.permutation(np.arange(len(uPep), dtype = int))
        self.uPep = uPep
        self.uHLA = uHLA

    def __setitem__(self, key, val):
        dict.__setitem__(self, key, val)
    def addPredictions(self, method, hlas, peptides, cpus=1, verbose=False):
        """Run all neccessary predictions and add results to the cache without updating existing predictions
        Will batch peptides of each length.
        Will attempt to remove invalid peptides.
        Returns number of predictions added (counting only those that were new to the cache)"""
        if self.predictionMethod == '':
            self.predictionMethod = method
        if not method == self.predictionMethod:
            print('METHOD does not match existing method name for this cache')

        neededHLAs = set()
        neededPeptides = set()
        for h, pep in itertools.product(hlas, peptides):
            if (h, pep) not in self:
                neededHLAs.add(h)
                neededPeptides.add(pep)

        hlas = list(neededHLAs)
        """Remove bad peptides"""
        mers = [m for m in neededPeptides if isvalidmer(m)]

        kmers = {k: [m for m in mers if len(m)==k] for k in [8, 9, 10, 11, 12, 13, 14, 15]}

        for k in list(kmers.keys()):
            nAdded = 0
            if len(hlas) > 0 and len(kmers[k]) > 0:
                resDf = iedbPredict(method, hlas, kmers[k], cpus=cpus, verbose=verbose)
                self.update({(re.sub(self.repAsteriskPattern, '_', h), p):v for h, p, v in zip(resDf['hla'], resDf['peptide'], resDf['pred'])})
                nAdded += resDf.shape[0]
        return nAdded
    def addPredictionValues(self, hlas, peptides, values):
        """Add predictions as hla, peptide and values without running any predictor
        (basically just a dict update)"""
        self.update({(re.sub(self.repAsteriskPattern, '_', h), p):v for h, p, v in zip(hlas, peptides, values)})
    def dumpToFile(self, fn):
        with open(fn, 'w') as fh:
            for k, v in self.items():
                fh.write('%s,%s,%1.6f\n' % (k[0], k[1], v))
    def addFromFile(self, fn):
        """Add predictions from a file: hla,peptide,prediction
        Returns number of predictions added (all those in file w/o checking for duplicates)"""
        predDf = pd.read_csv(fn, names = ['hla', 'peptide', 'pred'], header = None)
        self.update({(h, p):v for h, p, v in zip(predDf['hla'], predDf['peptide'], predDf['pred'])})
        return predDf.shape[0]
    def slice(self, hlas, peptides):
        """Return a new hlaPredCache() with a subset of the predictions,
        identified by hlas and peptides"""
        out = hlaPredCache(warn = self.warn)
        out.update({(h, pep):self[(h, pep)] for h, pep in itertools.product(hlas, peptides)})
        return out

class RandCache(dict):
    """Starts as an empty hlaPredCache.
    As predictions are requested, random predictions are added to the cache.
    The effect is that you have random consistent/repeatable predictions.
    
    Good for testing code without actual predictions and for generating null
    distributions (i.e. "noisy predictor null")"""

    def __init__(self):
        dict.__init__(self)
        self.repAsteriskPattern = re.compile('\*')

        self.predictionMethod = 'random'
        self.name = 'RandCache'
        self.warn = False

    def __str__(self):
        return '%s (%s)' % (self.name, self.predictionMethod)

    def __getitem__(self, key):
        """Calls the getItem() method with default options"""
        return self.getItem(key)
    def getItem(self, key, useRand = False):
        """Returns the requested prediction."""
        try:
            """First, try to get the prediction"""
            val = dict.__getitem__(self, key)
        except KeyError:
            hla, peptide = key
            if not isvalidmer(peptide):
                val = np.nan
            else:
                hla = re.sub(self.repAsteriskPattern, '_', hla)
                try:
                    val = dict.__getitem__(self, (hla, peptide))
                except KeyError:
                    dict.__setitem__(self, (hla, peptide), self._generateNewPrediction())
                    val = dict.__getitem__(self, (hla, peptide))
        return val
    def _generateNewPrediction(self):
        return np.abs(11 - stats.expon.rvs(0, 1.5, size = 1))[0]
    def getRand(self, key):
        """Here to preserve the interface, but does nothing functionally different"""
        return self.getItem(key, useRand = False)
    def permutePeptides(self, seed = None):
        """Here to preserve the interface, but does nothing functionally different"""
    def __setitem__(self, key, val):
        dict.__setitem__(self, key, val)
    def addPredictions(self, method, hlas, peptides):
        """Here to preserve the interface, but does nothing."""
        pass
    def addPredictionValues(self, hlas, peptides, values):
        """Here to preserve the interface, but does nothing."""
        pass
    def dumpToFile(self, fn):
        with open(fn, 'w') as fh:
            for k, v in self.items():
                fh.write('%s,%s,%1.6f\n' % (k[0], k[1], v))
    def addFromFile(self, fn):
        """Here to preserve the interface, but does nothing."""
        pass
    def slice(self, hlas, peptides):
        """Return a new hlaPredCache() with a subset of the predictions,
        identified by hlas and peptides"""
        out = hlaPredCache(warn = self.warn)
        out.update({(h, pep):self[(h, pep)] for h, pep in itertools.product(hlas, peptides)})
        return out