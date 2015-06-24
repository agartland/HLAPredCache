import unittest
import numpy as np

from cache import hlaPredCache
from predict import iedbPredict
from helpers import *

class TestHelpers(unittest.TestCase):
    def setUp(self):
        self.gag = 'MGARASVLSGGELDRWEKIRLRPGGKKKYKLKHIVWASRELERFAVNPGLLETSEGCRQILGQLQPSLQTGSEELRSLYNTVATLYCVHQRIEIKDTKEALDKIEEEQ'
        self.badpep = 'MGARXASZVLSGGE'
        self.gappedpep = 'MGAR-ASZVL-SGGE'
        self.h = 'A*0201'
    def tearDown(self):
        pass
    def test_getmers(self):
        mers = getMers(self.gag, nmers = [9])
        self.assertEqual(mers[0], 'MGARASVLS')

class TestCache(unittest.TestCase):
    def test_load(self):
        ba = hlaPredCache(baseFn = 'data/test', kmers = [9], warn = False, oldFile = False)
    def test_get(self):
        ba = hlaPredCache(baseFn = 'data/test', kmers = [9], warn = False, oldFile = False)
        self.assertEqual(ba[('A*2601','MGPGQVLFR')],10.3372161729)
        self.assertEqual(ba[('A*0203','ASRKLGDRG')],10.6185304083)
        self.assertEqual(ba[('A*0201','ASRKLGDRG')],10.7537776369)
        
        """ Actually it returns nan and does not raise an exception."""
        """
        with self.assertRaises(KeyError):
            test = ba[('A*2601','AGPGQVLFR')]
        """
        self.assertTrue(np.isnan(ba[('A*2601','AGPGQVLFR')]))
    def test_slice(self):
        ba = hlaPredCache(baseFn = 'data/test', kmers = [9], warn = False, oldFile = False)
        
        self.assertEqual(ba[('A*0203','ASRKLGDRG')],10.6185304083)
        self.assertEqual(ba[('A*0201','ASRKLGDRG')],10.7537776369)
        
        ba_slice = ba.slice(hlas = ['A*2601', 'A*3201', 'A*0201'], peptides = ['MGPGQVLFR', 'GSSSQVSRN','ASRKLGDRG'])
        self.assertEqual(ba_slice[('A*2601','MGPGQVLFR')],10.3372161729)
        self.assertEqual(ba_slice[('A*0201','ASRKLGDRG')],10.7537776369)
        #self.assertAlmostEqual(value, expected, places = 3)
        
        self.assertTrue(np.isnan(ba_slice[('A*0203','ASRKLGDRG')]))

class TestIEDBWrap(unittest.TestCase):
    def setUp(self):
        self.gag = 'MGARASVLSGGELDRWEKIRLRPGGKKKYKLKHIVWASRELERFAVNPGLLETSEGCRQILGQLQPSLQTGSEELRSLYNTVATLYCVHQRIEIKDTKEALDKIEEEQ'
        self.badpep = 'MGARXASZVLSGGE'
        self.gappedpep = 'MGAR-ASZVL-SGGE'
        self.h = 'A*0201'
        self.hlas = ['A*0201', 'A*0203', 'B*5701', 'A*2402']
    def test_dummy(self):
        mers = getMers(self.gag, nmers = [9])
        df = iedbPredict(method = 'RAND', hlas = self.hlas, peptides = mers[:10])
        self.assertEqual(df.shape[0], len(self.hlas) * 10)
        self.assertEqual(df['method'].iloc[0], 'RAND')
    def test_predict(self):
        mers = getMers(self.gag, nmers = [9])
        df = iedbPredict(method = 'netmhcpan', hlas = self.hlas, peptides = mers[:10])
        self.assertEqual(df.shape[0], len(self.hlas) * 10)
        self.assertEqual(df['method'].iloc[0], 'netmhcpan')
class TestIEDBSrc(unittest.TestCase):
    pass


if __name__ == '__main__':
    unittest.main()