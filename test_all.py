import unittest
import numpy as np

from cache import hlaPredCache, RandCache
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
        mers = getMers(self.gag, nmer = [9])
        self.assertEqual(mers[0], 'MGARASVLS')

class TestCache(unittest.TestCase):
    def setUp(self):
        self.gag = 'MGARASVLSGGELDRWEKIRLRPGGKKKYKLKHIVWASRELERFAVNPGLLETSEGCRQILGQLQPSLQTGSEELRSLYNTVATLYCVHQRIEIKDTKEALDKIEEEQ'
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
        self.assertFalse(ba.warn)

        self.assertEqual(ba[('A*0203','ASRKLGDRG')],10.6185304083)
        self.assertEqual(ba[('A*0201','ASRKLGDRG')],10.7537776369)
        
        ba_slice = ba.slice(hlas = ['A*2601', 'A*3201', 'A*0201'], peptides = ['MGPGQVLFR', 'GSSSQVSRN','ASRKLGDRG'])
        self.assertEqual(ba_slice[('A*2601','MGPGQVLFR')],10.3372161729)
        self.assertEqual(ba_slice[('A*0201','ASRKLGDRG')],10.7537776369)
        #self.assertAlmostEqual(value, expected, places = 3)

        self.assertFalse(ba_slice.warn)
        self.assertTrue(np.isnan(ba_slice[('A*0203','ASRKLGDRG')]))
    def test_add(self):
        ba = hlaPredCache()
        mers = getMers(self.gag, nmer = [9])
        nAdded = ba.addPredictions('netmhcpan', ['A*2601', 'A*3201', 'A*0201'], mers)
        """This peptide is a known A*02 binder"""
        self.assertTrue(ba[('A*0201','SLYNTVATL')] < np.exp(6))
        self.assertEqual(nAdded, 3 * len(mers))

class TestIEDBWrap(unittest.TestCase):
    def setUp(self):
        self.gag = 'MGARASVLSGGELDRWEKIRLRPGGKKKYKLKHIVWASRELERFAVNPGLLETSEGCRQILGQLQPSLQTGSEELRSLYNTVATLYCVHQRIEIKDTKEALDKIEEEQ'
        self.badpep = 'MGARXASZVLSGGE'
        self.gappedpep = 'MGAR-ASZVL-SGGE'
        self.h = 'A*0201'
        self.hlas = ['A*0201', 'A*0203', 'B*5701', 'A*2402']
        self.methods = ['ann', 'comblib_sidney2008', 'consensus', 'IEDB_recommended', 'netmhcpan', 'smm', 'smmpmbec', 'pickpocket', 'netmhccons']
    def test_dummy(self):
        mers = getMers(self.gag, nmer = [9])
        df = iedbPredict(method = 'RAND', hlas = self.hlas, peptides = mers[:10])
        self.assertEqual(df.shape[0], len(self.hlas) * 10)
        self.assertEqual(df['method'].iloc[0], 'RAND')
    def test_predict(self):
        mers = getMers(self.gag, nmer = [9])
        df = iedbPredict(method = 'netmhcpan', hlas = self.hlas, peptides = mers[:10])
        self.assertEqual(df.shape[0], len(self.hlas) * 10)
        self.assertEqual(df['method'].iloc[0], 'netmhcpan')
    def test_ann(self):
        self.method_test(method = 'ann')
    def test_comblib_sidney2008(self):
        self.skipTest('Method supports limited MHC alleles.')
        self.method_test(method = 'comblib_sidney2008')
    def test_netmhcpan(self):
        self.method_test(method = 'netmhcpan')
    def test_smm(self):
        self.method_test(method = 'smm')
    def test_smmpmbec(self):
        self.method_test(method = 'smmpmbec')
    def test_pickpocket(self):
        self.method_test(method = 'pickpocket')
    def test_netmhccons(self):
        self.method_test(method = 'netmhccons')
    def method_test(self,method):
        mers = getMers(self.gag, nmer = [9])
        cols = ['pred','method','peptide','hla']
        df = iedbPredict(method = method, hlas = self.hlas, peptides = mers[:10])
        self.assertEqual(df.shape[0], len(self.hlas) * 10)
        self.assertEqual(df['method'].iloc[0], method)
        self.assertTrue(np.all([c in df.columns for c in cols]))
class TestIEDBSrc(unittest.TestCase):
    pass

class TestRandomCache(unittest.TestCase):
    def test_init(self):
        ba = RandCache()
        self.assertEqual(str(ba),'RandCache (random)')
    def test_get(self):
        ba = RandCache()
        self.assertEqual(len(ba), 0)
        x = ba[('A*2601','MGPGQVLFR')]
        self.assertEqual(x, ba[('A*2601','MGPGQVLFR')])
        self.assertEqual(len(ba), 1)

if __name__ == '__main__':
    unittest.main()