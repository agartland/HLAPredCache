import unittest


class TestTools(unittest.TestCase):
    def setUp(self):
        self.gag = 'MGARASVLSGGELDRWEKIRLRPGGKKKYKLKHIVWASRELERFAVNPGLLETSEGCRQILGQLQPSLQTGSEELRSLYNTVATLYCVHQRIEIKDTKEALDKIEEEQ'
        self.badpep = 'MGARXASZVLSGGE'
        self.gappedpep = 'MGAR-ASZVL-SGGE'
        self.h = 'A*0201'
    def tearDown(self):
        pass
    def test_getmer(self):
        pass

class TestCache(unittest.TestCase)
    def test_load(self):
        ba = hlaPredCache(baseFn = 'data/test', kmers = [9], warn = True, oldFile = False)
    def test_get(self):
        ba = hlaPredCache(baseFn = 'data/test', kmers = [9], warn = True, oldFile = False)
        self.assertEqual(ba[('A*2601','MGPGQVLFR')],10.3372161729)
        self.assertEqual(ba[('A*0203','ASRKLGDRG')],10.6185304083)
        
        with self.assertRaises(KeyError):
            test = ba[('A*2601','AGPGQVLFR')]
    def test_slice(self):
        ba = hlaPredCache(baseFn = 'data/test', kmers = [9], warn = True, oldFile = False)
        
        self.assertEqual(ba[('A*0203','ASRKLGDRG')],10.6185304083)
        ba.slice(hlas = ['A*2601', 'A*3201', 'A*0201'], peptide = ['MGPGQVLFR', 'GSSSQVSRN','ASRKLGDRG'])
        self.assertEqual(ba[('A*2601','MGPGQVLFR')],10.3372161729)
        self.assertEqual(ba[('A*0201','ASRKLGDRG')],10.6185304083)
        #self.assertAlmostEqual(value, expected, places = 3)

        with self.assertRaises(KeyError):
            test = ba[('A*0203','ASRKLGDRG')]