#! /usr/bin/python

import os

from setupinfo import *
from seqpredictor import *


def test_netmhcpan():
    method = 'netmhcpan'
    mhc    = 'HLA-A*02:01'
    length = 9
    input_protein_mhc = None
    sequence = 'MLVMAPRTVLLLLSAALALTETWAGSHSMRYFYTSVSRPGRGEPRFISVGYVDDTQFVRFDSDAASPREEPRAPWIEQEGPEYWDRNTQIYKAQAQTDRESLRNLRGYYNQSEAGSHTLQSMYGCDVGPDGRLLRGHDQYAYDGKDYIALNEDLRSWTAADTAAQITQRKWEAAREAEQRRAYLEGECVEWLRRYLENGKDKLERADPPKTHVTHHPISDHEATLRCWALGFYPAEITLTWQRDGEDQTQDTELVETRPAGDRTFQKWAAVVVPSGEEQRYTCHVQHEGLPKPLTLRWEPSSQSTVPIVGIVAGLAVLAVVVIGAVVAAVMCRRKSSGGKGGSYSQAACSDSAQGSDVSLTA'
    setup = SetupInfo()
    path_method = setup.path_method
    path_data = os.path.join(setup.path_data, method)
    predictor = NetMHCpanPredictor(path_method, path_data)
    predictor.initialize(mhc, length, input_protein_mhc=input_protein_mhc)
    scores = predictor.predict_sequence(sequence)
    print scores[0:5]
    assert False  # Check (1) there are multiple scores, (2) they are all numbers.

def test_netmhcpan_mhc():
    '''If userprovided mhc is true,'''
    method = 'netmhcpan'
    mhc = 'HLA-A*02:01'
    length = 9

    sequence = 'MLVMAPRTVLLLLSAALALTETWAGSHSMRYFYTSVSRPGRGEPRFISVGYVDDTQFVRFDSDAASPREEPRAPWIEQEGPEYWDRNTQIYKAQAQTDRESLRNLRGYYNQSEAGSHTLQSMYGCDVGPDGRLLRGHDQYAYDGKDYIALNEDLRSWTAADTAAQITQRKWEAAREAEQRRAYLEGECVEWLRRYLENGKDKLERADPPKTHVTHHPISDHEATLRCWALGFYPAEITLTWQRDGEDQTQDTELVETRPAGDRTFQKWAAVVVPSGEEQRYTCHVQHEGLPKPLTLRWEPSSQSTVPIVGIVAGLAVLAVVVIGAVVAAVMCRRKSSGGKGGSYSQAACSDSAQGSDVSLTA'
    input_protein_mhc = Proteins()
    input_protein_mhc.add_protein(sequence)

    setup = SetupInfo()
    path_method = setup.path_method
    path_data = os.path.join(setup.path_data, method)
    predictor = NetMHCpanPredictor(path_method, path_data)
    predictor.initialize(mhc, length, input_protein_mhc=input_protein_mhc)
    scores = predictor.predict_sequence(sequence)
    print scores[0:5]
    assert False


