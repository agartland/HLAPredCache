#! /usr/bin/python
import sys
sys.path.append('./src')
from predict_binding import *

def test_ann():
    method_test = 'ann'
    (mhc_test, peptide_length_test, peptide_test, score_test) = simple_prediction(input_sequence = 'AAAVVVLLLGGG', mhc='HLA A*0201', length=9, method=method_test)
    print method_test, mhc_test, peptide_length_test, peptide_test, score_test
    assert False

def test_smm():
    method_test = 'smm'
    (mhc_test, peptide_length_test, peptide_test, score_test) = simple_prediction(input_sequence = 'AAAVVVLLLGGG', mhc='HLA A*0201', length=9, method=method_test)
    print method_test, mhc_test, peptide_length_test, peptide_test, score_test
    assert False

def test_arb():
    method_test = 'arb'
    (mhc_test, peptide_length_test, peptide_test, score_test) = simple_prediction(input_sequence = 'AAAVVVLLLGGG', mhc='HLA A*0201', length=9, method=method_test)
    print method_test, mhc_test, peptide_length_test, peptide_test, score_test
    assert False

def test_smmpmbec():
    method_test = 'smmpmbec'
    (mhc_test, peptide_length_test, peptide_test, score_test) = simple_prediction(input_sequence = 'AAAVVVLLLGGG', mhc='HLA A*0201', length=9, method=method_test)
    print method_test, mhc_test, peptide_length_test, peptide_test, score_test
    assert False

def test_comblibsidney2008():
    method_test = 'comblib_sidney2008'
    (mhc_test, peptide_length_test, peptide_test, score_test) = simple_prediction(input_sequence = 'AAAVVVLLLGGG', mhc='HLA A*0201', length=9, method=method_test)
    print method_test, mhc_test, peptide_length_test, peptide_test, score_test
    assert False

def test_netmhcpan():
    method_test = 'netmhcpan'
    (mhc_test, peptide_length_test, peptide_test, score_test) = simple_prediction(input_sequence = 'AAAVVVLLLGGG', mhc='HLA A*0201', length=9, method=method_test)
    print method_test, mhc_test, peptide_length_test, peptide_test, score_test
    assert False

def test_netmhcpan_usermhc():
    method_test = 'netmhcpan'
    usermhc_seq = 'MLVMAPRTVLLLLSAALALTETWAGSHSMRYFYTSVSRPGRGEPRFISVGYVDDTQFVRFDSDAASPREEPRAPWIEQEGPEYWDRNTQIYKAQAQTDRESLRNLRGYYNQSEAGSHTLQSMYGCDVGPDGRLLRGHDQYAYDGKDYIALNEDLRSWTAADTAAQITQRKWEAAREAEQRRAYLEGECVEWLRRYLENGKDKLERADPPKTHVTHHPISDHEATLRCWALGFYPAEITLTWQRDGEDQTQDTELVETRPAGDRTFQKWAAVVVPSGEEQRYTCHVQHEGLPKPLTLRWEPSSQSTVPIVGIVAGLAVLAVVVIGAVVAAVMCRRKSSGGKGGSYSQAACSDSAQGSDVSLTA'
    (mhc_test, peptide_length_test, peptide_test, score_test) = simple_prediction(input_sequence = 'AAAVVVLLLGGG', mhc='HLA A*0201', length=9, method=method_test, usermhc_seq = usermhc_seq)
    print method_test, mhc_test, peptide_length_test, peptide_test, score_test
    assert False

def test_consensus():
    method_test = 'consensus'
    (mhc_test, peptide_length_test, peptide_test, score_test) = simple_prediction(input_sequence = 'AAAVVVLLLGGG', mhc='HLA A*0201', length=9, method=method_test)
    print method_test, mhc_test, peptide_length_test, peptide_test, score_test
    assert False



