#! /usr/bin/python

from util import *

def test_predictor_selection():
    setupinfo = SetupInfo() # assumes version 20090901.
    #path_data = setupinfo.path_data
    #pred_method = 'ann'
    #species = 'human'
    #mhc = 'HLA A*0201'
    #length = '9'

    method_list = ['ann', 'smm', 'comblib_sidney2008']  # This is the predictors used to build 'consensus'
    #psb = PredictorSelectionB(path_data, pred_method, species, mhc, length)
    psb = PredictorSelectionB(setupinfo)
    model_list = psb.get_model_list_combined()
    print 'debug model_list ', len(model_list), model_list[0:5]
    content = []
    for mhc_length in model_list:
        print mhc_length, psb.is_model_available('ann', mhc_length)

#    for name in model_list:
#        method_list_available = psb.get_available_methods(name, method_list)
#        line = name+'\t'+str(method_list_available)
#        content.append(line)
#        #if (len(method_list_available) < 3):
#        #    print 'name ', name, method_list_available
#    content.sort()
#    for line in content:
#        print line
    assert False