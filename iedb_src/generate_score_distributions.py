#! /usr/bin/python

import cPickle

from predict_binding import *
from generate_mhc_list import get_dic_model_list

'''For each method, generate prediction score distributions for 1million randomly selected peptides.
   These score distributions will be used by consensus predictor when combining scores.'''

def get_peptide_list(peptide_length=9):
    #fname = '/home/life/resource_biology/uniprot/peptide_list.txt'
    fname = '/home/life/resource_biology/uniprot/peptide_list_15.txt'
    f=open(fname,'r')
    lines=f.readlines()
    f.close()
    peptide_list = [line[0:peptide_length] for line in lines]
    return peptide_list

def get_dic_peptide_list():
    dic={}
    dic[8]  = get_peptide_list(peptide_length=8)
    dic[9]  = get_peptide_list(peptide_length=9)
    dic[10] = get_peptide_list(peptide_length=10)
    dic[11] = get_peptide_list(peptide_length=11)
    dic[12] = get_peptide_list(peptide_length=12)
    return dic

def get_method_list():
    method_list = ['ann','smm', 'smmpmbec', 'arb', 'comblib_sidney2008']
    #method_list = ['ann','smm','comblib_sidney2008']
    #method_list = ['smm','comblib_sidney2008']
    return method_list

def generate_score_distribution():
    dic_score_distributions = {}
    dic_peptide_list = get_dic_peptide_list()
    method_list = get_method_list()

    for method in method_list:
        version = '20090901'
        ms = MethodSet()
        method_index = ms.get_method_index(method)
        species_mhc_length_list = get_mhc_list(method_index)

        for (species, mhc, length) in species_mhc_length_list:
            print method, (species, mhc, length)
            peptide_list = dic_peptide_list[length]
            ###peptide_list = dic_peptide_list[length][0:1000]  # This is to speed up the test.
            scores = debug_predict_peptide_list(peptide_list, mhc, length, method)
            #print 'debug ', scores[0:5]
            scores = list(scores)
            scores.sort()
            scores_len = len(scores)
            scores_thousand = scores[0:scores_len:1000]  # scores_thousand should contain one thousand elements.
            ###scores_thousand = scores
            key = (method, mhc, length)
            #print 'scores_thousand ', scores_thousand[0:5]
            dic_score_distributions[key] = scores_thousand
            print key, len(scores_thousand), scores_thousand[0:5]
    return dic_score_distributions

def generate_score_distribution_training():
    '''Makes predictions for the training data. Stores returned scores to get (mean,std).'''
    sys.path.append('/home/life/workspace/consensus/src')
    from util_benchmark import BindingDataGeneral

    fname_bd = '/home/life/workspace/sharedlibs/retrieve_binding_data/binding_data/data_20090901/binding_data_mhc_i_final_20090813.txt'
    #fname_bd = ''
    dic_column = {'mhc':1, 'peptide_length':2, 'cv':3, 'sequence':4, 'inequality':5, 'meas':6}
    bd = BindingDataGeneral()
    bd.set_dic_column(dic_column)
    bd.read_data(fname_bd)

    dic_score_distributions = {}
    method_list = get_method_list()

    for method in method_list:
        version = '20090901'
        ms = MethodSet()
        method_index = ms.get_method_index(method)
        species_mhc_length_list = get_mhc_list(method_index)

        for (species, mhc, length) in species_mhc_length_list:
            key = bd.get_key(mhc,str(length))
            data_all        = bd.get_data_subset(key)
            data_binders    = bd.filter_binders(data_all)
            data_nonbinders = bd.filter_nonbinders(data_all)

            peptide_list            = [r[3].strip() for r in data_all]
            peptide_list_binders    = [r[3].strip() for r in data_binders]
            peptide_list_nonbinders = [r[3].strip() for r in data_nonbinders]

            print 'debug peptide_list',  peptide_list[0]

            print method, (species, mhc, length)
            #peptide_list = dic_peptide_list[length]
            ###peptide_list = dic_peptide_list[length][0:1000]  # This is to speed up the test.
            scores            = debug_predict_peptide_list(peptide_list, mhc, length, method)
            scores_binders    = debug_predict_peptide_list(peptide_list_binders, mhc, length, method)
            scores_nonbinders = debug_predict_peptide_list(peptide_list_nonbinders, mhc, length, method)


            #key = (method, mhc, length)
            method_key = (method, key)
            dic_row={}
            dic_row['all']        = list(scores)
            dic_row['binders']    = list(scores_binders)
            dic_row['nonbinders'] = list(scores_nonbinders)

            dic_score_distributions[method_key] = dic_row
            print method_key, len(scores), len(scores_binders), len(scores_nonbinders)
    return dic_score_distributions



def make_predictions_exclude():
    '''Makes predictions for the training data. Stores returned scores to get (mean,std).'''
    sys.path.append('/home/life/workspace/consensus/src')
    from util_benchmark import BindingDataGeneral

    #fname_bd = '/home/life/workspace/sharedlibs/retrieve_binding_data/binding_data/data_20090901/binding_data_mhc_i_final_20090813.txt'
    fname_bd = '/home/life/workspace/consensus_b/data/bd_exclude_20101007.txt'
    dic_column = {'mhc':0, 'peptide_length':1, 'cv':2, 'sequence':3, 'inequality':4, 'meas':5}
    bd = BindingDataGeneral()
    bd.set_dic_column(dic_column)
    bd.read_data(fname_bd)

    dic_score_distributions = {}
    method_list = get_method_list()

    for method in method_list:
        version = '20090901'
        ms = MethodSet()
        method_index = ms.get_method_index(method)
        species_mhc_length_list = get_mhc_list(method_index)

        for (species, mhc, length) in species_mhc_length_list:
            key = bd.get_key(mhc,str(length))
            if bd.has_key(key)==True:
                data_all        = bd.get_data_subset(key)

                peptide_list            = [r[3].strip() for r in data_all]


                print 'debug peptide_list',  peptide_list[0]

                print method, (species, mhc, length)
                #peptide_list = dic_peptide_list[length]
                ###peptide_list = dic_peptide_list[length][0:1000]  # This is to speed up the test.
                scores            = debug_predict_peptide_list(peptide_list, mhc, length, method)


                #key = (method, mhc, length)
                method_key = (method, key)
                dic_row={}
                dic_row['all']        = list(scores)


                dic_score_distributions[method_key] = dic_row
                print method_key, len(scores)
    return dic_score_distributions

def write_score_distribution(dic_score_distributions,fname):
    #dic_score_distributions = generate_score_distribution()
    #f=open('dic_distribution_scores.cPickle','w')
    f=open(fname,'w')
    cPickle.dump(dic_score_distributions, f)
    f.close()

def generate_model_list_consensus():
    dic_model_list = get_dic_model_list()
    model_list = []
    method_list = ['ann','smm','comblib_sidney2008']
    [model_list.extend(dic_model_list[method]) for method in method_list]
    model_list = list(set(model_list))
    model_list.sort()
    for modelname in model_list:
        available_method_list = []
        for method in method_list:
            if (modelname in dic_model_list[method]):
                available_method_list.append(method)
        print modelname, '\t', available_method_list



if __name__ == '__main__':
#    dic_score_distributions = generate_model_list_consensus()
#    write_score_distribution(dic_score_distributions, 'dic_distribution_scores.cpickle')


    dic_score_distributions = generate_score_distribution_training()
    write_score_distribution(dic_score_distributions, 'dic_distribution_scores_training.cpickle')

    dic_score_distributions = make_predictions_exclude()
    write_score_distribution(dic_score_distributions, 'dic_distribution_scores_exclude.cpickle')