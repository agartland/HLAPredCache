#! /usr/bin/python

#import sys
import os

from util import get_standard_mhc_name

def get_method_list():
    method_list = ['ann','smm','arb', 'comblib_sidney2008']
    return method_list

def get_dic_model_list():
    dic = {}
    method_list = get_method_list()
    for method in method_list:
        fname = os.path.join('./data/MHCI_mhcibinding20090901/', method, 'model_list.txt')
        f=open(fname,'r')
        lines=f.readlines()
        model_list = [line.split()[0].strip() for line in lines]
        dic[method] = model_list
    return dic

def write_mhc_list():
    dic = get_dic_model_list()
    method_list = get_method_list()
    for method in method_list:
        print 'A list of available (MHC, PeptideLength) for ', method
        print 'MHC\tPeptideLength'
        model_list = dic[method]
        for model_name in model_list:
            (mhc, peptide_length) = get_standard_mhc_name(model_name)
            #temp = model_name.strip().split('-')
            #peptide_length = temp[-1]
            #mhc = '-'.join(temp[0:-1])
            print mhc,'\t', peptide_length
        print '\n'


if __name__ == '__main__':
    write_mhc_list()
