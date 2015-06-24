#! /usr/bin/python

import os

# For how many (MHC,length) do we have predictors?

def get_method_list():
    method_list = ['ann','arb','smm','smmpmbec','comblib_sidney2008']
    return method_list

def get_stats():
    mhc_length_list_list = []
    dir_prefix = '/home/life/workspace/sharedlibs/batch_mode_iedb_tools/mhc_i/data/MHCI_mhcibinding20090901/'
    method_list = get_method_list()
    for method in method_list:
        fname = os.path.join(dir_prefix, method, 'model_list.txt')
        f=open(fname,'r')
        lines=f.readlines()
        f.close()
        mhc_length_list = [line.split()[0].strip() for line in lines]
        mhc_length_list_list.extend(mhc_length_list)
    mhc_length_list_list_unique = list(set(mhc_length_list_list))
    mhc_list = ['-'.join(mhc_length.split('-')[0:-1]) for mhc_length in mhc_length_list_list_unique]
    mhc_list_unique = list(set(mhc_list)); mhc_list_unique.sort()
    print len(mhc_length_list_list_unique), len(mhc_length_list_list)
    print mhc_length_list_list_unique
    print 'len(mhc_list_unique) ', len(mhc_list_unique), len(mhc_list)
    print mhc_list_unique


if __name__ == '__main__':
    get_stats()