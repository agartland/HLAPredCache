# Overview

HLAPredCache is a Python dictionary for storing HLA binding predictions for fast, programatic access. The dictionary has been augmented with methods that are specific to this use.

# Installation with IEDB backend

Download the set of [MHC_I binding tools](http://tools.immuneepitope.org/mhci/download/) from iedb.org.

Unzip and configure the tools.
```
tar -zxvf IEDB_MHC_I-2.13.tar.gz
cd mhc_i
./configure
```

TODO: How does this configuration set the correct paths for the cloned code?

Clone this repository (which includes updated versions of the IEDB Python source code.
```
git clone https://github.com/agartland/HLAPredCache
```

Test the installation.
```
python HLAPredCache/test_all.py
```
# Usage
```
import HLAPredCache

gag = 'MGARASVLSGGELDRWEKIRLRPGGKKKYKLKHIVWASRE'
mers = HLAPredCache.helpers.getMers(gag, nmers = [9])
ba = HLAPredCache.hlaPredCache()
ba.addPredictions(method = 'netmhcpan', alleles = ['A*0201','A*2402'], peptides = mers)

print ba[('A*0201','MGARASVLS')] 
```