#!/usr/bin/env python
"""
IEDB HLA:peptide binding prediction with multi-processing.

Example:
./iedb_predict.py --method netmhcpan --pep data/test.9.mers --hla data/test.my.hla --out data/test.9.out2 --verbose --cpus 4

TODO:
    (1) All predictions must fit in memory in a pd.DataFrame before being written to a file.
"""

import argparse
import pandas as pd
from multiprocessing import Pool
import parmap
import sys
import logging
import numpy as np
import mhctools

sys.path.append('/home/agartlan/gitrepo/')
from HLAPredCache.predict import iedbPredict

def parseArgs():
    parser = argparse.ArgumentParser(description='Predict HLA:peptide binding using multi-processing (one k)')
    parser.add_argument('--method', metavar='METHOD', type = str,
                       help='prediction method')
    parser.add_argument('--pep', metavar='PEPTIDE_FILE', type = str,
                       help='peptide filename')
    parser.add_argument('--hla', metavar='HLA_FILE', type = str,
                       help='HLA filename')
    parser.add_argument('--out', metavar='OUTPUT_FILE', type = str,
                       help='output filename')
    parser.add_argument('--cpus', metavar='CPUS', type = int, default = 1,
                       help='number of CPUs to utilize')
    parser.add_argument('--verbose', help='print status for each HLA allele', action = "store_true")
    args = parser.parse_args()
    return args

def predictHLA(h, method, peptides, verbose):
    try:
        outDf = iedbPredict(method, [h], peptides)
    except:
        exc_type, exc_value, exc_traceback = sys.exc_info()
        logging.warning('Prediction with allele %s generated exception %s: %s', h, exc_type, exc_value)
        return None

    if verbose:
        logging.info('Completed binding prediction for %s with %d peptides', h, len(peptides))
    return outDf

def predictHLA_mhctools(h, method, peptides, verbose):
    try:
        pred = mhctools.NetMHCpan4(alleles=[h], program_name='/fh/fast/gilbert_p/agartlan/gitrepo/iedb/mhc_i/method/netMHCpan-4.0/netMHCpan')
        outDf = pred.predict_peptides(peptides).to_dataframe()
    except:
        exc_type, exc_value, exc_traceback = sys.exc_info()
        logging.warning('Prediction with allele %s generated exception %s: %s', h, exc_type, exc_value)
        return None

    if verbose:
        logging.info('Completed binding prediction for %s with %d peptides', h, len(peptides))
    return outDf

def generatePredictions(method, hlas, peptides, cpus=1, verbose=False):
    """Does not work because peptides is also an iterator...."""
    if verbose:
        """Create console handler and set level to debug"""
        logging.basicConfig(level=logging.INFO, format='%(levelname)s:%(asctime)s:%(message)s')
        logging.info('HLA prediction initialized for %d HLA allele(s) using method %s on %d CPU(s)', len(hlas), method, cpus)

    if cpus > 1:
        result = parmap.map(predictHLA_mhctools, hlas, method, peptides, verbose, pool=Pool(processes=cpus))
    else:
        result = parmap.map(predictHLA_mhctools, hlas, method, peptides, verbose, parallel=False)

    """Remove None's"""
    outDf = pd.concat([r for r in result if not r is None], axis=0)

    """Take the log of the prediction if neccessary."""
    if outDf.affinity.max() > 100:
        outDf.loc[:, 'pred'] = np.log(outDf.affinity)

    if verbose:
        logging.info('Completed %d predictions (expected %d)', outDf.shape[0], len(hlas) * len(peptides))
    return outDf
   
if __name__ ==  '__main__':
    args = parseArgs()

    with open(args.hla, 'r') as fh:
        hlas = [h.strip() for h in fh]

    with open(args.pep, 'r') as fh:
        peptides = [p.strip() for p in fh]

    results = []
    lengths = sorted(list(set([len(p) for p in peptides])))
    for L in lengths:
        tmp = [p for p in peptides if len(p) == L]
        results.append(generatePredictions(args.method, hlas, tmp, cpus=args.cpus, verbose=args.verbose))
    predDf = pd.concat(results, axis=0)
    predDf.to_csv(args.out, index=False)

    if args.verbose:
        logging.info('Wrote %d predictions to file %s', predDf.shape[0], args.out)
