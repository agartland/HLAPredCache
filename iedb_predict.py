#!/usr/bin/env python
import argparse
import pandas as pd
from multiprocessing import Pool
from functools import partial
import sys
import logging

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

def predictHLA(h):
    args = parseArgs()

    with open(args.pep, 'r') as fh:
        peptides = [p.strip() for p in fh]

    try:
        outDf = iedbPredict(args.method, [h], peptides)
    except:
        exc_type, exc_value, exc_traceback = sys.exc_info()
        #lines = traceback.format_exception(exc_type, exc_value, exc_traceback)
        #print ''.join('!! ' + line for line in lines)  # Log it or whatever here
        logging.warning('Prediction with allele %s generated exception %s: %s', h, exc_type, exc_value)
        return None

    if args.verbose:
        logging.info('Completed binding prediction for %s with %d peptides', h, len(peptides))
    return outDf

if __name__ ==  '__main__':
    args = parseArgs()

    with open(args.hla, 'r') as fh:
        hlas = [h.strip() for h in fh]

    with open(args.pep, 'r') as fh:
        peptides = [p.strip() for p in fh]

    if args.verbose:
        """Create console handler and set level to debug"""
        logging.basicConfig(level = logging.INFO, format = '%(levelname)s:%(asctime)s:%(message)s')
        logging.info('HLA prediction initialized for %d HLA allele(s) using method %s on %d CPU(s)', len(hlas), args.method, args.cpus)

    if args.cpus > 1:
        pool = Pool(processes = args.cpus)
        result = pool.map(predictHLA, hlas)
    else:
        result = map(predictHLA, hlas)

    """Remove None's"""
    outDf = pd.concat([r for r in result if not r is None], axis = 0)
    
    outDf.to_csv(args.out)

    if args.verbose:
        logging.info('Completed %d predictions (expected %d) and wrote to file %s', outDf.shape[0], len(hlas) * len(peptides), args.out)
