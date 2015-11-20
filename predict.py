import numpy as np
import pandas as pd
import itertools
from multiprocessing import Pool
import parmap
import logging
import sys

import iedb_src.predict_binding as iedb_predict

__all__ = ['iedbPredict']

def convertHLAToIEDB(h):
    """Takes format A*1234 or A_1234 and returns A*12:34"""
    return 'HLA-' + h[:4].replace('_','*') + ':' + h[4:]
def convertHLABack(h):
    """Takes format HLA-A*12:34 and returns A_1234"""
    return h[4:].replace('*','_').replace(':','')

def iedbPredict(method, hlas, peptides, cpus=1, verbose=False):
    """Generate HLA:peptide binding affinity (log-IC50) predictions using
    the tools distributed by IEDB.

    Predictions are computed for all HLA:peptide combinations.

    Parameters
    ----------
    method : string
        Prediction method (e.g. netmhcpan, smm, ann)
        If RAND is specified then random predictions are returned.
    hlas : list
        List of HLA alleles in the format A_0201 or A*0201
    peptides : list of strings
        List of peptides, required to be all be of the same length.
    cpus : int
        Number of cores to use in parallelizing the predictions.

    Returns
    -------
    df : pd.DataFrame
        Columns: method, hla, peptide, core, pred"""

    if verbose:
        """Create console handler and set level to debug"""
        logging.basicConfig(level=logging.INFO, format='%(levelname)s:%(asctime)s:%(message)s')
        logging.info('HLA prediction initialized for %d HLA allele(s) using method %s on %d CPU(s)', len(hlas), method, cpus)


    cols = ['method','hla','peptide','core','pred']
    
    if method == 'RAND':
        results = dict(method=[], hla=[], peptide=[], core=[], pred=[])
        for h,pep in itertools.product(hlas,peptides):
            results['method'].append('RAND')
            results['hla'].append(h)
            results['peptide'].append(pep)
            results['core'].append(pep)
            results['pred'].append(np.random.rand())
        resDf = pd.DataFrame(results, columns=cols)
    else:
        if cpus > 1:
            result = parmap.map(_predictOneHLA, hlas, method, peptides, verbose, pool=Pool(processes=cpus))
        else:
            result = parmap.map(_predictOneHLA, hlas, method, peptides, verbose, parallel=False)

        """Remove None's"""
        resDf = pd.concat([r for r in result if not r is None], axis=0)

        """Take the log of the prediction if neccessary."""
        if resDf.pred.max() > 100:
            resDf['pred'] = np.log(resDf.pred)

        if verbose:
            logging.info('Completed %d predictions (expected %d)', resDf.shape[0], len(hlas) * len(peptides))
    return resDf


def _predictOneHLA(h, method, peptides, verbose):
    cols = ['method','hla','peptide','core','pred']
    try:
        resDf = iedb_predict.Prediction().predict(method, convertHLAToIEDB(h), peptides)
        resDf['hla'] = resDf.allele.map(convertHLABack)
        resDf['method'] = method
        resDf['core'] = resDf.peptide
        resDf['pred'] = resDf.ic50
        resDf = resDf[cols]
        if verbose:
            logging.info('Completed binding prediction for %s with %d peptides', h, len(peptides))
    except:
        exc_type, exc_value, exc_traceback = sys.exc_info()
        logging.warning('Prediction with allele %s generated exception %s: %s', h, exc_type, exc_value)
        resDf = None
    return resDf