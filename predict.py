import numpy as np
import pandas as pd
import itertools

import time
from  tempfile import NamedTemporaryFile, mktemp
import subprocess
import os

try:
    import mhc_bindings.mhc_bindings as mb
except ImportError:
    mb = None

__all__ = ['iedbPredict',
           'predictHLABinding']

def iedbPredict(method,hlas,peptides):
    results = dict(method = [], hla = [], peptide = [], core = [], pred = [])
    if method == 'RAND':
        for h,pep in itertools.product(hlas,peptides):
            results['method'].append('RAND')
            results['hla'].append(h)
            results['peptide'].append(pep)
            results['core'].append(pep)
            results['pred'].append(np.random.rand())
    else:
        pass

    resDf = pd.DataFrame(results, columns = ['method','hla','peptide','core','pred'])
    return resDf

def predictHLABinding(method,hlas,peptides,useTempFiles=False,verbose=False,dumbledore=False,query=True,force=False,no_cache=False,cpus=1):
    """Use the new predictors and return a resDf with
    columns: method, hla, peptide, core, missing"""

    """Ensure that hlas and peptides are unique"""
    hlas = list(set(hlas))
    peptides = list(set(peptides))

    if verbose:
        print 'Predicting %d mers with %d HLAs...' % (len(peptides),len(hlas)),
    """Run all predictions and add results to the cache"""
    if not (useTempFiles or dumbledore):
        results = mb.fetch_return(methods = [method],
                                mhc_strs = hlas,
                                peptide_strs = peptides,
                                peptide_length = len(peptides[0]),
                                force = force,
                                query = query,
                                verb = verbose,
                                no_cache = no_cache,
                                cpus = cpus,
                                db_passwd = 'R1ghtN0w',login='agartlan')
        """Add the missing field"""
        #results=[(m,h,pep,c,pred,False) for m,h,pep,c,pred in results]

        resDf = pd.DataFrame(results, columns = ['method','hla','peptide','core','pred'])
    else:
        """Write peptides to a tempfile"""
        with NamedTemporaryFile(delete=False, mode='w', prefix='agtmp', suffix='.mers') as pepFh:
            pepFn = pepFh.name
            for pep in peptides:
                pepFh.write('%s\n' % pep)

        """Write HLAs to a tempfile"""
        with NamedTemporaryFile(delete=False,mode='w',prefix='agtmp',suffix='.hla') as hlaFh:
            hlaFn=hlaFh.name
            for h in hlas:
                hlaFh.write('%s\n' % h)

        """Create an out tempfile"""
        """with NamedTemporaryFile(delete=False,mode='w',prefix='agtmp',suffix='.out') as outFh:
            outFn=outFh.name"""
        outFn=mktemp(suffix='.out', prefix='agtmp')


        """Get predictions via tempfiles using the old style dumbledore predictor"""
        if dumbledore:
            cpus = len(hlas) if len(hlas)<10 else 10
            columnNames=['method','hla','peptide','ic50']
            """/home/agartlan/epd72/bin/python ~/nethome/scripts/predictorDB/fetch_energies.py -m NetMHCpan -L 9 ~/nethome/data/epitope_mapping/shuff.9.mers ~/nethome/data/HLA/common.hla ~/nethome/data/epitope_mapping/shuff.netmhcpan.9.out"""
            cmd=['/home/agartlan/epd72/bin/python',
                 '/home/agartlan/nethome/scripts/predictorDB/fetch_energies.py',
                 '-m%s' % method,
                 '-L%d' % len(peptides[0]),
                 '-j%d' % cpus,
                 pepFn,
                 hlaFn,
                 outFn]

            """Blocking call to the predictor"""
            print ' '.join(cmd)
            mhc_bindingsProc = subprocess.call(cmd)

            if mhc_bindingsProc == 0:
                """Read in the out file"""
                resDf = pd.read_csv(outFn,names=columnNames)
                resDf['hla']=resDf.hla.map(partial(re.sub,'[*]','_'))
                resDf = resDf.rename_axis(mapper={'ic50':'pred'},axis=1)
            else:
                resDf = pd.DataFrame()
        
        else:
            """Get predictions via tempfiles using the new predictor cache written by Dave Swan (works on compusrv and rhinos)"""
            cpus = len(hlas) if len(hlas)<6 else 6
            """Requires that mhc_bindings.sh be on the PATH"""
            cmd=['mhc_bindings.sh',
                 '-b%s' % method,
                 '-a%s' % hlaFn,
                 '-s%s' % pepFn,
                 '-L%d' % len(peptides[0]),
                 '-j%d' % cpus,
                 '-o%s' % outFn]
            if verbose:
                cmd.append('--verbose')

            """Blocking call to the predictor"""
            mhc_bindingsProc=subprocess.call(cmd)

            if mhc_bindingsProc == 0:
                """Read in the out file"""
                resDf = pd.read_csv(outFn)
                resDf = resDf.rename_axis(axis=1, mapper={'Value':'pred'})
                resDf = resDf.rename_axis(axis=1, mapper=string.lower)
            else:
                resDf = pd.DataFrame()

        """Reomve all the temporary files"""
        os.unlink(pepFn)
        os.unlink(hlaFn)
        os.unlink(outFn)

    expectedN = len(hlas)*len(peptides)
    actN = resDf.shape[0]
    if actN == expectedN:
        if verbose:
            print 'SUCESS: added %d predictions' % actN
    else:
        print 'WARNING: Of %d expected, missing %d predictions (added %d)' % (expectedN,expectedN - actN,actN)
        expectedKeys = [(h,p) for h,p in itertools.product(hlas,peptides)]
        if expectedN<actN:
            print 'Extra peptides:'
            for i,row in resDf.iterrows():
                if not (row['hla'],row['peptide']) in expectedKeys:
                    print '\t%s : %s' % (row['hla'],row['peptide'])
        else:
            print 'Missing peptides:'
            for h,p in expectedKeys:
                if ((resDf.hla == h) & (resDf.peptide==p)).sum()<1:
                    print '\t%s : %s' % (h,p)
        
        print 'Dropping duplicates...'
        duplicatesInd = resDf.duplicated(cols=['hla','peptide'])
        print resDf.loc[duplicatesInd]
        duplicatesInd = resDf.duplicated(cols=['hla','peptide'],take_last=True)
        print resDf.loc[duplicatesInd]
        resDf = resDf.drop_duplicates(cols=['hla','peptide'])
        actN = resDf.shape[0]
        if actN == expectedN:
            print 'DUPLICATES: added %d predictions' % actN
        else:
            print 'NOT DUPLICATES: Of %d expected, missing %d predictions (added %d)' % (expectedN,expectedN - actN,actN)
    return resDf
