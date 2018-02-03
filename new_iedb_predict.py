try:
    import skbio
except ImportError:
    pass
import numpy as np
import pandas as pd
import time
import requests
import io

__all__ = ['generateMersFromNT',
           'checkHLAs',
           'iedbPepPredict']
           
def generateMersFromNT(seqList, L=[8, 9, 10, 11]):
    """Take a list of nucleotide sequences,
    break into all peptides with lengths in L,
    expanding degenerate bases after peptides are generated, as needed.

    Parameters
    ----------
    seqList : list
        List of skbio.Sequence sequence objects
    L : list
        List of peptide lengths

    Returns
    -------
    mers : list
        Unique peptides in seqList"""
    L = np.array(L)
    Lnt = L*3
    mers = set()
    for i, seq in enumerate(seqList):
        print('Working on seq {} of {}'.format(i+1, len(seqList)))
        dna = skbio.sequence.DNA(seq).degap()
        for l in Lnt:
            for starti in range(0, len(dna)-l+1-1):
                dnaMer = dna[starti:starti + l]
                if dnaMer.has_degenerates():
                    for ex_mer in dnaMer.expand_degenerates():
                        mers.add(str(ex_mer.translate()))
                else:
                    mers.add(str(dnaMer.translate()))
    """Filter out peptides wit stop codon"""
    mers = [p for p in mers if not '*' in p]
    return mers


def checkHLAs(uHLAs, method='netmhcpan', lengths=[8, 9, 10, 11], verbose=False):
    testPeps = ['PAPGQEPRD', 'VHAGPAAPGQ', 'KYSRKKLTERS', 'HDKAHSMG']
    goodHLAs = []
    badHLAs = []

    for h in uHLAs:
        isBad = False
        for L in lengths:
            peptides = [p for p in testPeps if len(p) == L]
            try:
                resDf = iedbPepPredict([h], peptides, method=method)
            except:
                badHLAs.append((h, L))
                isBad = True
                continue
            if resDf.shape[0] != len(peptides):
                badHLAs.append((h, L))
                isBad = True
            elif resDf['ic50'].isnull().any():
                badHLAs.append((h, L))
                isBad = True
        if not isBad:
            goodHLAs.append(h)
            if verbose:
                print('{} good'.format(h))
        else:
            if verbose:
                print('{} bad'.format(h))
    return goodHLAs, badHLAs

def iedbPepPredict(hlas, peptides, method='netmhcpan', timeit=False):
    iedbURL = 'http://tools-cluster-interface.iedb.org/tools_api/mhci/'
    lengths = [len(p) for p in peptides]
    assert len(set(lengths)) == 1

    L = len(peptides[0])

    data = {'method':method,
            'email_address':'agartland@gmail.com',
            'sequence_text':' '.join(peptides),
            'sequence_format':'auto',
            'allele': ','.join(hlas),
            'length': ','.join([str(L)] * len(hlas))}
    startT = time.time()
    result = requests.post(iedbURL, data=data)
    endT = time.time()

    if result.content[:19].decode('utf-8') == 'Invalid allele name':
        print(result.content.decode('utf-8'))
        raise
    df = pd.read_csv(io.BytesIO(result.content), sep='\t')
    if timeit:
        print('Prediction of %d peptides took %1.1f minutes.' % (len(peptides), (endT - startT) / 60))
    
    if df.shape[0] != len(hlas) * len(peptides):
        print('Length of predictions ({}) is not as expected from number of HLAs ({}) and peptides ({})'.format(df.shape[0], len(hlas), len(peptides)))
    return df
