import numpy as np
import string
import re

__all__ = ['BADAA',
            'AALPHABET',
            'convertHLAAsterisk',
            'isvalidmer',
            'isvalidHLA',
            'rankEpitopes',
            'rankKmers',
            'rankMers',
            'getIC50',
            'getMers',
            'getMerInds',
            'grabKmer',
            'grabKmerInds',
            'findpeptide',
            'grabOverlappingKmer',
            'overlappingMers']

BADAA = '-*BX#Z'
AALPHABET = 'ACDEFGHIKLMNPQRSTVWY'


def convertHLAAsterisk(hlas):
    """Replace the * with _ in each HLA allele"""
    repAsteriskPattern = re.compile('\*')
    return [re.sub(repAsteriskPattern, '_', h) for h in hlas]

def isvalidmer(mer):
    if not mer is None:
        return not re.search('[%s]' % BADAA,mer)
    else:
        return False

def isvalidHLA(h):
    if h[0] in 'AB':
        return True
    else:
        return False

def rankEpitopes(ba, hlaList, peptide, nmer = [8,9,10,11], peptideLength = None):
    """Breaks peptide into kmers (all nmer lengths)
    and rank all (hla, kmer) pairs by predicted IC50 in hlaPredCache ba

    IDENTICAL to rankKmers but may have different performance?

    Can be used to find the most likely optimal epitope in a peptide sequence.

    Predictions that are not found in ba get a temporary prediction of 15 log-nM

    Parameters
    ----------
    ba : hlaPredCache
        dict-like container of all (hla, kmer) IC50 values
    hlaList : list
        HLA alleles to be used as keys in ba
    peptide : str
        AA sequence
    nmer : list
        Integers indicating optimal lengths to be tested as kmers.
    peptideLength : int or None
        If a number is specified then a number of '.' padded kmers are included
        so that there are always garaunteed to be a certain number of kmers and results

    Returns
    -------
    ranks : ndarray int
        Zero-based rankings of kmers based on predicted IC50 (lowest IC50, lowest rank)
    sorti : ndarray int
        Index that can be used to sort the returned arrays
    kmers : ndarray object
        Array of kmer strings in order by getMers() (can be sorted by rank with sorti)
    ic50 : ndarray float
        Predicted log-IC50 (log-nM) with the HLA allele with the lowest IC50
    hla : ndarray object
        Array of HLA alleles that were the best predicted binder to each kmer"""

    merList = getMers(peptide, nmer, peptideLength)
    kmers = np.empty((len(merList),len(hlaList)), dtype=object)
    ic50 = np.ones((len(merList),len(hlaList))) * 15
    hla = np.empty((len(merList),len(hlaList)), dtype=object)
    for i,m in enumerate(merList):
        for j,h in enumerate(hlaList):
            kmers[i,j] = m
            hla[i,j] = h
            tmp = ba[(h,m)]
            if not np.isnan(tmp):
                ic50[i,j] = tmp
    kmers = kmers.flatten()
    ic50 = ic50.flatten()
    hla = hla.flatten()
    sorti = ic50.argsort()
    ranks = np.empty(len(ic50), int)
    ranks[sorti] = np.arange(len(ic50))
    return (ranks,sorti,kmers,ic50,hla)

def rankKmers(ba, hlaList, peptide, nmer = [8,9,10,11], peptideLength = None):
    """Breaks peptide into kmers (all nmer lengths)
    and rank all (hla, kmer) pairs by predicted IC50 in hlaPredCache ba

    IDENTICAL to rankEpitopes but may have different performance?

    Can be used to find the most likely optimal epitope in a peptide sequence.

    Predictions that are not found in ba get a temporary prediction of 15 log-nM

    Parameters
    ----------
    ba : hlaPredCache
        dict-like container of all (hla, kmer) IC50 values
    hlaList : list
        HLA alleles to be used as keys in ba
    peptide : str
        AA sequence
    nmer : list
        Integers indicating optimal lengths to be tested as kmers.
    peptideLength : int or None
        If a number is specified then a number of '.' padded kmers are included
        so that there are always garaunteed to be a certain number of kmers and results

    Returns
    -------
    ranks : ndarray int
        Zero-based rankings of kmers based on predicted IC50 (lowest IC50, lowest rank)
    sorti : ndarray int
        Index that can be used to sort the returned arrays
    kmers : ndarray object
        Array of kmer strings in order by getMers() (can be sorted by rank with sorti)
    ic50 : ndarray float
        Predicted log-IC50 (log-nM) with the HLA allele with the lowest IC50
    hla : ndarray object
        Array of HLA alleles that were the best predicted binder to each kmer"""
    kmers = getMers(peptide, nmer, peptideLength)
    result = rankMers(ba, hlaList, kmers)
    return (result[0],result[1],kmers,result[2],result[3])
    
def rankMers(ba, hlaList, merList):
    """Ranks all (hla, mer) pairs by predicted IC50 found in hlaPredCache, ba

    Can be used to find the most likely optimal epitope from a list.

    Predictions that are not found in ba get a temporary prediction of 15 log-nM

    Parameters
    ----------
    ba : hlaPredCache
        dict-like container of all (hla, kmer) IC50 values
    hlaList : list
        HLA alleles to be used as keys in ba
    merList : list
        Peptide sequences to be tests with each HLA allele

    Returns
    -------
    ranks : ndarray int
        Zero-based rankings of kmers based on predicted IC50 (lowest IC50, lowest rank)
    sorti : ndarray int
        Index that can be used to sort the returned arrays
    kmers : ndarray object
        Array of kmer strings in order by getMers() (can be sorted by rank with sorti)
    ic50 : ndarray float
        Predicted log-IC50 (log-nM) with the HLA allele with the lowest IC50
    hla : ndarray object
        Array of HLA alleles that were the best predicted binder to each kmer"""

    ic50 = np.ones((len(merList)))*15
    hla = np.empty(len(merList), dtype=object)
    for i,m in enumerate(merList):
        if not '.' in m:
            ic50[i],hla[i],dumpmer = getIC50(ba, hlaList, m)
    sorti = ic50.argsort()
    ranks = np.empty(len(ic50), int)
    ranks[sorti] = np.arange(len(ic50))
    return (ranks,sorti,ic50,hla)

def getIC50(ba, hlaList, mer, nmer = [8,9,10,11]):
    """Return the IC50 from ba of the mer and its affinity with the most avid HLA in hlaList.
    Or if len(pep)>11, return that of the most avid kmer

    Parameters
    ----------
    ba : hlaPredCache
        dict-like container of all (hla, kmer) IC50 values
    hlaList : list
        HLA alleles to be used as keys in ba
    mer : string
        Peptide sequences to be tests with each HLA allele
    nmer : list
        Integers indicating optimal lengths to be tested as kmers."""

    if ba is None:
        raise NameError('Did not load IC50 values into ba!')
    #minimum IC50 over the HLAs
    if len(mer) <= 11:
        allPairs = [(ba[(h,mer)],h,mer) for h in hlaList]
    #minimum IC50 over all the mers and all the HLAs
    else:
        allPairs = [getIC50(ba,hlaList,m) for m in getMers(mer,nmer)]
    return min(allPairs, key = lambda x: x[0])


def getMers(seq, nmer = [8, 9 , 10, 11], seqLength = None):
    """Takes a AA sequence (string) and turns it into a list of 8, 9, 10, 11 mers
    
    The seq will be padded with one or more '.' if it is shorter than seqLength
    These indices will match the peptides created by getMers()

    Paramters
    ---------
    seq : str
        Peptide sequence.
    nmer : list
        List of k's for the creation of all kmers.
    seqLength : int
        Minimum length of seq ('.' used for padding before applying the process)
        Useful for garaunteeing that a certain number of kmers will be in the list.

    Returns
    -------
    mers : list
        All peptides of length nmer contained by seq"""
    if not seqLength is None:
        if len(seq) > seqLength:
            seq = seq[:seqLength]
        elif len(seq) < seqLength:
            seq = string.ljust(seq, seqLength, '.')

    mers = []
    for n in nmer:
        mers.extend([seq[i:i+n] for i in range(len(seq)-n+1)])
    return mers

def getMerInds(seq, nmer = [8, 9 , 10, 11], seqLength = None):
    """Takes a AA sequence (string) and turns it into a list of 8, 9, 10, 11 mers
    
    The seq will be padded with one or more '.' if it is shorter than seqLength
    These indices will match the peptides created by getMers()

    Paramters
    ---------
    seq : str
        Peptide sequence.
    nmer : list
        List of k's for the creation of all kmers.
    seqLength : int
        Minimum length of seq ('.' used for padding before applying the process)
        Useful for garaunteeing that a certain number of kmers will be in the list.

    Returns
    -------
    mers : list
        All peptides of length nmer contained by seq
    mers : list
        Seq indices for mers"""
    if not seqLength is None:
        if len(seq) > seqLength:
            seq = seq[:seqLength]
        elif len(seq) < seqLength:
            seq = string.ljust(seq, seqLength, '.')

    mers = []
    inds = []
    for n in nmer:
        mers.extend([seq[i:i+n] for i in range(len(seq)-n+1)])
        inds.extend([np.arange(n)+i for i in range(len(seq)-n+1)])
    return mers,inds

def itermer(seq, k=9, gapped=True, yield_inds=False):
    """Generator over all k-mers in seq.
    There are [len(seq) - k + 1] k-mers in seq.

    Parameters
    ----------
    seq : str
        Sequence which will be broken into kmers.
    k : int
        Length of peptides to return.
    gapped : bool
        If True (default), yield the k-mer including gaps.
        If False,  yield the "non-gapped" k-mer from grabKmer
    return_inds : bool
        If True, also yield an array of indices from grabKmerInds

    Yields
    ------
    mer : str
        If gapped, then a k-length peptide starting at starti from seq.
        If seq[starti] is a gap then returns None.
        If not gapped then all gaps are removed before taking the k-length peptide
            (if there aren't k AAs then return is None)
    inds : nd.array (optional)
        An array of indices for the mer"""

    for i in range(len(seq) - k + 1):
        g,ng = grabKmer(seq, i, k=k)
        if gapped:
            mer = g
        else:
            mer = ng
        if yield_inds:
            ginds,nginds = grabKmerInds(seq, i, k=k)
            if gapped:
                inds = ginds
            else:
                inds = nginds
            yield (mer,inds)
        else:
            yield (mer,)

def grabKmer(seq, starti, k=9):
    """Grab the kmer from seq starting at position starti with length k
    Return the gapped and non-gapped kmer

    If seq[starti] is a gap then the non-gapped kmer is None.
    If there are not enough non-gap AA to return after starti then it returns None

    Parameters
    ----------
    seq : str
        Sequence from which peptide will be grabbed.
    starti : int
        Starting position of the kmer (zero-based indexing)
    k : int
        Length of the peptide to return.

    Returns
    -------
    gapped : str
        A k-length peptide starting at starti from seq.
    nonGapped : str
        A k-length peptide starting at starti from seq.
        If seq[starti] is a gap then returns None.
        If not then all gaps are removed before taking the k-length peptide
            (if there aren't k AAs then return is None)"""
    if not isinstance(starti,int):
        starti = int(starti)

    if (starti+k-1) <= (len(seq)-1) and starti >= 0:
        tmp = seq[starti:]
        full = tmp[:k]
        if full[0] == '-':
            return None,None
        elif '-' in full:
            ng = tmp.replace('-','')
            if len(ng) >= k:
                ng = ng[:k]
            else:
                ng = None
        else:
            ng = full
        return full,ng
    else:
        return None,None

def grabKmerInds(seq, starti, k = 9):
    """Grab the kmer from seq starting at position starti with length k
    Return the indices of the gapped and non-gapped kmers

    i.e. indices are such that seq[ind] == kmer

    If seq[starti] is a gap then the non-gapped kmer is None.
    If there are not enough non-gap AA to return after starti then it returns None

    Parameters
    ----------
    seq : str
        Sequence from which peptide will be grabbed.
    starti : int
        Starting position of the kmer (zero-based indexing)
    k : int
        Length of the peptide to return.

    Returns
    -------
    gapped : ndarray
        A k-length vector starting with starti containing the indices for the kmer
    nonGapped : ndarray
        A k-length vector starting at starti.
        If seq[starti] is a gap then returns None.
        If not then all gaps are removed before taking the k-length peptide
            (if there aren't k AAs then return is None)"""
    if not isinstance(starti,int):
        starti = int(starti)

    if (starti+k-1) <= (len(seq)-1) and starti >= 0:
        tmp = np.arange(starti,len(seq))
        full = tmp[:k]
        """If it starts with a gap then it is invalid (arbitary rule)"""
        if seq[starti] == '-':
            return None,None
        elif '-' in seq[starti:starti+k]:
            """If there's a gap somewhere else then go through one by one adding non-gapped indices"""
            ng = []
            for sitei in tmp:
                if not seq[sitei] == '-':
                    ng.append(sitei)
                """If we get to k non-gapped AAs then return full,ng"""
                if len(ng) == k:
                    return full,ng
            """If we get to then end of the seq then return ng=None"""
            return full,None
        else:
            """If there are no gaps anywhere then just return k indices starting with starti"""
            return full,full
    else:
        """If its an invalid request then return None,None"""
        return None,None

def findpeptide(pep, seq, returnEnd = False):
    """Find pep in seq ignoring gaps but returning a start position that counts gaps
    pep must match seq exactly (otherwise you should be using pairwise alignment)

    Parameters
    ----------
    pep : str
        Peptide to be found in seq.
    seq : str
        Sequence to be searched.
    returnEnd : bool
        Flag to return the end position such that:
        seq[startPos:endPos] = pep

    Returns
    -------
    startPos : int
        Start position (zero-indexed) of pep in seq or -1 if not found"""

    ng = seq.replace('-','')
    ngInd = ng.find(pep)
    ngCount = 0
    pos = 0
    """Count the number of gaps prior to the non-gapped position. Add them to it to get the gapped position"""
    while ngCount < ngInd or seq[pos] == '-':
        if not seq[pos] == '-':
            ngCount += 1
        pos += 1
    startPos = ngInd + (pos - ngCount)

    if returnEnd:
        if startPos == -1:
            endPos = -1
        else:
            count = 0
            endPos = startPos
            while count < len(pep):
                if not seq[endPos] == '-':
                    count += 1
                endPos += 1
        return startPos,endPos
    else:
        return startPos

def grabOverlappingKmer(seq,sitei, pos = 0, k = 9):
    """Grab the kmer from seq for which it is in the pos position at sitei
    Return the gapped and non-gapped kmer

    This is a generalization of grabKmer for pos = 0

    If seq[sitei] is a gap then the non-gapped kmer is None.
    If there are not enough non-gap AA to return before/after sitei then it returns None

    Parameters
    ----------
    seq : str
        Sequence from which peptide will be grabbed.
    sitei : int
        Key position of the kmer (zero-based indexing)
    pos : int
        The position of the key sitei in the kmer.
    k : int
        Length of the peptide to return.

    Returns
    -------
    gapped : str
        A k-length peptide that overlaps sitei
    nonGapped : str
        A k-length peptide that overlaps sitei
        If seq[sitei] is a gap then returns None.
        If not then all gaps are removed before taking the k-length peptide
            (if there aren't k AAs then return is None)"""
    aaRight = k - pos
    aaLeft = pos
    if seq[sitei] == '-':
        return None,None

    if (sitei + aaRight) <= len(seq) and (sitei - aaLeft) >= 0:
        if pos<k:
            rh = seq[sitei:]
            fullRH = rh[:aaRight]
            if '-' in fullRH:
                ngRH = rh.replace('-','')
                if len(ngRH) >= aaRight:
                    ngRH = ngRH[:aaRight]
                else:
                    ngRH = None
            else:
                ngRH = fullRH
        else:
            fullRH = ''
            ngRH = ''

        if pos>0:
            lh = seq[:sitei]
            fullLH = lh[-aaLeft:]
            if '-' in fullLH:
                ngLH = lh.replace('-','')
                if len(ngLH) >= aaLeft:
                    ngLH = ngLH[-aaLeft:]
                else:
                    ngLH = None
            else:
                ngLH = fullLH
        else:
            fullLH = ''
            ngLH = ''
        full = fullLH + fullRH
        #print aaLeft,fullLH,",", aaRight,fullRH

        if ngLH is None or ngRH is None:
            ng = None
        else:
            ng = ngLH + ngRH
        return full,ng
    else:
        return None,None


def overlappingMers(seq, sitei, nmer = [8, 9, 10, 11], padding = 0):
    """Create a list of kmers that overlap sitei in seq

    Returns parallel lists of the mers, start positions and lengths

    Parameters
    ----------
    seq : str
    sitei : int
        Zero-based index into seq
    nmer : list
        Lengths of kmers to consider
    padding : int
        Allow kmer to be within padding.
        Defalut is no padding (must overlap)

    Returns
    -------
    mers : list
        List of overlapping peptides
    starti : list
        List of start positions"""
    
    def _overlappingMersNoPadding(seq, sitei, nmer):
        mers = []
        starti = []
        for k in nmer:
            for posi in range(k):
                ng = grabOverlappingKmer(seq, sitei, pos=posi, k=k)[1]
                if not ng is None:
                    mers.append(ng)
                    starti.append(sitei-posi)
                    #print sitei, posi, k, ng
        mers,uniqi = np.unique(mers,return_index = True)
        starti = np.array(starti)[uniqi]
        return mers, starti

    mers,starti = _overlappingMersNoPadding(seq, sitei, nmer = nmer)
    if padding > 0:
        for padi in (np.arange(padding) + 1):
            for tmpSitei in [sitei+padi, sitei-padi]:
                tmpMers, tmpStarti = _overlappingMersNoPadding(seq, tmpSitei, nmer)
                mers = np.concatenate((mers,tmpMers))
                starti = np.concatenate((starti,tmpStarti))
    mers,uniqi = np.unique(mers,return_index = True)
    starti = np.array(starti)[uniqi]
    return mers,starti