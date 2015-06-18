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
            'findpeptide']

BADAA='-*BX#Z'
AALPHABET='ACDEFGHIKLMNPQRSTVWY'


def convertHLAAsterisk(hlas):
    """Replace the * with _ in each HLA allele"""
    repAsteriskPattern = re.compile('\*')
    return [re.sub(repAsteriskPattern,'_',h) for h in hlas]

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

def rankEpitopes(ba,hlaList,peptide,nmer=[8,9,10,11],peptideLength=None):
    """Breaks the mer into kmers and then ranks all the (hla,kmer) pairs
    Returns (ranks,sorti,kmers,ic50,hla)"""
    merList = getMers(peptide,nmer,peptideLength)
    kmers = empty((len(merList),len(hlaList)),dtype=object)
    ic50 = ones((len(merList),len(hlaList)))*15
    hla = empty((len(merList),len(hlaList)),dtype=object)
    for i,m in enumerate(merList):
        for j,h in enumerate(hlaList):
            kmers[i,j]=m
            hla[i,j]=h
            tmp=ba[(h,m)]
            if not isnan(tmp):
                ic50[i,j]=tmp
    kmers=kmers.flatten()
    ic50=ic50.flatten()
    hla=hla.flatten()
    sorti = ic50.argsort()
    ranks = empty(len(ic50), int)
    ranks[sorti] = arange(len(ic50))
    return (ranks,sorti,kmers,ic50,hla)

def rankKmers(ba,hlaList,peptide,nmer=[8,9,10,11],peptideLength=16):
    """Breaks the mer into kmers and then ranks the kmers"""
    kmers=getMers(peptide,nmer,peptideLength)
    result=rankMers(ba,hlaList,kmers)
    return (result[0],result[1],kmers,result[2],result[3])
    
def rankMers(ba,hlaList,merList):
    """
    Rank the peptides in merList based on the min. IC50 across HLAs in hlaList
    Top rank is 0
    Peptides not in BA get a 15
    """
    ic50 = np.ones((len(merList)))*15
    hla = np.empty(len(merList),dtype=object)
    for i,m in enumerate(merList):
        if not '.' in m:
            ic50[i],hla[i],dumpmer=getIC50(ba,hlaList,m)
    sorti = ic50.argsort()
    ranks = np.empty(len(ic50), int)
    ranks[sorti] = np.arange(len(ic50))
    return (ranks,sorti,ic50,hla)

def getIC50(ba,hlaList,mer,nmers=[8,9,10,11]):
    """Return the IC50 from ba of the mer and its affinity with the most avid HLA in hlaList.
    Or if len(pep)>11, return that of the most avid kmer"""
    if ba is None:
        raise NameError('Did not load IC50 values into ba!')
    #minimum IC50 over the HLAs
    if len(mer)<=11:
        allPairs=[(ba[(h,mer)],h,mer) for h in hlaList]
    #minimum IC50 over all the mers and all the HLAs
    else:
        allPairs=[getIC50(ba,hlaList,m) for m in getMers(mer,nmers)]
    return min(allPairs,key=lambda x: x[0])


def getMers(seq,nmers=[8, 9 , 10, 11],seqLength=None):
    """Takes a AA sequence (string) and turns it into a list of 8, 9, 10, 11 mers 
    seq will be padded with '.' if shorter than seqLength"""
    if not seqLength is None:
        if len(seq) > seqLength:
            seq=seq[:seqLength]
        elif len(seq) < seqLength:
            seq=string.ljust(seq,16,'.')

    mers=[]
    for n in nmers:
        mers.extend([seq[i:i+n] for i in range(len(seq)-n+1)])
    return mers

def getMerInds(seq,nmers=[8, 9 , 10, 11],seqLength=None):
    """Takes a AA sequence (string) and turns it into a list of 8, 9, 10, 11 mers
    seq will be padded with '.' if shorter than seqLength"""
    if not seqLength is None:
        if len(seq) > seqLength:
            seq=seq[:seqLength]
        elif len(seq) < seqLength:
            seq=string.ljust(seq,16,'.')

    mers=[]
    inds=[]
    for n in nmers:
        mers.extend([seq[i:i+n] for i in range(len(seq)-n+1)])
        inds.extend([np.arange(n)+i for i in range(len(seq)-n+1)])
    return mers,inds

def grabKmer(seq,starti,k=9):
    """Grab the kmer from seq starting at position starti with length k
       Return the gapped and non-gapped kmer"""
    if (starti+k-1) <= (len(seq)-1) and starti>=0:
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
def grabKmerInds(seq,starti,k=9):
    """Grab the kmer from seq starting at position starti with length k
       Return the indices of the gapped and non-gapped kmer"""
    if (starti+k-1)<=(len(seq)-1) and starti>=0:
        tmp=np.arange(starti,len(seq))
        full=tmp[:k]
        """If it starts with a gap then it is invalid (arbitary rule)"""
        if seq[starti]=='-':
            return None,None
        elif '-' in seq[starti:starti+k]:
            """If there's a gap somewhere else then go through one by one adding non-gapped indices"""
            ng=[]
            for sitei in tmp:
                if not seq[sitei]=='-':
                    ng.append(sitei)
                """If we get to k non-gapped AAs then return full,ng"""
                if len(ng)==k:
                    return full,ng
            """If we get to then end of the seq then return ng=None"""
            return full,None
        else:
            """If there are no gaps anywhere then just return k indices starting with starti"""
            return full,full
    else:
        """If its an invalid request then return None,None"""
        return None,None

def findpeptide(pep,seq,returnEnd=False):
    """Find pep in seq ignoring gaps but returning a start position that counts gaps
       pep must match seq exactly (otherwise you should be using pairwise alignment)
       Returns startPos of pep in seq (zero-indexed) or -1 if not found"""
    ng=seq.replace('-','')
    ngInd=ng.find(pep)
    ngCount=0
    pos=0
    """Count the number of gaps prior to the non-gapped position. Add them to it to get the gapped position"""
    while ngCount<ngInd or seq[pos]=='-':
        if not seq[pos]=='-':
            ngCount+=1
        pos+=1
    startPos=ngInd+(pos-ngCount)

    if returnEnd:
        if startPos==-1:
            endPos=-1
        else:
            count=0
            endPos=startPos
            while count<len(pep):
                if not seq[endPos]=='-':
                    count+=1
                endPos+=1
        return startPos,endPos
    else:
        return startPos
