from pylab import *

from custom_legends import colorLegend
from objhist import *

__all__ = ['plotHLAFreq']

ic50Ticks=array([10, 25, 50, 100, 250, 500, 1000, 2500, 5000, 10e3], dtype=int)

def plotHLAFreq(hlaDf,twoDigit=True,groupColumn=None,stackColumn=None,loci='AB',cutoff=0.05,annotateCount=True):
    """Bar plot of HLA allele expression (optionally by group)
    Requires "hlas" column in hlaDf which is a list of alleles"""
    allHLAs=[]
    for hlas in hlaDf.hlas:
        allHLAs.extend(hlas)
    if twoDigit:
        allHLAs=[h[:4] for h in allHLAs]
    uHLA=list(set(allHLAs))

    if not groupColumn is None:
        uGroup=sorted(hlaDf[groupColumn].unique().tolist())
        groupS=hlaDf[groupColumn]
    else:
        uGroup=['All']
        groupS=pd.Series(data=[uGroup[0] for i in arange(hlaDf.shape[0])], index=hlaDf.index, name='Group')

    if not stackColumn is None:
        uStack=sorted(hlaDf[stackColumn].unique().tolist())
        stackS=hlaDf[stackColumn]
    else:
        uStack=['All']
        stackS=pd.Series(data=[uStack[0] for i in arange(hlaDf.shape[0])], index=hlaDf.index, name='Stack')

    counts=zeros((len(uGroup), len(uStack), len(uHLA)))

    for hlas, group, stack in zip(hlaDf.hlas, groupS, stackS):
        if twoDigit:
            tmp={h[:4] for h in hlas}
        else:
            tmp=set(hlas)
        for h in tmp:
            """For each person's unique alleles, increment the counter"""
            hlai=uHLA.index(h)
            groupi=uGroup.index(group)
            stacki=uStack.index(stack)
            counts[groupi, stacki, hlai]+=1
    sHLA=sorted(uHLA, key=lambda h: counts[:,:, uHLA.index(h)].sum(), reverse=True)
    mycm=array(['darkgreen', 'fuchsia', 'saddlebrown', 'lightseagreen', 'gold', 'royalblue', 'tomato', 'thistle', 'tan'])
    myalphas=linspace(1, 0.3, len(uStack))
    N=hlaDf.shape[0]
    groupN=objhist(groupS)
    w=1/(len(uGroup)+2)

    clf()
    for locusi, locus in enumerate(loci):
        subplot(len(loci), 1, locusi+1)
        xlabs=[]
        x=0
        other=zeros((len(uGroup), len(uStack)))
        for hla in sHLA:
            hlai=uHLA.index(hla)
            if hla[0]==locus:
                if counts[:,:, hlai].sum()/N > cutoff:
                    xlabs.append(hla)
                    for groupi, group in enumerate(uGroup):
                        for stacki, stack in enumerate(uStack):
                            bar(left=x+w*groupi-len(uGroup)*w/2,
                                height=1e2*counts[groupi, stacki, hlai]/groupN[group],
                                bottom=1e2*counts[groupi, :stacki, hlai].sum()/groupN[group],
                                width=w, color=mycm[groupi], alpha=myalphas[stacki])
                            if counts[groupi, stacki, hlai]>1 and annotateCount:
                                annotate('%d' % counts[groupi, stacki, hlai],
                                        xy=(x+w*groupi-len(uGroup)*w/2+w/2, 1e2*counts[groupi, :(stacki+1), hlai].sum()/groupN[group]),
                                        ha='center', va='top', xytext=(0, -4), textcoords='offset points')
                    x+=1
                else:
                    for groupi, group in enumerate(uGroup):
                        for stacki, stack in enumerate(uStack):
                            other[groupi, stacki]+=counts[groupi, stacki, hlai]
        xlabs.append('Other')
        for groupi, group in enumerate(uGroup):
            for stacki, stack in enumerate(uStack):
                bar(left=x+w*groupi,
                    height=1e2*other[groupi, stacki]/groupN[group],
                    bottom=1e2*other[groupi, :stacki].sum()/groupN[group],
                    width=w, color=mycm[groupi], alpha=myalphas[stacki])

        if len(uGroup)>1:
            sorti=argsort(uGroup)
            colorLegend(colors=mycm[:len(uGroup)][sorti], labels=[uGroup[i] for i in sorti], title=groupColumn)
        ylabel('% of participants', size='x-large')
        xticks(arange(len(xlabs)), xlabs, size='x-large')
        xlim((-0.5, len(xlabs)+1))


