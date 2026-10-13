"""P-independent common-SNP ranking, reusing phase patterns across CpGs."""
import numpy as np
from common import require

POPCOUNT=np.array([i.bit_count() for i in range(1024)],dtype='u1')
MULTI=2**32

def variant_keys(a):return a['pos'].astype(np.uint64)*16+a['ref'].astype(np.uint64)*4+a['alt']

def pool(parts):
    """Merge variant supports; a fragment/tissue never creates a new animal."""
    if not parts:return np.empty(0,dtype='<u8'),np.empty(0,dtype='<u2')
    k,ix=np.unique(np.concatenate([p[0] for p in parts]),return_inverse=True)
    masks=np.zeros(len(k),dtype='<u2')
    np.bitwise_or.at(masks,ix,np.concatenate([np.full(len(a),b,dtype='<u2') for a,b in parts]))
    return k,masks

def nearest(winners,positions):
    if not len(winners):return np.zeros(len(positions),dtype='<u8')
    # All tied variants at the same position share distance; use the first key.
    winners=np.sort(winners);pos=winners//16;winners=winners[np.r_[True,pos[1:]!=pos[:-1]]];pos=winners//16
    j=np.searchsorted(pos,positions);left=np.maximum(0,j-1);right=np.minimum(len(pos)-1,j)
    dl=np.abs(pos[left].astype(np.int64)-positions.astype(np.int64));dr=np.abs(pos[right].astype(np.int64)-positions.astype(np.int64))
    return winners[np.where(dl<=dr,left,right)]

def assign(positions,available_pairs,tested_pairs,phase_arrays,guard=lambda:None):
    """Pairs encode site_index<<32 | phase_dictionary_index, sorted and unique.

    Only available/tested PS membership is accepted; no counts or P enter here.
    The common one-PS case shares variant pools and then ranks each test mask.
    Rare multiple-PS sites use the exact same ranking without pooling counts.
    """
    n=len(positions);animals=len(phase_arrays);require(animals<=10,'Too many animals')
    rawsig=np.zeros((n,animals),dtype='<u8');tmask=np.zeros(n,dtype='<u2');multi={};testsets={}
    for i,(av,te) in enumerate(zip(available_pairs,tested_pairs)):
        require(np.all(np.isin(te,av,assume_unique=True)),'Tested PS missing from available opportunity')
        if not len(av):continue
        index=(av>>32).astype(np.int64);code=(av&0xffffffff).astype('<u4')
        require(np.all(code>0) and np.all(index<n),'Bad encoded opportunity')
        sites,start,count=np.unique(index,return_index=True,return_counts=True)
        one=count==1;rawsig[sites[one],i]=code[start[one]]
        for site,first,length in zip(sites[~one],start[~one],count[~one]):
            rawsig[site,i]=MULTI+int(site);multi[(int(site),i)]=tuple(map(int,code[first:first+length]))
        if len(te):
            ti=(te>>32).astype(np.int64);tc=(te&0xffffffff).astype('<u4');tmask[np.unique(ti)]|=1<<i
            for site in sites[~one]:
                lo,hi=np.searchsorted(ti,site,'left'),np.searchsorted(ti,site,'right');testsets[(int(site),i)]=set(map(int,tc[lo:hi]))
    sig,inv=np.unique(rawsig,axis=0,return_inverse=True);order=np.argsort(inv,kind='stable');starts=np.r_[0,np.cumsum(np.bincount(inv,minlength=len(sig)))]
    anchors=np.zeros(n,dtype='<u8');nt=np.zeros(n,dtype='u1');na=nt.copy();nc=np.zeros(n,dtype='<u4')
    cache={}
    def pkeys(i,ps):
        key=(i,ps)
        if key not in cache:
            a=phase_arrays[i];lo=np.searchsorted(a['ps'],ps,'left');hi=np.searchsorted(a['ps'],ps,'right')
            cache[key]=variant_keys(a[lo:hi])
        return cache[key]
    for z,row in enumerate(sig):
        if z%256==0:guard()
        ix=order[starts[z]:starts[z+1]]
        if np.any(row>=MULTI):
            require(len(ix)==1,'Multiple-PS signature must be site-specific');site=int(ix[0]);apart=[];tpart=[]
            for i,ps in enumerate(row):
                if ps==0:continue
                codes=multi[(site,i)] if ps>=MULTI else [int(ps)]
                for code in codes:
                    k=pkeys(i,code)
                    if not len(k):continue
                    apart.append((k,1<<i))
                    is_test=code in testsets.get((site,i),set()) if ps>=MULTI else bool(tmask[site]&(1<<i))
                    if is_test:tpart.append((k,1<<i))
            k,am=pool(apart);tk,tm=pool(tpart)
            if not len(k):continue
            test=np.zeros(len(k),dtype='<u2');test[np.searchsorted(k,tk)]=tm
            score=POPCOUNT[test].astype(int)*11+POPCOUNT[am];best=int(score.max());win=k[score==best]
            anchors[ix]=nearest(win,positions[ix]);nt[ix]=best//11;na[ix]=best%11;nc[ix]=len(k)
        else:
            parts=[(pkeys(i,int(ps)),1<<i) for i,ps in enumerate(row) if ps]
            k,am=pool([(a,b) for a,b in parts if len(a)])
            if not len(k):continue
            masks,mi=np.unique(am,return_inverse=True)
            for test in np.unique(tmask[ix]):
                where=ix[tmask[ix]==test];score=POPCOUNT[masks&test].astype(int)*11+POPCOUNT[masks]
                best=int(score.max());win=k[(score==best)[mi]]
                anchors[where]=nearest(win,positions[where]);nt[where]=best//11;na[where]=best%11;nc[where]=len(k)
    return anchors,nt,na,nc,dict(phase_patterns=len(sig),sites=n,multiple_PS_sites=int(np.any(rawsig>=MULTI,axis=1).sum()))
