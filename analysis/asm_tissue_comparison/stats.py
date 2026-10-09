"""Directional exact tails, Fisher partial conjunction, full-family log BH."""
import math
import numpy as np
from scipy.stats import hypergeom
from scipy.special import gammaln, logsumexp
from common import require

LOG_MIN=math.log(np.nextafter(0.,1.))
LOG2=math.log(2.)

def exact_log_tail(a,b,c,d,upper):
    """Finite hypergeometric sum in bounded memory; used only for underflow."""
    total=a+b+c+d;success=a+c;draw=a+b
    lo=max(0,draw-(total-success));hi=min(success,draw)
    lo,hi=(a,hi) if upper else (lo,a)
    answer=-np.inf
    constant=-(gammaln(total+1)-gammaln(draw+1)-gammaln(total-draw+1))
    for first in range(lo,hi+1,65536):
        x=np.arange(first,min(hi+1,first+65536),dtype=float)
        lp=(gammaln(success+1)-gammaln(x+1)-gammaln(success-x+1)
            +gammaln(total-success+1)-gammaln(draw-x+1)-gammaln(total-success-draw+x+1)+constant)
        answer=np.logaddexp(answer,logsumexp(lp))
    require(np.isfinite(answer),'Non-finite exact log tail')
    return min(0.,float(answer))

def directional_logp(counts):
    """Rows ALT M,U / REF M,U. Columns ALT>REF, REF>ALT; no two-sided P/2."""
    counts=np.asarray(counts)
    require(counts.ndim==2 and counts.shape[1]==4 and np.issubdtype(counts.dtype,np.integer) and np.all(counts>=0),'Invalid integer table')
    if not len(counts):return np.empty((0,2)),0
    a,b,c,d=counts.astype(np.int64).T;total=a+b+c+d;success=a+c;draw=a+b
    require(np.all((draw>0)&(c+d>0)),'Zero haplotype depth')
    pp=np.column_stack([hypergeom.sf(a-1,total,success,draw),hypergeom.cdf(a,total,success,draw)])
    require(np.all(np.isfinite(pp)&(pp>=0)&(pp<=1)),'Invalid hypergeometric tail')
    logs=np.zeros_like(pp);np.log(pp,out=logs,where=pp>0)
    fallback=np.argwhere(pp==0)
    for i,j in fallback:logs[i,j]=exact_log_tail(*map(int,counts[i]),upper=j==0)
    require(np.all(np.isfinite(logs)&(logs<=0)),'Invalid log tails')
    return logs,len(fallback)

def gamma_logsf(s,k):
    s=np.asarray(s,dtype=float);require(isinstance(k,int) and k>=1 and np.all(np.isfinite(s)&(s>=0)),'Invalid gamma tail')
    l=np.full(s.shape,-np.inf);np.log(s,out=l,where=s>0);term=np.zeros(s.shape);total=term.copy()
    for j in range(1,k):term=term+l-math.log(j);total=np.logaddexp(total,term)
    return np.minimum(0.,-s+total)

def partial_conjunction(logp,r):
    p=np.asarray(logp,dtype=float)
    require(p.ndim==2 and 1<=r<=p.shape[1] and np.all(np.isfinite(p)&(p<=0)),'Invalid PC input / N / r')
    chosen=np.sort(p,axis=1)[:,r-1:]
    return gamma_logsf(-chosen.sum(axis=1),p.shape[1]-r+1)

def bh_log(logp):
    p=np.asarray(logp,dtype=float);require(p.ndim==1 and np.all(np.isfinite(p)&(p<=0)),'Invalid complete BH family')
    if not len(p):return p.copy()
    order=np.argsort(p,kind='stable');value=np.minimum(0.,p[order]+math.log(len(p))-np.log(np.arange(1,len(p)+1)))
    value=np.minimum.accumulate(value[::-1])[::-1];out=np.empty_like(p);out[order]=value
    require(np.all(out>=p-1e-11),'BH q below P')
    return out

def adjusted_direction_logq(raw_logq):return np.minimum(0.,np.asarray(raw_logq)+LOG2)
def display(logp):return np.exp(np.maximum(LOG_MIN,np.asarray(logp)))
