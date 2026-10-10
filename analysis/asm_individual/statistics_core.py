"""Exact small-table Fisher arithmetic, bounded memoization, and ordinary BH."""
from functools import lru_cache
import math
import numpy as np
from scipy.stats import fisher_exact


def canonical(a,b,c,d):
    return min((a,b,c,d),(b,a,d,c),(c,d,a,b),(d,c,b,a),
               (a,c,b,d),(c,a,d,b),(b,d,a,c),(d,b,c,a))


@lru_cache(maxsize=200000)
def _fisher(t):
    a,b,c,d=t; n1=a+b;n2=c+d;total=n1+n2;col=a+c
    if total<=200:
        lo,hi=max(0,n1-(total-col)),min(n1,col)
        observed=math.comb(col,a)*math.comb(total-col,n1-a)
        numerator=sum(w for x in range(lo,hi+1)
            if (w:=math.comb(col,x)*math.comb(total-col,n1-x))<=observed)
        p=numerator/math.comb(total,n1)
    else:
        p=float(fisher_exact([[a,b],[c,d]],alternative='two-sided').pvalue)
    if not math.isfinite(p) or not 0<=p<=1: raise ValueError('Invalid Fisher P value')
    return max(math.nextafter(0.0,1.0),p),int(p==0)


def fisher(x):
    if len(x)!=4 or any(type(v) is not int or v<0 for v in x): raise ValueError('Fisher requires four nonnegative integers')
    a,b,c,d=x
    if not a+b or not c+d: raise ValueError('Both haplotypes need valid observations')
    if not a+c or not b+d: return 1.0,0
    return _fisher(canonical(a,b,c,d))


def bh(values):
    p=np.asarray(values,dtype=np.float64)
    if p.ndim!=1 or np.any(~np.isfinite(p)) or np.any((p<0)|(p>1)): raise ValueError('Invalid BH P values')
    m=p.size
    if not m:return p.copy()
    order=np.argsort(p,kind='stable')
    adjusted=np.minimum.accumulate((p[order]*(m/np.arange(1,m+1,dtype=np.float64)))[::-1])[::-1]
    out=np.empty_like(p);out[order]=np.minimum(adjusted,1.0)
    return out
