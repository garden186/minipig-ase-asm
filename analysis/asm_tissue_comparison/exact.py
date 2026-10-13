"""Exact equal-OR test on full fixed-margin support, with no sign conditioning."""
from functools import lru_cache
from math import comb, log


def checked(table):
    t=tuple(int(x) for x in table)
    if len(t)!=4 or any(x<0 for x in t) or tuple(table)!=t or min(t[0]+t[1],t[2]+t[3])<=0:
        raise ValueError('Four integer counts and positive allele depths required')
    return t


def sign(table):
    a,b,c,d=map(int,table);cross=a*(c+d)-c*(a+b)
    return (cross>0)-(cross<0)


def paired_difference(first,second):
    """Exact rational fraction contrast and sign, including exact ties."""
    a,b,c,d=map(int,first);e,f,g,h=map(int,second)
    n1=a*(c+d)-c*(a+b);d1=(a+b)*(c+d)
    n2=e*(g+h)-g*(e+f);d2=(e+f)*(g+h)
    numerator=n1*d2-n2*d1
    return numerator/(d1*d2),(numerator>0)-(numerator<0)


@lru_cache(maxsize=32768)
def pair_test(first,second,max_states=200000):
    a,b,c,d=checked(first);e,f,g,h=checked(second)
    na,nr,m=a+b,c+d,a+c;nb,ns,k=e+f,g+h,e+g;total=a+e
    l1,h1=max(0,m-nr),min(na,m);l2,h2=max(0,k-ns),min(nb,k)
    lo,hi=max(l1,total-h2),min(h1,total-l2)
    if hi-lo+1>max_states:raise ValueError('Exact support limit exceeded; no approximate P substituted')
    if not lo<=a<=hi:raise ValueError('Observed table absent from full conditional support')
    observed=comb(na,a)*comb(nr,c)*comb(nb,e)*comb(ns,g)
    y=total-lo
    weight=comb(na,lo)*comb(nr,m-lo)*comb(nb,y)*comb(ns,k-y)
    numerator=denominator=0
    for x in range(lo,hi+1):
        denominator+=weight
        if weight<=observed:numerator+=weight
        if x<hi:
            y=total-x
            top=(na-x)*(m-x)*y*(ns-k+y)
            bottom=(x+1)*(nr-m+x+1)*(nb-y+1)*(k-y+1)
            weight,remainder=divmod(weight*top,bottom)
            if remainder:raise ArithmeticError('Exact integer weight recurrence has remainder')
    return min(0.,log(numerator)-log(denominator)),hi-lo+1


def reference(first,second):
    """Independent direct-combination oracle; no recurrence or bound helper."""
    a,b,c,d=map(int,first);e,f,g,h=map(int,second)
    total=a+e;na,nr,m=a+b,c+d,a+c;nb,ns,k=e+f,g+h,e+g
    observed=comb(na,a)*comb(nr,c)*comb(nb,e)*comb(ns,g)
    weights=[]
    for x in range(na+1):
        y=total-x
        if not (0<=m-x<=nr and 0<=y<=nb and 0<=k-y<=ns):continue
        weights.append((x,comb(na,x)*comb(nr,m-x)*comb(nb,y)*comb(ns,k-y)))
    return dict(numerator=sum(w for _,w in weights if w<=observed),denominator=sum(w for _,w in weights),full_states=len(weights),weights=weights)
