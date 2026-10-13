"""Whole aligned paired input; observed signs are descriptive, never filters."""
import csv,gzip
import numpy as np
from common import require,IDS,TISSUES
from schema import PROFILE,SITE

PAIRS=[(0,1),(0,2),(1,2)]
PAIR_NAMES=['Heart_Kidney','Heart_Liver','Kidney_Liver']
FILTER=np.dtype([('pos','<u4'),('N','u1',(3,)),('positive','u1',(3,)),('negative','u1',(3,)),
    ('ties','u1',(3,)),('paired_N','u1',(3,)),('reason','u1',(3,))])
REASONS=['PASS','PAIRED_N_BELOW_2']
RESULT=np.dtype([('pos','<u4'),('pair','u1'),('N','u1'),('N_first','u1'),('N_second','u1'),
    ('sign_counts_first','u1',(3,)),('sign_counts_second','u1',(3,)),('D_sign_counts','u1',(3,)),('gene','<u4'),('anchor','<u8'),
    ('mean_D','<f8'),('median_D','<f8'),('min_D','<f8'),('max_D','<f8'),
    ('median_first_all','<f8'),('median_second_all','<f8'),('median_first_paired','<f8'),('median_second_paired','<f8'),
    ('singleton_animals','u1'),('log_pc','<f8',(9,))])
INDIVIDUAL=np.dtype([('logp','<f8',(10,)),('full_states','<u4',(10,)),('D','<f8',(10,)),('D_sign','i1',(10,))])
FAMILY=np.dtype([('chrom','u1'),('pos','<u4'),('tile','<u2'),('row','<u4'),('N','u1'),
    ('logp','<f8'),('logq_BH','<f8'),('logq_BY','<f8')])

def rows(path):
    op=gzip.open if str(path).endswith('.gz') else open
    with op(path,'rt',encoding='utf-8-sig',newline='') as f:yield from csv.DictReader(f,delimiter='\t')

def table(path,data,fields):
    op=gzip.open if str(path).endswith('.gz') else open
    with op(path,'wt',encoding='utf-8',newline='') as f:
        w=csv.DictWriter(f,fields,delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(data)

def orient(pr,sites):
    require(pr.dtype==PROFILE and sites.dtype==SITE,'Stage2 schema changed')
    require(not len(sites) or np.all(sites['pos'][1:]>sites['pos'][:-1]),'Duplicate/unsorted sites')
    ix=np.searchsorted(sites['pos'],pr['pos'])
    require(np.all(ix<len(sites)) and np.array_equal(sites['pos'][ix],pr['pos']),'Profile coordinate outside sites')
    require(np.all(pr['sample']<10)&np.all(pr['tissue']<3),'Unknown animal/tissue')
    good=pr['state']==0;r=pr[good];j=ix[good];a=r['sample'].astype(int);t=r['tissue'].astype(int)
    require(len(np.unique((j.astype(np.int64)*10+a)*3+t))==len(j),'Duplicate aligned PS per animal/tissue/CpG; no pooling')
    require(np.all(((r['GT']==1)&(r['alt_hap']==2))|((r['GT']==2)&(r['alt_hap']==1))),'GT/ALT mismatch')
    require(np.all(r['anchor_ps']>0)&np.all(sites['anchor'][j]>0),'Missing aligned anchor/PS')
    counts=np.full((len(sites),10,3,4),-1,dtype='<i4');phase=np.zeros((len(sites),10,3),dtype='<u4')
    links=np.full((len(sites),10,3),-1,dtype='<i4')
    x=np.column_stack([r[k] for k in 'abcd'])
    require(np.all(x<=np.iinfo(np.int32).max),'Count overflow')
    x=x.astype(np.int64);n1=x[:,0]+x[:,1];n2=x[:,2]+x[:,3]
    require(np.all((n1>=3)&(n2>=3)&(n1+n2>=10)),'Original depth eligibility changed')
    # int32 counts imply count cross products fit signed int64 when each total
    # is <= int32.max; refuse unexpected larger totals instead of overflowing.
    require(np.all(n1<=np.iinfo(np.int32).max)&np.all(n2<=np.iinfo(np.int32).max),'Allele total overflow')
    swap=r['alt_hap']==2;x[swap]=x[swap][:,[2,3,0,1]]
    counts[j,a,t]=x;phase[j,a,t]=r['anchor_ps'];links[j,a,t]=np.flatnonzero(good)
    for u,v in PAIRS:
        both=(phase[:,:,u]>0)&(phase[:,:,v]>0)
        require(np.all(phase[:,:,u][both]==phase[:,:,v][both]),'Different exact animal PS across tissues')
    require(np.array_equal((counts[:,:,:,0]>=0).sum(axis=1),sites['N']),'Aligned N differs from stage2')
    return counts,phase,links

def filters(counts,positions):
    """Eligibility uses paired availability only; retain mixed and exact-zero signs."""
    present=counts[:,:,:,0]>=0;x=counts.astype(np.int64)
    cross=x[:,:,:,0]*(x[:,:,:,2]+x[:,:,:,3])-x[:,:,:,2]*(x[:,:,:,0]+x[:,:,:,1])
    signs=np.sign(cross).astype(np.int8);signs[~present]=0
    f=np.zeros(len(counts),dtype=FILTER);f['pos']=positions;f['N']=present.sum(axis=1)
    f['positive']=(present&(signs>0)).sum(axis=1);f['negative']=(present&(signs<0)).sum(axis=1)
    f['ties']=(present&(signs==0)).sum(axis=1)
    for p,(u,v) in enumerate(PAIRS):
        f['paired_N'][:,p]=(present[:,:,u]&present[:,:,v]).sum(axis=1)
        f['reason'][:,p]=(f['paired_N'][:,p]<2).astype('u1')
    return f,signs,present

def deltas(counts):
    x=counts.astype(float);valid=x[:,:,:,0]>=0
    out=np.full(valid.shape,np.nan)
    out[valid]=x[:,:,:,0][valid]/(x[:,:,:,0][valid]+x[:,:,:,1][valid])-x[:,:,:,2][valid]/(x[:,:,:,2][valid]+x[:,:,:,3][valid])
    return out

def locate(positions,segments):
    require(not len(segments) or np.all(segments['start'][1:]>=segments['end'][:-1]),'Overlapping primary annotation segments')
    out=np.zeros(len(positions),dtype=np.uint32)
    if len(segments):
        p=positions.astype(np.int64)-1;j=np.searchsorted(segments['start'],p,side='right')-1
        ok=j>=0;ok[ok]&=p[ok]<segments['end'][j[ok]];out[ok]=segments['gene'][j[ok]]
    return out
