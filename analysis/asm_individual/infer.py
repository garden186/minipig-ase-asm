"""Individual two-sided ASM and complete unit/tissue/cohort BH families.

The manifest lists every GQ20 tile. A reference FAI and exactly ten animals by
three tissues are required; truncated chromosome sets are rejected.
"""
import argparse,csv,gzip,hashlib,json,math
from array import array
from collections import defaultdict
from pathlib import Path
import numpy as np
from eligibility import parse,COUNTS,SOURCE_FIELDS
from statistics_core import fisher,bh

def rows(p):
    with gzip.open(p,'rt',encoding='utf-8',newline='') as h:yield from csv.DictReader(h,delimiter='\t')
def sha(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as f:
        for x in iter(lambda:f.read(4*1024**2),b''):h.update(x)
    return h.hexdigest()
def coverage(tasks,lengths):
    by=defaultdict(list)
    for t in tasks:by[str(t['chrom'])].append((int(t['start0']),int(t['end0'])))
    if set(by)!=set(lengths):raise ValueError('Missing/extra autosomes')
    for chrom,intervals in by.items():
        last=0
        for a,b in sorted(intervals):
            if a!=last or b<=a:raise ValueError('Gap/overlap in tile manifest')
            last=b
        if last!=lengths[chrom]:raise ValueError('Truncated chromosome')
def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--manifest',type=Path,required=True);p.add_argument('--fai',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True);p.add_argument('--bh-memory-gib',type=float,default=16);a=p.parse_args()
    tasks=json.loads(a.manifest.read_text());groups=defaultdict(list)
    for t in tasks:groups[t['sample'],t['tissue']].append(t)
    animals={s for s,t in groups};tissues={t for s,t in groups}
    if len(animals)!=10 or tissues!={'Heart','Kidney','Liver'} or len(groups)!=30:raise ValueError('Final global30 requires the complete 10-by-3 manifest')
    lengths={r.split('\t')[0]:int(r.split('\t')[1]) for r in a.fai.read_text().splitlines() if r.split('\t')[0] in set(map(str,range(1,19)))}
    if len(lengths)!=18:raise ValueError('Reference requires autosomes 1..18')
    for gg in groups.values():coverage(gg,lengths)
    a.output.mkdir(parents=True,exist_ok=False);policy=dict(minimum_each_haplotype=3,minimum_combined_depth=10)
    ps={};source_hashes={}
    for key,gg in sorted(groups.items()):
        values=array('d');last=None
        for t in sorted(gg,key=lambda t:(int(t['chrom']),int(t['start0']))):
            path=Path(t['path']);source_hashes[str(path)]=sha(path)
            gate=path.parent/'validation.json'
            g=json.loads(gate.read_text())
            if g.get('status')!='PASS_GQ20_FRAGMENT_DEPTH' or g.get('errors') or g['sha256'][path.name]!=source_hashes[str(path)]:raise ValueError('Unvalidated count input')
            for row in rows(path):
                k,x,n1,n2,status=parse(row,t,policy);current=(int(t['chrom']),*k)
                if last is not None and current<=last:raise ValueError('Duplicate or unsorted exact-PS rows')
                last=current
                if status=='TESTED':values.append(fisher(x)[0])
        ps[key]=np.asarray(values,dtype=float)
        print('Tested',key,len(values),flush=True)
    keys=sorted(ps);total=sum(len(v) for v in ps.values())
    if total*100>a.bh_memory_gib*1024**3:raise MemoryError('Complete-family BH forecast exceeds configured memory')
    q_all=bh(np.concatenate([ps[k] for k in keys]));offset=0;globalq={}
    for k in keys:globalq[k]=q_all[offset:offset+len(ps[k])];offset+=len(ps[k])
    qt={}
    for t in sorted(tissues):
        kk=[k for k in keys if k[1]==t];q=bh(np.concatenate([ps[k] for k in kk]));offset=0
        for k in kk:qt[k]=q[offset:offset+len(ps[k])];offset+=len(ps[k])
    extra=['test_status','fisher_two_sided_p','BH_q_within_animal','BH_q_tissue10','BH_q_all30','H1_methylation','H2_methylation','delta_M']
    for key in keys:
        qw=bh(ps[key]);i=0;out=a.output/(key[0]+'_'+key[1]+'.tsv.gz')
        with gzip.open(out,'wt',encoding='utf-8',newline='') as h:
            w=csv.DictWriter(h,fieldnames=SOURCE_FIELDS+extra,delimiter='\t',lineterminator='\n');w.writeheader()
            for t in sorted(groups[key],key=lambda t:(int(t['chrom']),int(t['start0']))):
                if sha(t['path'])!=source_hashes[t['path']]:raise ValueError('Input changed')
                for row in rows(t['path']):
                    _,x,n1,n2,status=parse(row,t,policy)
                    row.update({k:'.' for k in extra});row['test_status']=status
                    row.update(H1_methylation=x[0]/n1 if n1 else '.',H2_methylation=x[2]/n2 if n2 else '.',delta_M=x[0]/n1-x[2]/n2 if n1 and n2 else '.')
                    if status=='TESTED':
                        row.update(fisher_two_sided_p=ps[key][i],BH_q_within_animal=qw[i],BH_q_tissue10=qt[key][i],BH_q_all30=globalq[key][i]);i+=1
                    w.writerow(row)
        if i!=len(ps[key]):raise ValueError('Family denominator changed')
    (a.output/'summary.json').write_text(json.dumps(dict(status='PASS_COMPLETE_GLOBAL30',tested=total,policy=policy,inputs=source_hashes),indent=2)+'\n')
if __name__=='__main__':main()
