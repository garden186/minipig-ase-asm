"""One sequential decompression per source, source-preserving chromosome cache."""
import gzip
import hashlib
import io
import math
from collections import Counter
from pathlib import Path
import numpy as np
from common import *

DTYPE=np.dtype([('chrom','u1'),('pos','<u4'),('ps','<u4'),('a','<u4'),('b','<u4'),
                ('c','<u4'),('d','<u4'),('p','<f8'),('q_unit','<f8'),('q_tissue','<f8'),
                ('q_all','<f8'),('source_row','<u4')])

class DigestReader:
    def __init__(self, raw): self.raw=raw; self.hash=hashlib.sha256()
    def read(self,n=-1):
        b=self.raw.read(n); self.hash.update(b); return b
    def seekable(self): return False
    def tell(self): return self.raw.tell()

def measure(fields,index,spec):
    def get(k): return fields[index[k]]
    require(get('sample')==spec['sample'] and get('tissue')==spec['tissue'],'Sample identity mismatch')
    require(get('genotype_disrupted')=='0' and get('test_status')=='TESTED','Masked/untested input row')
    chrom=int(get('chrom').removeprefix('chr')); pos=int(get('cpg_pos1')); ps=get('ps')
    require(1<=chrom<=18 and 0<pos<2**32 and ps not in ('','.','NA','None'),'Invalid coordinate/PS')
    counts=tuple(int(get(k)) for k in COUNTS)
    require(all(0<=v<2**32 for v in counts),'Count overflow/negative count')
    a,b,c,d=counts; n1=a+b; n2=c+d
    require(min(n1,n2)>=3 and n1+n2>=10,'Eligibility differs from GQ20 each3+total10')
    require((int(get('H1_depth')),int(get('H2_depth')))==(n1,n2),'Depth differs from counts')
    for k,want in [('H1_methylation',a/n1),('H2_methylation',c/n2),('delta_M',a/n1-c/n2)]:
        require(math.isclose(float(get(k)),want,abs_tol=1e-12),'Effect/count mismatch: '+k)
    p=float(get('fisher_two_sided_p')); qs=[float(get(k)) for k in QCOLS]
    require(math.isfinite(p) and 0<p<=1,'Invalid or unrecorded underflow P')
    require(all(math.isfinite(q) and p<=q+1e-14<=1+1e-14 for q in qs),'Invalid stored q')
    return chrom,pos,ps,counts,p,qs

def process_unit(task):
    spec,parent,identity,guards=task
    reuse=completed(parent,identity)
    if reuse:
        verify_file(spec['path'],spec['sha256'])
        return str(reuse)
    out=attempt(parent); guard=Guard(out,**guards); guard()
    source=Path(spec['path']); before=snapshot(source)
    expected=int(spec['rows']); require(0<expected<2**32,'Invalid source row count')
    ar=np.empty(expected,dtype=DTYPE); pss=[]; psmap={}; texts={}; histogram=Counter(); header=None
    aliases={}; count=0
    try:
        with source.open('rb') as raw:
            hasher=DigestReader(raw)
            with io.TextIOWrapper(gzip.GzipFile(fileobj=hasher,mode='rb'),encoding='utf-8',newline='') as src:
                header=src.readline().rstrip('\r\n').split('\t')
                require(len(set(header))==len(header),'Duplicate source columns')
                index={k:i for i,k in enumerate(header)}
                require(all(k in index for k in ['sample','tissue','chrom','cpg_pos1','ps','genotype_disrupted',
                    'test_status','H1_depth','H2_depth','H1_methylation','H2_methylation','delta_M',
                    'fisher_two_sided_p',*COUNTS,*QCOLS]),'Missing source columns')
                require(any(k.startswith('annotation_') for k in header),'Source annotation absent')
                for line in src:
                    require(count<expected,'More source rows than registered')
                    fields=line.rstrip('\r\n').split('\t')
                    require(len(fields)==len(header),'Malformed source row')
                    chrom,pos,ps,counts,p,qs=measure(fields,index,spec)
                    alias=fields[index['chrom']]
                    require(chrom not in aliases or aliases[chrom]==alias,'Ambiguous chromosome aliases')
                    aliases[chrom]=alias
                    if ps not in psmap: psmap[ps]=len(pss); pss.append(ps)
                    ar[count]=(chrom,pos,psmap[ps],*counts,p,*qs,count+1)
                    if chrom not in texts:
                        texts[chrom]=gzip.open(out/f'chr{chrom}.source.tsv.gz','wt',encoding='utf-8',newline='',compresslevel=1)
                        texts[chrom].write('source_row\t'+'\t'.join(header)+'\n')
                    texts[chrom].write(str(count+1)+'\t'+line.rstrip('\r\n')+'\n')
                    histogram[chrom]+=1; count+=1
                    if count%100000==0:
                        guard(); print('CACHE='+spec['sample']+'-'+spec['tissue']+' ROWS='+str(count),flush=True)
            while hasher.read(1024**2): pass
            require(hasher.hash.hexdigest()==spec['sha256'],'Full compressed source hash mismatch')
        require(snapshot(source)==before,'Source changed during scan')
        require(count==expected,'Source row count differs')
    finally:
        for f in texts.values(): f.close()
    guard()
    order=np.lexsort((ar['ps'],ar['pos'],ar['chrom']))
    ar=ar[order]; del order
    if len(ar)>1:
        duplicate=(ar['chrom'][1:]==ar['chrom'][:-1]) & (ar['pos'][1:]==ar['pos'][:-1]) & (ar['ps'][1:]==ar['ps'][:-1])
        require(not bool(duplicate.any()),'Duplicate same-CpG exact-PS source rows')
    files=[]
    for chrom,n in sorted(histogram.items()):
        part=ar[ar['chrom']==chrom]
        np.save(out/f'chr{chrom}.npy',part,allow_pickle=False)
        check=np.load(out/f'chr{chrom}.npy',mmap_mode='r',allow_pickle=False)
        require(np.array_equal(part,check),'Serialized cache differs')
        files.extend([f'chr{chrom}.npy',f'chr{chrom}.source.tsv.gz'])
    write(out/'phase_dictionary.json',pss); files.append('phase_dictionary.json')
    return seal(parent,out,identity,files,dict(sample=spec['sample'],tissue=spec['tissue'],rows=count,
        chromosome_rows=dict(histogram),source=spec,source_snapshot=before,header=header,
        counts_P_q_preserved=True,annotation='all original columns retained in chromosome source sidecars',
        columns=DTYPE.descr,policy=POLICY))
