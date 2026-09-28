"""Compact genome indexes. No methylation P values are calculated here."""
from __future__ import annotations
import collections, csv, gzip, io, math, mmap, time
from pathlib import Path
import numpy as np
from common import *

AVAIL=np.dtype([('pos','<u4'),('ps','<u4')])
VAR=np.dtype([('pos','<u4'),('ref','u1'),('alt','u1'),('gt','u1'),('gq','<f8'),('ps','<u4'),('dup','u1')])
SPLIT=np.dtype([('pos','<u4'),('h1','u1'),('h2','u1'),('dup','u1')])
ANCHOR=np.dtype(VAR.descr+[('status','u1'),('alt_hap','u1'),('gq_source','u1'),('conversion','u1')])
STATUS=['ELIGIBLE_ANCHOR','DUPLICATE_PHASED_POSITION','NOT_HETEROZYGOUS','NOT_PHASED','MISSING_PS',
        'REFERENCE_MISMATCH','DUPLICATE_SNPSPLIT_POSITION','SNPSPLIT_MISSING','SNPSPLIT_ORDER_MISMATCH',
        'DUPLICATE_GENOTYPE_POSITION','GENOTYPE_DISCORDANT','GQ_MISSING','GQ_BELOW_20']
GT=['.','0|1','1|0','0/1','1/0','0/0','1/1','OTHER','0|0','1|1']
BASE={b:i for i,b in enumerate('ACGT')}
NAN=float('nan')

class HashReader:
    def __init__(self,raw):self.raw=raw;self.hash=hashlib.sha256()
    def read(self,n=-1):
        b=self.raw.read(n);self.hash.update(b);return b
    def seekable(self):return False
    def tell(self):return self.raw.tell()

@contextlib.contextmanager
def checked_text(spec):
    before=snapshot(spec['path'])
    with open(spec['path'],'rb') as raw:
        h=HashReader(raw)
        stream=gzip.GzipFile(fileobj=h,mode='rb') if spec['path'].endswith('.gz') else h
        if stream is h:
            # SNPsplit is plain text; hashing line chunks avoids an in-memory copy.
            text=io.TextIOWrapper(raw,encoding='utf-8',newline='')
            yield text
        else:
            with io.TextIOWrapper(stream,encoding='utf-8',newline='') as text:yield text
            while h.read(4*1024**2):pass
            require(h.hash.hexdigest()==spec['sha256'],'Compressed input hash mismatch: '+spec['path'])
    if not spec['path'].endswith('.gz'): verify(spec['path'],spec['sha256'])
    require(snapshot(spec['path'])==before,'Input changed during reading: '+spec['path'])

def unit_job(job):
    spec=job['spec'];parent=Path(job['parent']);ident=job['identity']
    start=time.monotonic();cache=spec['cache'];verify(cache['gate'],cache['gate_sha256'])
    expected_gate=read(cache['gate']);require(not expected_gate['errors'],'Source cache failed')
    for n,h in cache['sha256'].items():
        require(expected_gate['sha256'].get(n)==h,'Cache manifest disagrees with frozen receipt')
        verify(Path(cache['path'])/n,h)
    out,reuse=begin(parent,ident)
    if reuse:return str(out)
    guard(out)
    cp=read(Path(cache['path'])/'phase_dictionary.json')
    require(len(set(cp))==len(cp),'Duplicate cache PS dictionary values')
    raw=spec['raw'];nexpected=raw['rows']
    ar=np.empty(nexpected,dtype=[('chrom','u1'),*AVAIL.descr]);ps=[];pmap={};n=0
    with checked_text(raw) as f:
        header=f.readline().rstrip('\r\n').split('\t');ix={k:i for i,k in enumerate(header)}
        require(len(ix)==len(header) and all(k in ix for k in ['sample','tissue','chrom','cpg_pos1','ps']), 'Raw header mismatch')
        for line in f:
            v=line.rstrip('\r\n').split('\t');require(len(v)==len(header),'Malformed raw row')
            require(v[ix['sample']]==spec['sample'] and v[ix['tissue']]==spec['tissue'],'Raw unit mismatch')
            c=int(v[ix['chrom']].removeprefix('chr'));pos=int(v[ix['cpg_pos1']]);s=v[ix['ps']]
            require(1<=c<=18 and 0<pos<2**32 and s not in ('','.','NA','None'),'Invalid raw coordinate/PS')
            if s not in pmap:pmap[s]=len(ps);ps.append(s)
            require(n<nexpected,'More raw rows than registered');ar[n]=(c,pos,pmap[s]);n+=1
            if n%500000==0:guard(out);print('RAW='+job['name']+' ROWS='+str(n),flush=True)
    require(n==nexpected,'Raw row count mismatch');write(out/'phase_dictionary.json',ps)
    require(all(x in pmap for x in cp),'Tested PS absent from raw availability')
    mapping=np.array([pmap[x] for x in cp],dtype='<u4');np.save(out/'cache_to_raw_ps.npy',mapping,allow_pickle=False)
    totals=0;hist={};bins=np.array([0,3,5,10,20,50,100,1000,2**33],dtype=np.uint64)
    for c in range(1,19):
        guard(out);part=ar[ar['chrom']==c];a=np.empty(len(part),dtype=AVAIL)
        a['pos']=part['pos'];a['ps']=part['ps'];order=np.lexsort((a['ps'],a['pos']));a=a[order]
        keys=(a['pos'].astype(np.uint64)<<32)|a['ps']
        require(not len(keys) or np.all(keys[1:]>keys[:-1]),'Duplicate raw CpG/PS')
        np.save(out/f'chr{c}.available.npy',a,allow_pickle=False)
        tested=np.load(Path(cache['path'])/f'chr{c}.npy',mmap_mode='r',allow_pickle=False)
        require(all(k in tested.dtype.names for k in ['chrom','pos','ps','a','b','c','d','p','q_all','source_row']),'Unknown cache dtype')
        require(np.all(tested['chrom']==c) and np.all(tested['ps']<len(cp)),'Cache chromosome/PS invalid')
        tk=(tested['pos'].astype(np.uint64)<<32)|mapping[tested['ps']]
        q=np.searchsorted(keys,tk);require(np.all(q<len(keys)) and np.array_equal(keys[q],tk),'Tested row not in raw PS index')
        n1=tested['a'].astype(np.uint64)+tested['b'];n2=tested['c'].astype(np.uint64)+tested['d']
        require(np.all((n1>=3)&(n2>=3)&(n1+n2>=10)),'Tested eligibility changed')
        totals+=len(tested);hist[str(c)]=dict(raw_rows=len(a),tested_rows=len(tested),
            H1_depth_hist=np.histogram(n1,bins)[0].tolist(),H2_depth_hist=np.histogram(n2,bins)[0].tolist())
    require(totals==spec['annotated']['rows'],'Tested row count mismatch')
    return seal(parent,out,ident,dict(sample=spec['sample'],tissue=spec['tissue'],raw_rows=n,
        tested_rows=totals,phase_sets=len(ps),chromosomes=hist,depth_bins=bins.tolist(),seconds=time.monotonic()-start,
        max_rss=psutil.Process().memory_info().rss,sources=spec,
        note='Original count/P/q arrays remain read-only; raw index contains only CpG/PS opportunity'))

def parse_gq(s):
    try:x=float(s)
    except (ValueError,TypeError):return NAN
    return x if math.isfinite(x) and x>=0 else NAN

def gt_code(s):
    return GT.index(s) if s in GT else 7

class Sink:
    """Bounded buffers to append compact rows, then publish sorted NPY per chr."""
    def __init__(self,out,dtype):
        self.out=Path(out);self.dtype=dtype;self.buff={};self.used={};self.count=collections.Counter()
    def add(self,c,record):
        if c not in self.buff:self.buff[c]=np.empty(20000,dtype=self.dtype);self.used[c]=0
        self.buff[c][self.used[c]]=record;self.used[c]+=1;self.count[c]+=1
        if self.used[c]==len(self.buff[c]):self.flush(c)
    def flush(self,c):
        if self.used[c]:
            with (self.out/f'chr{c}.bin').open('ab') as f:self.buff[c][:self.used[c]].tofile(f)
            self.used[c]=0
    def finish(self):
        for c in list(self.buff):self.flush(c)
        for c in range(1,19):
            p=self.out/f'chr{c}.bin';a=np.fromfile(p,dtype=self.dtype) if p.exists() else np.empty(0,dtype=self.dtype)
            if len(a):
                order=np.lexsort((a['alt'],a['ref'],a['pos'])) if 'ref' in a.dtype.names else np.argsort(a['pos'],kind='stable')
                a=a[order]
                # Count duplicate positions including duplicate records of different alleles.
                same=a['pos'][1:]==a['pos'][:-1];a['dup'][:-1]|=same;a['dup'][1:]|=same
            np.save(self.out/f'chr{c}.npy',a,allow_pickle=False)
            if p.exists():p.unlink() # Only a literal internal spool, never a source.

def scan_job(job):
    parent=Path(job['parent']);out,reuse=begin(parent,job['identity'])
    if reuse:return str(out)
    guard(out);spec=job['spec'];role=job['role'];start=time.monotonic();sink=Sink(out,SPLIT if role=='snpsplit' else VAR)
    ps=['.'];pmap={'.':0};counts=collections.Counter();aliases={};last=None;seen_chrom=set();oldchrom=None
    duplicates=collections.defaultdict(list)
    def coordinate(chrom,pos):
        nonlocal last,oldchrom
        cc=chrom.removeprefix('chr')
        if cc not in {str(c) for c in range(1,19)}:return None
        c=int(cc);require(c not in aliases or aliases[c]==chrom,'Ambiguous contig alias');aliases[c]=chrom
        require(0<pos<2**32,'Invalid variant position')
        if c!=oldchrom:
            require(c not in seen_chrom,'Contig repeated in input');seen_chrom.add(c);oldchrom=c;last=None
        require(last is None or pos>=last,'Input not coordinate sorted');last=pos;return c
    with checked_text(spec) as f:
        if role=='snpsplit':
            group=[];loc=None
            def emit_group(g):
                for c,pos,h1,h2 in g:sink.add(c,(pos,h1,h2,int(len(g)>1)))
            for line in f:
                if line.startswith('#') or not line.strip():continue
                counts['records']+=1;v=line.rstrip('\r\n').split('\t');require(len(v)>=5,'Malformed SNPsplit row')
                pos=int(v[2]);c=coordinate(v[1],pos)
                if c is None:continue
                ab=v[-1].upper().split('/');require(len(ab)==2 and all(x in BASE for x in ab),'Non-SNV SNPsplit allele')
                new=(c,pos)
                if new!=loc:emit_group(group);group=[];loc=new
                group.append((c,pos,BASE[ab[0]],BASE[ab[1]]))
                if counts['records']%500000==0:guard(out)
            emit_group(group)
        else:
            col=None;format_cache={};group=[];loc=None
            def emit_group(g):
                duplicate=len(g)>1
                if duplicate:duplicates[loc[0]].append(loc[1])
                for item in g:
                    if item is not None:sink.add(item[0],(*item[1:],int(duplicate)))
            for line in f:
                if line.startswith('##'):continue
                if line.startswith('#CHROM'):
                    h=line.rstrip('\r\n').split('\t');require(h[9:].count(job['sample'])==1,'VCF sample identity missing/ambiguous')
                    col=h.index(job['sample']);continue
                if line.startswith('#'):continue
                require(col is not None,'VCF header missing');v=line.rstrip('\r\n').split('\t');require(len(v)>col,'Truncated VCF record')
                counts['records']+=1;pos=int(v[1]);c=coordinate(v[0],pos)
                if c is None:continue
                new=(c,pos)
                if new!=loc:emit_group(group);group=[];loc=new
                ref,alt=v[3].upper(),v[4].upper()
                if ref not in BASE or alt not in BASE or ref==alt:
                    counts['non_biallelic_SNV']+=1;group.append(None);continue
                names=format_cache.setdefault(v[8],v[8].split(':'));d=dict(zip(names,v[col].split(':')))
                s=d.get('PS','.') if role=='phased' else '.'
                if s in ('','NA'):s='.'
                if s not in pmap:pmap[s]=len(ps);ps.append(s)
                group.append((c,pos,BASE[ref],BASE[alt],gt_code(d.get('GT','.')),parse_gq(d.get('GQ')),pmap[s]))
                if counts['records']%500000==0:guard(out);print('SCAN='+job['name']+' ROWS='+str(counts['records']),flush=True)
            require(col is not None,'VCF header missing');emit_group(group)
    sink.finish();write(out/'phase_dictionary.json',ps)
    for c in range(1,19):
        np.save(out/f'chr{c}.duplicate_positions.npy',np.array(duplicates[c],dtype='<u4'),allow_pickle=False)
    return seal(parent,out,job['identity'],dict(sample=job['sample'],role=role,records=dict(counts),
        chromosome_rows=dict(sink.count),seconds=time.monotonic()-start,source=spec,
        note='Biallelic SNV index only; non-SNVs counted, duplicate-position evidence retained'))

def keys(a):return a['pos'].astype(np.uint64)*16+a['ref'].astype(np.uint64)*4+a['alt']

def lookup(a,k):
    """Index, found; duplicate sources flagged separately before deduplication."""
    ak=keys(a);ix=np.searchsorted(ak,k)
    found=ix<len(a);found[found]&=ak[ix[found]]==k[found]
    return ix,found

def unique_positions(a):
    if not len(a):return a
    return a[np.r_[True,a['pos'][1:]!=a['pos'][:-1]]]

def reference_codes(spec,c,positions):
    lines={}
    with open(spec['index_path'],encoding='utf-8') as f:
        for line in f:
            v=line.rstrip().split('\t');name=v[0].removeprefix('chr')
            require(name not in lines,'Ambiguous FASTA aliases');lines[name]=tuple(map(int,v[1:5]))
    size,offset,width,stride=lines[str(c)]
    require(np.all((positions>0)&(positions<=size)),'Variant outside reference contig')
    pos=positions.astype(np.int64)-1;addresses=offset+(pos//width)*stride+pos%width
    with open(spec['path'],'rb') as f:
        with mmap.mmap(f.fileno(),0,access=mmap.ACCESS_READ) as mm:
            bases=np.frombuffer(mm,dtype='u1');value=bases[addresses].copy();del bases
    lut=np.full(256,255,dtype='u1')
    for b,k in BASE.items():lut[ord(b)]=k;lut[ord(b.lower())]=k
    return lut[value]

def join_arrays(p,g,s,refcodes,genotype_duplicate_positions=None):
    """Vector implementation of the previously reviewed phase.decorate rules."""
    out=np.empty(len(p),dtype=ANCHOR)
    for k in p.dtype.names:out[k]=p[k]
    status=np.zeros(len(p),dtype='u1');pgood=np.isfinite(p['gq']);gq=p['gq'].copy()
    gi,gfound=lookup(g,keys(p));gg=np.zeros(len(p),dtype=VAR)
    if len(g):gg[gfound]=g[gi[gfound]]
    phet=np.isin(p['gt'],[1,2,3,4]);ghet=np.isin(gg['gt'],[1,2,3,4])
    agree=gfound&((phet&ghet)|(np.isin(p['gt'],[5,8])&np.isin(gg['gt'],[5,8]))|(np.isin(p['gt'],[6,9])&np.isin(gg['gt'],[6,9])))
    source=np.where(pgood,1,np.where(agree&np.isfinite(gg['gq']),2,0)).astype('u1')
    fallback=(~pgood)&(source==2);gq[fallback]=gg['gq'][fallback]
    su=unique_positions(s);si=np.searchsorted(su['pos'],p['pos']);sf=si<len(su)
    sf[sf]&=su['pos'][si[sf]]==p['pos'][sf]
    ss=np.zeros(len(p),dtype=SPLIT)
    if len(su):ss[sf]=su[si[sf]]
    h1=np.where(np.isin(p['gt'],[1,3,5]),p['ref'],p['alt'])
    h2=np.where(np.isin(p['gt'],[1,3,6]),p['alt'],p['ref'])
    gdup=gfound&(gg['dup']>0)
    if genotype_duplicate_positions is not None:gdup|=np.isin(p['pos'],genotype_duplicate_positions)
    conditions=[p['dup']>0,~phet,~np.isin(p['gt'],[1,2]),p['ps']==0,p['ref']!=refcodes,
                sf&(ss['dup']>0),~sf,sf&((ss['h1']!=h1)|(ss['h2']!=h2)),
                gdup,gfound&~agree,source==0,(source>0)&(gq<20)]
    for code,mask in enumerate(conditions,1):status[(status==0)&mask]=code
    out['status']=status;out['gq']=gq;out['gq_source']=source
    out['alt_hap']=np.where(p['gt']==2,1,np.where(p['gt']==1,2,0))
    out['conversion']=phet&(((np.minimum(p['ref'],p['alt'])==1)&(np.maximum(p['ref'],p['alt'])==3))|((np.minimum(p['ref'],p['alt'])==0)&(np.maximum(p['ref'],p['alt'])==2)))
    return out

def join_job(job):
    out,reuse=begin(job['parent'],job['identity'])
    if reuse:return str(out)
    guard(out);c=job['chrom'];arrays=[]
    for role in ['phased','genotype','snpsplit']:
        parent=Path(job['paths'][role]);v=read(parent/'validation.json');name=f'chr{c}.npy'
        verify(parent/name,v['sha256'][name]);arrays.append(np.load(parent/name,mmap_mode='r',allow_pickle=False))
    gp=Path(job['paths']['genotype']);name=f'chr{c}.duplicate_positions.npy'
    verify(gp/name,read(gp/'validation.json')['sha256'][name])
    gd=np.load(gp/name,allow_pickle=False)
    result=join_arrays(*arrays,reference_codes(job['reference'],c,arrays[0]['pos']),gd)
    np.save(out/'variants.npy',result,allow_pickle=False)
    states={STATUS[int(k)]:int(v) for k,v in zip(*np.unique(result['status'],return_counts=True))}
    eligible=result[result['status']==0]
    # PS-sorted index supports CpG-PS lookups without repeated full scans.
    order=np.lexsort((eligible['alt'],eligible['ref'],eligible['pos'],eligible['ps']))
    np.save(out/'eligible_by_ps.npy',eligible[order],allow_pickle=False)
    return seal(job['parent'],out,job['identity'],dict(sample=job['sample'],chrom=c,rows=len(result),
         eligible=len(eligible),states=states,phase_dictionary=str(Path(job['paths']['phased'])/'phase_dictionary.json'),
         conversion_ambiguous_eligible=int(eligible['conversion'].sum()),
         note='CT/AG are WGS orientation markers only; no new WGBS fragment assignment'))
