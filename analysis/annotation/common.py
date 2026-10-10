"""Immutable annotation contracts and bounded, source-preserving IO."""
import csv
import gzip
import hashlib
import json
import os
from pathlib import Path
import shutil
import sys
import time

VERSION = '0.1.0'
ROOT = Path(__file__).resolve().parent
GTF_SHA256 = '1d94edba462adb702ee9bf16a68e455684d48a5b07fc3b12871661bc40079cbe'
IDS = ('0235','0242','0276','0309','0326','0349','0355','0377','0384','0385')
TISSUES = ('Heart','Kidney','Liver')
POLICY = dict(assembly='Sscrofa11.1',ensembl_release=115,
    gene_scope='protein_coding',transcript_scope='All transcripts of protein-coding genes',
    promoter_tss_inclusive=[-2000,500],canonical_tag='Ensembl_canonical',
    coordinate='Reference plus-strand C anchor, 1-based; internal half-open intervals',
    primary_order=['PROMOTER','FIVE_PRIME_UTR','CDS','THREE_PRIME_UTR','INTRON','OTHER_GENE_BODY'],
    alternatives='Retained separately; no longest-transcript fallback',
    inference='No new tests, filters, BH adjustment, ASE linkage or causality inference',
    scope='18 autosomes; preserve all input CpG-phase-set rows including untested rows')
EXTRA = ['annotation_profile_id','annotation_status','annotation_gene_id','annotation_gene_name',
    'annotation_primary_feature','annotation_gene_count','annotation_transcript_count',
    'annotation_canonical_promoter','annotation_any_transcript_promoter',
    'annotation_dyad_crosses_segment_boundary']
REQUIRED = ['sample','tissue','chrom','cpg_pos1','ps','genotype_disrupted',
    'H1_methylated','H1_unmethylated','H2_methylated','H2_unmethylated','H1_depth','H2_depth',
    'test_status','fisher_two_sided_p','BH_q_within_animal','ASM_status']

def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
    return h.hexdigest()

def digest(obj):return hashlib.sha256(json.dumps(obj,sort_keys=True,separators=(',',':')).encode()).hexdigest()
def read_json(p):return json.loads(Path(p).read_text(encoding='utf-8'))
def write_json(p,obj):
    p=Path(p); tmp=p.with_name(p.name+'.writing')
    tmp.write_text(json.dumps(obj,indent=2,sort_keys=True)+'\n',encoding='utf-8');os.replace(tmp,p)

def open_text(p,mode='rt'):
    return gzip.open(p,mode,encoding='utf-8',newline='') if str(p).endswith('.gz') else open(p,mode.replace('t',''),encoding='utf-8',newline='')

def rows(p):
    with open_text(p) as f:yield from csv.DictReader(f,delimiter='\t')

def writer(p,fields):
    f=gzip.open(p,'xt',encoding='utf-8',newline='',compresslevel=1)
    w=csv.DictWriter(f,fieldnames=fields,delimiter='\t',lineterminator='\n',extrasaction='raise');w.writeheader()
    return f,w

def contained(root,rel):
    root=Path(root).resolve();p=(root/rel).resolve()
    if not p.is_relative_to(root):raise ValueError('Path escaped its registered root')
    return p

def verify_release():
    m=read_json(ROOT/'release_manifest.json')
    for name,h in m['sha256'].items():
        if sha(contained(ROOT,name))!=h:raise ValueError('Release checksum differs: '+name)
    return sha(ROOT/'release_manifest.json')

def check_gate(p):
    g=read_json(p)
    if not str(g.get('status','')).startswith('PASS') or g.get('errors')!=[]:
        raise ValueError('Source requires a completed PASS gate with zero errors: '+str(p))
    return g

class Guard:
    def __init__(self,output,max_gib=8,min_free_gib=20,max_hours=48):
        self.output=Path(output);self.max_gib=max_gib;self.min_free_gib=min_free_gib
        self.max_hours=max_hours;self.start=time.monotonic();self.peak=None
    def __call__(self):
        if sys.platform!='win32':
            import resource
            rss=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*(1 if sys.platform=='darwin' else 1024)
            self.peak=max(self.peak or 0,rss)
            if rss>self.max_gib*1024**3:raise MemoryError('Annotation process memory limit reached')
        if shutil.disk_usage(self.output).free<self.min_free_gib*1024**3:raise OSError('Free disk reserve reached')
        if time.monotonic()-self.start>self.max_hours*3600:raise TimeoutError('Annotation runtime limit reached')

def available_gib():
    p=Path('/proc/meminfo')
    if not p.exists():return None
    for line in p.read_text().splitlines():
        if line.startswith('MemAvailable:'):return int(line.split()[1])/1024**2
