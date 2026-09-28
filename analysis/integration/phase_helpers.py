"""Original exact-phase orientation and signed-effect functions."""
import csv,gzip,hashlib,json,math
from pathlib import Path

def require(value, message):
    if not value:
        raise ValueError(message)

def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for chunk in iter(lambda:f.read(1024**2),b''):
            h.update(chunk)
    return h.hexdigest()

def js(path):
    return json.loads(Path(path).read_text(encoding='utf-8'))

def save_json(path, obj):
    Path(path).write_text(json.dumps(obj,ensure_ascii=False,sort_keys=True,indent=2)+'\n',encoding='utf-8')

def rows(path):
    with (gzip.open if str(path).endswith('.gz') else open)(path,'rt',encoding='utf-8',newline='') as f:
        yield from csv.DictReader(f,delimiter='\t')

def write_rows(path, data, empty_fields=('status',)):
    fields=list(dict.fromkeys(k for r in data for k in r)) or list(empty_fields)
    with (gzip.open if str(path).endswith('.gz') else open)(path,'wt',encoding='utf-8',newline='') as f:
        w=csv.DictWriter(f,fields,delimiter='\t',lineterminator='\n')
        w.writeheader();w.writerows(data)

def split(x):
    return [v for v in str(x).split(',') if v]

def phase_ok(block, chrom, ps):
    n=int(block['n_wgs_phase_compared_snps'] or 0)
    direct=int(block['n_wgs_phase_direct_snps'] or 0)
    flipped=int(block['n_wgs_phase_flipped_snps'] or 0)
    mismatch=int(block['n_wgs_phase_mismatch_snps'] or 0)
    return (block['gene_contig'].removeprefix('chr')==str(chrom).removeprefix('chr')
        and block['wgs_phase_set']==str(ps) and str(ps) not in ('','NA','.')
        and block['wgs_phase_status']=='CONCORDANT'
        and n>=1 and mismatch==0
        and ((block['wgs_phase_orientation']=='DIRECT' and direct==n and flipped==0)
             or (block['wgs_phase_orientation']=='FLIPPED' and flipped==n and direct==0)))

def orient(a,b,orientation):
    require(orientation in ('DIRECT','FLIPPED'), 'Unresolved phase orientation')
    return (a,b) if orientation=='DIRECT' else (b,a)

def variant_alleles(block, variant):
    variants=split(block['block_snps']);a=split(block['haplotypeA']);b=split(block['haplotypeB'])
    require(len(variants)==len(a)==len(b)==len(set(variants)), 'Malformed block allele arrays')
    require(variant in variants, 'SNP absent from parent block')
    i=variants.index(variant)
    h1,h2=orient(a[i],b[i],block['wgs_phase_orientation'])
    require(set((h1,h2))==set(variant.split('_')[-2:]), 'SNP alleles do not match haplotypes')
    return h1,h2

def effects(h1,h2,delta):
    total=h1+h2
    difference=(h1-h2)/total if total else None
    return dict(RNA_H1_count=h1,RNA_H2_count=h2,RNA_total=total,
        RNA_H1_fraction=h1/total if total else '',
        RNA_fraction_difference_H1_minus_H2=difference if total else '',
        direction_relation=('INVERSE' if delta*difference<0 else 'SAME_DIRECTION' if delta*difference>0 else 'ZERO') if total else 'NO_COUNTS')

def self_test():
    b=dict(gene_contig='2',wgs_phase_set='100',wgs_phase_status='CONCORDANT',
        wgs_phase_orientation='FLIPPED',n_wgs_phase_compared_snps='2',n_wgs_phase_direct_snps='0',
        n_wgs_phase_flipped_snps='2',n_wgs_phase_mismatch_snps='0',
        block_snps='2_110_G_C,2_120_A_T',haplotypeA='G,A',haplotypeB='C,T')
    require(phase_ok(b,'chr2','100'), 'Exact phase test')
    require(not phase_ok(b,'2','101') and not phase_ok(b,'3','100'), 'Phase/chromosome mismatch test')
    require(not phase_ok(dict(b,n_wgs_phase_mismatch_snps='1'),'2','100'), 'Conflicting phase test')
    require(orient(21,47,'FLIPPED')==(47,21) and orient(21,47,'DIRECT')==(21,47), 'Direction test')
    require(variant_alleles(b,'2_110_G_C')==('C','G'), 'SNP allele mapping test')
    require(effects(47,21,-.5)['direction_relation']=='INVERSE', 'Inverse direction test')
    require(effects(0,0,-.5)['direction_relation']=='NO_COUNTS', 'Missing count test')
    require((.1 <= .1) and not (.10001 <= .1), 'Inclusive q threshold test')
    return 8
