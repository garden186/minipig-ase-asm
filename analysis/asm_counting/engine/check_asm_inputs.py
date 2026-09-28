"""Independently check the single-CpG contingency-table export."""
import math
from common import rows
from reassign import ASM_FIELDS


def check(directory,min_gq=None):
    masks={tuple(r[k] for k in ('sample','tissue','chrom','cpg_pos1')):r['genotype_disrupted']
           for r in rows(directory/'cpg_assignment.tsv.gz')}
    expected={}
    for r in rows(directory/'cpg_haplotype_by_ps.tsv.gz'):
        key=tuple(r[k] for k in ('sample','tissue','chrom','cpg_pos1','ps'))
        h=tuple(int(r[f'FRAGMENT_{hap}_{call}']) for hap in ('H1','H2') for call in ('METHYLATED','UNMETHYLATED'))
        if sum(h):expected[key]=h
    seen=set()
    for r in rows(directory/'cpg_asm_inputs.tsv.gz'):
        if list(r)!=ASM_FIELDS:raise ValueError('ASM input export schema differs')
        key=tuple(r[k] for k in ('sample','tissue','chrom','cpg_pos1','ps'))
        if key in seen or key not in expected:raise ValueError('Duplicate or unknown ASM input key')
        seen.add(key)
        h=tuple(int(r[k]) for k in ('H1_methylated','H1_unmethylated','H2_methylated','H2_unmethylated'))
        if h!=expected[key] or r['genotype_disrupted']!=masks[key[:4]]:raise ValueError('ASM counts or genotype mask differ')
        depths=[sum(h[:2]),sum(h[2:])]
        if depths!=[int(r['H1_depth']),int(r['H2_depth'])]:raise ValueError('ASM depths differ')
        fractions=[h[0]/depths[0] if depths[0] else None,h[2]/depths[1] if depths[1] else None]
        for field,value in zip(('H1_methylation_fraction','H2_methylation_fraction','delta_m'),
            fractions+[fractions[0]-fractions[1] if all(depths) else None]):
            if value is None:
                if r[field]!='.':raise ValueError('Missing haplotype measurement became a number')
            elif not math.isclose(float(r[field]),value,abs_tol=1e-12,rel_tol=1e-12):raise ValueError('ASM fraction arithmetic differs')
        status='BOTH_HAPLOTYPES_OBSERVED' if all(depths) else 'ONE_HAPLOTYPE_UNOBSERVED'
        if r['observation_status']!=status:raise ValueError('ASM opportunity status differs')
        if min_gq is None:
            if r['assignment_GQ_cutoff']!='NONE':raise ValueError('Unexpected GQ cutoff')
        elif float(r['assignment_GQ_cutoff'])!=min_gq:raise ValueError('GQ cutoff differs')
    if seen!=set(expected):raise ValueError('ASM export omitted observed fragment-CpG-PS rows')
    return len(seen)
