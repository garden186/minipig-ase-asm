"""Original exact-PS/GQ20 count eligibility checks."""
COUNTS=['H1_methylated', 'H1_unmethylated', 'H2_methylated', 'H2_unmethylated']
SOURCE_FIELDS=['sample', 'tissue', 'chrom', 'cpg_pos1', 'ps', 'genotype_disrupted', 'H1_methylated', 'H1_unmethylated', 'H2_methylated', 'H2_unmethylated', 'H1_depth', 'H2_depth', 'baseline_H1_depth', 'baseline_H2_depth', 'assignment_GQ_cutoff', 'observation_status']

def parse(r,t,p):
    if set(r)!=set(SOURCE_FIELDS): raise ValueError('Unexpected GQ20 row schema')
    if tuple(r[k] for k in ('sample','tissue','chrom'))!=tuple(str(t[k]) for k in ('sample','tissue','chrom')):
        raise ValueError('GQ20 row sample/tissue/chromosome differs')
    pos=int(r['cpg_pos1']);ps=r['ps']
    if not t['start0']<=pos-1<t['end0'] or ps in ('','.'): raise ValueError('Invalid core coordinate or exact phase set')
    x=tuple(int(r[k]) for k in COUNTS);n1=x[0]+x[1];n2=x[2]+x[3]
    if min(x)<0 or n1!=int(r['H1_depth']) or n2!=int(r['H2_depth']): raise ValueError('Invalid haplotype depth arithmetic')
    if min(int(r['baseline_H1_depth']),int(r['baseline_H2_depth']))<0: raise ValueError('Negative baseline depth')
    if r['genotype_disrupted'] not in ('0','1') or r['assignment_GQ_cutoff']!='20': raise ValueError('Invalid GQ/mask policy')
    depth_status='BOTH_HAPLOTYPES_OBSERVED' if min(n1,n2) else 'ONE_HAPLOTYPE_UNOBSERVED' if n1+n2 else 'NO_VALID_HAPLOTYPE_AFTER_GQ20'
    if r['observation_status']!=depth_status: raise ValueError('GQ20 assignment-depth status differs')
    status=('CPG_GENOTYPE_DISRUPTED' if r['genotype_disrupted']=='1' else
            'INSUFFICIENT_HAPLOTYPE_DEPTH' if min(n1,n2)<p['minimum_each_haplotype'] else
            'INSUFFICIENT_COMBINED_DEPTH' if n1+n2<p['minimum_combined_depth'] else 'TESTED')
    return (pos,ps),x,n1,n2,status
