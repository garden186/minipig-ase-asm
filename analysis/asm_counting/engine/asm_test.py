"""Explicit, bounded pilot Fisher tests with one combined 30-unit BH family."""
from __future__ import annotations
import argparse
import math
from pathlib import Path
from common import ANIMALS, TISSUES, read_json, rows, sha, write_table, write_json, stamp, verify_package
from reassign import ASM_FIELDS


def fisher_exact(a,b,c,d):
    """Probability-ordering two-sided test, with stable integer hypergeometric weights."""
    if any(type(v) is not int or v<0 for v in (a,b,c,d)):
        raise ValueError('Fisher counts must be nonnegative integers')
    n1,n2,col = a+b,c+d,a+c
    if not n1 or not n2: raise ValueError('Both haplotypes need observations')
    total=n1+n2
    # The bounded diagnostic cap prevents accidental expensive genome-scale use.
    if total>10000: raise ValueError('Pilot Fisher depth exceeds the documented 10000-fragment bound')
    lo,hi=max(0,n1-(total-col)),min(n1,col)
    observed=math.comb(col,a)*math.comb(total-col,n1-a)
    numerator=sum(w for x in range(lo,hi+1)
                  if (w:=math.comb(col,x)*math.comb(total-col,n1-x)) <= observed)
    value=numerator/math.comb(total,n1)
    # Never output a floating-point zero as an exact P value.
    return max(math.nextafter(0.0,1.0),min(1.0,value))


def bh(values):
    if any(not math.isfinite(x) or not 0<=x<=1 for x in values):
        raise ValueError('Invalid P value')
    result=[0.0]*len(values);bound=1.0
    ordered=sorted(range(len(values)),key=lambda i:values[i])
    for rank in range(len(ordered),0,-1):
        i=ordered[rank-1];bound=min(bound,values[i]*len(values)/rank);result[i]=bound
    return result


def test_rows(source_rows, minimum):
    if type(minimum) is not int or minimum<1: raise ValueError('Specify a positive minimum haplotype depth')
    result=[];indices=[];ps=[];seen=set()
    for r in source_rows:
        key=tuple(r[k] for k in ('sample','tissue','chrom','cpg_pos1','ps'))
        if key in seen: raise ValueError('Duplicate sample-tissue-CpG-PS hypothesis')
        seen.add(key)
        h=[int(r[k]) for k in ('H1_methylated','H1_unmethylated','H2_methylated','H2_unmethylated')]
        if any(x<0 for x in h) or [sum(h[:2]),sum(h[2:])] != [int(r['H1_depth']),int(r['H2_depth'])]:
            raise ValueError('Invalid haplotype count arithmetic')
        if r['genotype_disrupted'] not in ('0','1'): raise ValueError('Invalid genotype mask')
        reason='CPG_GENOTYPE_DISRUPTED' if r['genotype_disrupted']=='1' else 'INSUFFICIENT_HAPLOTYPE_DEPTH' if min(sum(h[:2]),sum(h[2:]))<minimum else 'TESTED'
        p=fisher_exact(*h) if reason=='TESTED' else '.'
        result.append([r[k] for k in ASM_FIELDS]+[reason,p,'.','EXPLORATORY_PILOT_ONLY'])
        if reason=='TESTED':indices.append(len(result)-1);ps.append(p)
    for i,q in zip(indices,bh(ps)): result[i][-2]=q
    return result,len(ps)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--runs',type=Path,nargs='+',required=True,help='Completed extraction or same-cutoff reassignment runs covering all 30 units')
    p.add_argument('--min-haplotype-depth',type=int,required=True)
    p.add_argument('--allow-pilot',action='store_true',help='Acknowledge that these P/q values are bounded diagnostics, not genome-wide results')
    p.add_argument('--output-root',type=Path,required=True)
    args=p.parse_args();package_sha=verify_package()
    if not args.allow_pilot: raise ValueError('This release implements bounded diagnostic tests only; specify --allow-pilot')
    all_rows=[];units=set();cutoffs=set();sources=[]
    for root in args.runs:
        gate=read_json(root/'validation.json')
        if gate['errors']: raise ValueError('Input run failed validation')
        if gate['status']=='PASS_COHORT_CPG_ASSIGNMENT_PILOT':
            request=read_json(root/'run_request.json');cohort=request['units'];cutoff=None
            if request['package_sha256']!=package_sha:raise ValueError('Extraction release differs')
        elif gate['status']=='PASS_COHORT_GQ_REASSIGNMENT' and gate['source_mode']=='pilot':
            cohort=gate['units'];cutoff=gate['min_gq']
            if gate['package_sha256']!=package_sha:raise ValueError('Reassignment release differs')
        else: raise ValueError('Unexpected analysis scope or completion status')
        cutoffs.add(cutoff)
        for unit in cohort:
            key=unit['sample'],unit['tissue']
            if key in units: raise ValueError('Overlapping input cohorts')
            units.add(key)
        for item in gate['regions']:
            directory=root/item['directory']
            if sha(directory/'validation.json')!=item['validation_sha256']: raise ValueError('Region registry changed')
            rg=read_json(directory/'validation.json')
            if rg['errors']: raise ValueError('Region failed')
            name='cpg_asm_inputs.tsv.gz'
            if sha(directory/name)!=rg['sha256'][name]: raise ValueError('ASM input table changed')
            all_rows.extend(rows(directory/name))
            if len(all_rows)>200000: raise ValueError('Pilot table budget exceeded')
        sources.append(dict(directory=str(root.resolve()),validation_sha256=sha(root/'validation.json')))
    if units!={(a,t) for a in ANIMALS for t in TISSUES} or len(cutoffs)!=1:
        raise ValueError('Require all 30 units exactly once with one common GQ policy')
    all_rows.sort(key=lambda r:(r['sample'],r['tissue'],int(r['chrom']),int(r['cpg_pos1']),r['ps']))
    result,ntested=test_rows(all_rows,args.min_haplotype_depth)
    output=args.output_root/('cpg_ASM_pilot_diagnostic_v0.3.0_'+stamp());output.mkdir(parents=True,exist_ok=False)
    write_table(output/'cpg_ASM_pilot_diagnostic.tsv.gz',ASM_FIELDS+['test_status','fisher_two_sided_p','combined_30_unit_BH_q','scope'],result)
    write_json(output/'validation.json',dict(status='PASS_BOUNDED_FISHER_ARITHMETIC',errors=[],sources=sources,
        package_sha256=package_sha,minimum_per_haplotype=args.min_haplotype_depth,min_gq=next(iter(cutoffs)),
        hypotheses_tested=ntested,input_rows=len(result),multiple_testing_family='All eligible CpG-PS tests from all 30 units in these three pilot intervals; no effect prefilter.',
        significance_threshold_adopted=False,genome_wide_inference=False,
        absent_phase_rows='No assigned fragment methylation observation; not a negative ASM result.',
        genotype_mask_limit='The inherited disruption flag is not proof of confidently diploid CpG preservation.',
        sha256={'cpg_ASM_pilot_diagnostic.tsv.gz':sha(output/'cpg_ASM_pilot_diagnostic.tsv.gz')}))
    print('COMPLETE='+str(output))


if __name__=='__main__':main()
