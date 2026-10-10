"""Per-animal physical-CpG tissue unions under nominal P and global30 BH.

Final composition policy: six displayed features, with canonical exon overlap
reported separately; all categories remain in the denominator.
"""
import argparse,csv,gzip,json
from collections import defaultdict,Counter
from pathlib import Path
from reference import Index
FEATURES=['Promoter','5′ UTR','CDS','3′ UTR','Exon','Intron']
def category(payload):
    summary=payload['summary'];canon=set()
    for g in payload['genes']:
        if g['gene_id']==summary['gene_id']:canon.update(g['canonical_features'].split(';'))
    primary=summary['primary_feature']
    label={'PROMOTER':'Promoter','FIVE_PRIME_UTR':'5′ UTR','CDS':'CDS','THREE_PRIME_UTR':'3′ UTR','INTRON':'Intron'}.get(primary,'No direct annotation' if primary=='NO_DIRECT_PROTEIN_CODING_RELATION' else 'Other gene')
    return label,int('EXON' in canon)
def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--inputs',nargs='+',type=Path,required=True);p.add_argument('--index',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    unions=defaultdict(set);units=set()
    for path in a.inputs:
        with (gzip.open(path,'rt',encoding='utf-8') if path.suffix=='.gz' else path.open(encoding='utf-8')) as h:
            for r in csv.DictReader(h,delimiter='\t'):
                units.add((r['sample'],r['tissue']))
                if r['test_status']!='TESTED':continue
                key=(r['chrom'].removeprefix('chr'),int(r['cpg_pos1']))
                if float(r['fisher_two_sided_p'])<.05:unions['Fisher P',r['sample']].add(key)
                if float(r['BH_q_all30'])<=.05:unions['Global BH',r['sample']].add(key)
    animals={s for s,t in units}
    if len(animals)!=10 or units!={(s,t) for s in animals for t in ['Heart','Kidney','Liver']}:raise ValueError('Final composition requires all 30 units')
    a.output.mkdir(parents=True,exist_ok=False);idx=Index(a.index,verify=True);cache={};result=[];totals=[]
    try:
        for criterion in ['Fisher P','Global BH']:
            for animal in sorted(animals):
                counts=Counter();exons=0;sites=unions[criterion,animal]
                for chrom,pos in sites:
                    key=chrom,pos
                    if key not in cache:cache[key]=category(idx.payload(idx.locate(chrom,pos-1)))
                    label,exon=cache[key];counts[label]+=1;exons+=exon
                for f in FEATURES:
                    n=exons if f=='Exon' else counts[f]
                    result.append(dict(criterion=criterion,sample=animal,feature=f,count=n,total_CpGs=len(sites),percent=100*n/len(sites) if sites else ''))
                totals.append(dict(criterion=criterion,sample=animal,total_CpGs=len(sites),exclusive_categories=dict(counts)))
    finally:idx.close()
    with (a.output/'composition.tsv').open('w',encoding='utf-8',newline='') as h:
        w=csv.DictWriter(h,list(result[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(result)
    (a.output/'denominators.json').write_text(json.dumps(totals,indent=2)+'\n')
if __name__=='__main__':main()
