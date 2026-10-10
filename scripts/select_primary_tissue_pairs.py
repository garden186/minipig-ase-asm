"""Select final pair-level ASE enrichment without the retired omnibus gate."""
from pathlib import Path
import argparse,csv,gzip,json

def selected(row):
    q=row.get('tissue_enrichment_q_global_BH')
    return q not in (None,'','NA') and float(q)<=.05 and str(row.get('is_primary_recurrent')).lower() in ('1','true','yes')

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--contrasts',type=Path,required=True)
    p.add_argument('--metadata',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    a=p.parse_args()
    m=json.loads(a.metadata.read_text(encoding='utf-8'))
    if m['parameters']['min_informative_animals']!=5:raise ValueError('Primary analysis requires minimum five informative animals')
    if m['parameters']['fdr']!=.05:raise ValueError('Primary FDR must equal 0.05')
    with (gzip.open if a.contrasts.suffix=='.gz' else open)(a.contrasts,'rt',encoding='utf-8-sig',newline='') as f:
        reader=csv.DictReader(f,delimiter='\t');fields=reader.fieldnames;rows=list(reader)
    if a.contrasts.name not in ('threshold_15_tissue_contrasts_all.tsv.gz','threshold_15_tissue_contrasts_all.tsv'):
        raise ValueError('Supply the original named threshold-15 contrast output')
    summary=m['threshold_summaries']['15']
    if len(rows)!=summary['n_cross_tissue_candidate_contrasts']:
        raise ValueError('Incomplete contrast table or mismatched metadata')
    if sum(r['is_contrast_evaluable']=='1' for r in rows)!=summary['n_tissue_contrasts_tested']:
        raise ValueError('Evaluable contrast denominator differs from metadata')
    hits=[r for r in rows if selected(r)]
    keys=[(r['gene_id'],r['target_tissue']) for r in hits]
    if len(keys)!=len(set(keys)):raise ValueError('Duplicate gene-tissue pair')
    with a.output.open('x',encoding='utf-8',newline='') as f:
        w=csv.DictWriter(f,fields,delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(hits)
    print(json.dumps(dict(selected_pairs=len(hits),complete_contrasts=len(rows),omnibus_gate=False)))

if __name__=='__main__':main()
