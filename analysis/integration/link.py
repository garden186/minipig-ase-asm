"""Join prepared promoter CpG profiles to ASE blocks using original phase/transcript rules.

Input schemas are described in docs/INTEGRATION.md; no new significance test.
"""
import argparse,json,sys
from pathlib import Path
import phase_helpers as L
from link_engine import join
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'annotation'))
from reference import Index
def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ['cpgs','blocks','sites','genes','annotation-index','output']:p.add_argument('--'+name,type=Path,required=True)
    a=p.parse_args();a.output.mkdir(parents=True,exist_ok=False);L.self_test()
    cpgs=list(L.rows(a.cpgs));blocks=list(L.rows(a.blocks));sites=list(L.rows(a.sites))
    genes={(r['sample'],r['tissue'],r['gene_id']):r for r in L.rows(a.genes)}
    index=Index(a.annotation_index,verify=True)
    try:bp,sp,audit,attempts=join(cpgs,blocks,sites,index,genes)
    finally:index.close()
    # This new interface records actual inputs, rather than claiming a historical source receipt.
    for group in [bp,sp]:
        for row in group:row['phase_evidence']='INPUT_V026_BLOCK; SEE_INPUT_HASHES'
    for name,data in [('block_links',bp),('site_links',sp),('cpg_audit',audit),('block_attempts',attempts)]:L.write_rows(a.output/(name+'.tsv'),data)
    L.save_json(a.output/'input_hashes.json',{name:L.sha(getattr(a,name)) for name in ['cpgs','blocks','sites','genes']})
if __name__=='__main__':main()
