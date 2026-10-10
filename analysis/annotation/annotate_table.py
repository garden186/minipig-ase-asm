"""Append original Ensembl115 annotations to a CpG TSV with chrom/cpg_pos1."""
import argparse,csv,gzip
from pathlib import Path
from common import EXTRA
from reference import Index
def open_text(p,mode):return gzip.open(p,mode+'t',encoding='utf-8',newline='') if str(p).endswith('.gz') else Path(p).open(mode,encoding='utf-8',newline='')
def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ['input','output','index']:p.add_argument('--'+name,type=Path,required=True)
    a=p.parse_args();idx=Index(a.index,verify=True)
    try:
        with open_text(a.input,'r') as h,open_text(a.output,'x') as out:
            reader=csv.DictReader(h,delimiter='\t');fields=list(reader.fieldnames)
            if not {'chrom','cpg_pos1'}<=set(fields):raise ValueError('Missing coordinates')
            writer=csv.DictWriter(out,fields+[k for k in EXTRA if k not in fields],delimiter='\t',lineterminator='\n');writer.writeheader()
            for r in reader:
                ann=idx.annotation(r['chrom'].removeprefix('chr'),int(r['cpg_pos1']))
                for k,v in ann.items():
                    if k in r and str(r[k])!=str(v):raise ValueError('Existing annotation differs: '+k)
                writer.writerow(dict(r,**ann))
    finally:idx.close()
if __name__=='__main__':main()
