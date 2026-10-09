"""Recount validated fragment evidence using the original GQ >= 20 policy."""
import argparse,json,sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parent/'gq20'))
from recount import count
from bootstrap import sha
def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--tiles',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    paths=json.loads(a.tiles.read_text());a.output.mkdir(parents=True,exist_ok=False);records=[]
    for i,source in enumerate(paths):
        source=Path(source);out=a.output/f'tile_{i:06d}';gate=count(source,out,sha(source/'validation.json'))
        records.append(dict(**gate['task'],path=str((out/'cpg_depth_by_phase_set.tsv.gz').resolve())))
    (a.output/'counts_manifest.json').write_text(json.dumps(records,indent=2)+'\n')
if __name__=='__main__':main()
