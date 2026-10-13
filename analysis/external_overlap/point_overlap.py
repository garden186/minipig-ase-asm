"""Exact 1-based CpG point overlap with registered 0-based half-open BED peaks.

Outputs source-replicate counts and threshold membership, not enrichment P values.
"""
import argparse,csv,gzip,hashlib,json
from collections import defaultdict
from bisect import bisect_left
from pathlib import Path
def rows(p):
    with (gzip.open(p,'rt',encoding='utf-8') if str(p).endswith('.gz') else Path(p).open(encoding='utf-8')) as h:yield from csv.DictReader(h,delimiter='\t')
def sha(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda:f.read(4*1024**2),b''):h.update(b)
    return h.hexdigest()
def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--points',type=Path,required=True);p.add_argument('--peaks',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    points=list(rows(a.points));by=defaultdict(list)
    for c,pos in sorted({(r['chrom'].removeprefix('chr'),int(r['cpg_pos1'])) for r in points}):
        if pos<1:raise ValueError('Invalid CpG coordinate')
        by[c].append(pos-1)
    hits=defaultdict(set);tracks={};sources=[]
    for spec in rows(a.peaks):
        path=Path(spec['path']);track=spec['track'];minimum=int(spec['minimum_replicates'])
        if minimum<1 or (track in tracks and tracks[track]!=minimum):raise ValueError('Inconsistent replicate threshold')
        tracks[track]=minimum
        if sha(path)!=spec['sha256']:raise ValueError('Peak identity differs')
        sources.append(dict(spec))
        with (gzip.open(path,'rt') if path.suffix=='.gz' else path.open()) as h:
            for line in h:
                if not line.strip() or line.startswith(('#','track','browser')):continue
                v=line.split();chrom=v[0].removeprefix('chr');start,end=int(v[1]),int(v[2])
                if not 0<=start<end:raise ValueError('Invalid BED bounds')
                positions=by.get(chrom,[])
                for pos in positions[bisect_left(positions,start):bisect_left(positions,end)]:hits[chrom,pos+1,track].add(spec['replicate'])
    a.output.mkdir(parents=True,exist_ok=False)
    with (a.output/'point_overlap.tsv').open('w',encoding='utf-8',newline='') as h:
        fields=list(points[0])+[f'{t}_{suffix}' for t in tracks for suffix in ['replicates','overlap']]
        w=csv.DictWriter(h,fields,delimiter='\t',lineterminator='\n');w.writeheader()
        for row in points:
            r=dict(row)
            for track,minimum in tracks.items():
                n=len(hits[r['chrom'].removeprefix('chr'),int(r['cpg_pos1']),track]);r[track+'_replicates']=n;r[track+'_overlap']=int(n>=minimum)
            w.writerow(r)
    (a.output/'inputs.json').write_text(json.dumps(dict(points_sha256=sha(a.points),sources=sources,coordinates='CpG pos1-1 in [BED start,end)'),indent=2)+'\n')
if __name__=='__main__':main()
