#!/usr/bin/env python3
"""Portable sequential tiled execution of the original fragment counting engine.

Linux; no change to read filters, conversion handling, mate grouping or exact PS.
Use --chromosomes/--stop-after-tiles for an explicit small first run.
"""
import argparse,csv,json,sys
from pathlib import Path
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE/'engine'))
from common import ROLES,sha,signature,write_json,read_json
from run import preflight,complete_region
from genome_support import tiles,fresh_indexes,prepare_tile_splits,HALO
from worker import run as count_region

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--manifest',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--chromosomes',default=','.join(map(str,range(1,19))))
    p.add_argument('--tile-bp',type=int,default=1000000)
    p.add_argument('--stop-after-tiles',type=int)
    p.add_argument('--samtools',default='samtools');p.add_argument('--bcftools',default='bcftools')
    a=p.parse_args();a.output=a.output.resolve();a.output.mkdir(parents=True,exist_ok=True)
    with a.manifest.open(encoding='utf-8-sig',newline='') as h:data=list(csv.DictReader(h,delimiter='\t'))
    if not data:raise ValueError('Empty manifest')
    keys=[(r['sample'],r['tissue']) for r in data]
    if len(keys)!=len(set(keys)):raise ValueError('Duplicate sample/tissue')
    for r in data:
        for role in ROLES:
            z=Path(r[role]).expanduser()
            r[role]=str((a.manifest.parent/z if not z.is_absolute() else z).resolve())
        r['expected_bam_bytes']=str(Path(r['bam']).stat().st_size)
    code={str(f.relative_to(HERE)):sha(f) for f in sorted(HERE.rglob('*.py'))}
    identity=signature(dict(data=data,code=code,chromosomes=a.chromosomes,tile_bp=a.tile_bp))
    record=a.output/'run_identity.json'
    if record.exists() and read_json(record)['identity']!=identity:raise ValueError('Output has other inputs/code; use a new output directory')
    write_json(record,dict(identity=identity,code=code))
    checked=preflight(data,a,a.output)
    lengths={c:checked['lengths'][c] for c in a.chromosomes.split(',')}
    geometry=tiles(lengths,a.tile_bp)
    aliases=fresh_indexes(checked['units'],checked['tools'],a.output,lambda:None)
    prepared={};done=[]
    for unit in checked['units']:
        if unit['sample'] not in prepared:
            prepared[unit['sample']]=prepare_tile_splits(unit,geometry,a.output,lambda:None)
        directory,hashes=prepared[unit['sample']]
        for chrom,start,end in geometry:
            name=f'{chrom}_{start}_{end}.tsv'
            query={role:aliases[unit['sample'],role] for role in ('phased_vcf','genotype_vcf')}
            task=dict(sample=unit['sample'],tissue=unit['tissue'],chrom=chrom,start0=start,end0=end,
                chrom_length=lengths[chrom],inputs=unit['inputs'],bam_device=unit['bam_device'],
                split_subset=str(directory/name),split_subset_sha256=hashes[name],
                query_paths={k:v['alias'] for k,v in query.items()},query_index_sha256={k:v['index_sha256'] for k,v in query.items()},
                tools=checked['tools'],write_detail=False,samtools_threads=0,max_seconds=21600,max_records=50000000,
                max_memory_gib=8,tile=True,adaptive_context=True,operational_context_halo=HALO,package_sha256=signature(code))
            task['task_id']=signature(task)
            out=a.output/'tiles'/unit['sample']/unit['tissue']/f'{chrom}_{start}_{end}'
            out.parent.mkdir(parents=True,exist_ok=True)
            if not complete_region(out,task['task_id']):
                if out.exists():raise ValueError('Incomplete tile exists; inspect it and use a fresh output directory: '+str(out))
                count_region(task,out)
            if not complete_region(out,task['task_id']):raise ValueError('Tile validation failed')
            done.append(str(out));write_json(a.output/'completed_tiles.json',done)
            print('Validated tiles:',len(done),flush=True)
            if a.stop_after_tiles and len(done)>=a.stop_after_tiles:return
    write_json(a.output/'completion.json',dict(status='PASS_ALL_REQUESTED_TILES',scope=lengths,units=keys,tiles=len(done)))

if __name__=='__main__':main()
