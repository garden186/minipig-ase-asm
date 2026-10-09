"""Stage 1: full indexes and a bounded orientation pilot; no inferential statistics."""
from __future__ import annotations
import os
for variable in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS']:
    os.environ[variable]='1'
import argparse, collections, datetime, platform, sys, tarfile, traceback
import numpy as np
from common import *
from indexing import unit_job,scan_job,join_job
from pilot import pilot_job


def job(name,parent,identity,estimate,**kwargs):
    return dict(name=name,parent=str(parent),identity=digest([identity,name]),memory_estimate=int(estimate),**kwargs)

def execute(config,out,release_hash,workers=None,fixture=False):
    out=Path(out).resolve();out.mkdir(parents=True,exist_ok=True)
    require(config.get('fixture',False)==fixture,'Fixture/production configuration mismatch')
    identity=digest(dict(config=config,release=release_hash))
    ids=config['ids'];tissues=config['tissues'];chroms=config['chromosomes'];stages={};times={}
    def dispatch(fn,jobs,stage):
        if fixture:
            # The synthetic fixture has nine count rows; production estimates
            # would test host free RAM instead of the implementation.
            for j in jobs:j['memory_estimate']=64*1024**2
        return parallel(fn,jobs,out,stage,workers)
    if not fixture:
        require(ids==IDS and tissues==TISSUES and chroms==list(range(1,19)),'Unexpected production scope')
        require(len(config['units'])==30 and len(config['wgs'])==10,'Incomplete production inputs')
    with lock(out):
        if (out/'run_identity.json').exists():require(read(out/'run_identity.json')['identity']==identity,'Output belongs to another release/input set')
        else:write(out/'run_identity.json',dict(identity=identity,release=release_hash,config=config))
        guard(out)
        # Authentication is integrated into each large input scan. Resume verifies
        # the derived payloads; downstream stages use those authenticated payloads.
        reference=config['reference'];ref_before={}
        for path,h in [(reference['path'],reference['sha256']),(reference['index_path'],reference['index_sha256'])]:
            require(not Path(path).resolve().is_relative_to(out),'Source is inside output')
            ref_before[path]=verify(path,h)
        for spec in config['units']:
            for role in ['raw','annotated']:
                if 'gate' in spec[role]:verify(spec[role]['gate'],spec[role]['gate_sha256'])
            for path in [spec['raw']['path'],spec['cache']['path']]:require(not Path(path).resolve().is_relative_to(out),'Source inside output')
        for spec in config['wgs']:
            for role in ['phased','genotype','snpsplit']:
                require(not Path(spec[role]['path']).resolve().is_relative_to(out),'Source inside output')
                require(Path(spec[role]['path']).is_file(),'Missing source: '+spec[role]['path'])
        write(out/'machine.json',dict(python=sys.version,executable=sys.executable,platform=platform.platform(),
            numpy=np.__version__,psutil=psutil.__version__,resources=resources(),pid=os.getpid(),fixture=fixture,
            note='Worker count follows available CPU affinity/quota and memory admission. No timing guarantee.'))
        jobs=[]
        for spec in config['units']:
            name=spec['sample']+'-'+spec['tissue'];estimate=max(512*1024**2,spec['raw']['rows']*80+256*1024**2)
            jobs.append(job(name,out/'availability'/name,identity,estimate,spec=spec))
        started=time.monotonic();paths=dispatch(unit_job,jobs,'RAW_AVAILABILITY_AND_COUNT_CACHE')
        stages['units']={read(Path(p)/'summary.json')['sample']+'-'+read(Path(p)/'summary.json')['tissue']:p for p in paths}
        times['units']=time.monotonic()-started;write(out/'result_locations.json',stages)
        summaries=[read(Path(p)/'summary.json') for p in paths]
        require(sum(s['raw_rows'] for s in summaries)==config['expected']['raw_rows'],'Full raw total mismatch')
        require(sum(s['tested_rows'] for s in summaries)==config['expected']['tested_rows'],'Full tested total mismatch')
        jobs=[]
        for spec in config['wgs']:
            for role in ['phased','genotype','snpsplit']:
                name=spec['sample']+'-'+role
                jobs.append(job(name,out/'variant_sources'/name,identity,768*1024**2,spec=spec[role],sample=spec['sample'],role=role))
        started=time.monotonic();paths=dispatch(scan_job,jobs,'WGS_SOURCE_INDEXES')
        stages['variant_sources']={read(Path(p)/'summary.json')['sample']+'-'+read(Path(p)/'summary.json')['role']:p for p in paths}
        times['variant_sources']=time.monotonic()-started;write(out/'result_locations.json',stages)
        jobs=[]
        for s in ids:
            for c in chroms:
                name=f'{s}-chr{c}'
                jobs.append(job(name,out/'phase_index'/name,identity,512*1024**2,sample=s,chrom=c,reference=reference,
                    paths={r:stages['variant_sources'][s+'-'+r] for r in ['phased','genotype','snpsplit']}))
        started=time.monotonic();paths=dispatch(join_job,jobs,'EXACT_PHASE_ALLELE_INDEXES')
        stages['phase_index']={read(Path(p)/'summary.json')['sample']+'-chr'+str(read(Path(p)/'summary.json')['chrom']):p for p in paths}
        times['phase_index']=time.monotonic()-started;write(out/'result_locations.json',stages)
        jobs=[]
        for c in chroms:
            name=f'chr{c}'
            # Count/cache files are mapped, but include their sizes in the RAM
            # estimate so process parallelism remains conservative during scans.
            size=sum((Path(s['cache']['path'])/f'chr{c}.npy').stat().st_size for s in config['units'])
            estimate=max(GIB,3*size+512*1024**2)
            jobs.append(job(name,out/'pilot'/name,identity,estimate,chrom=c,ids=ids,tissues=tissues,units=config['units'],
                unit_paths=stages['units'],join_paths=stages['phase_index'],per_stratum=config['pilot_per_stratum'],
                frozen_anchors=str(asset_path(config['frozen_anchors'])) if not fixture else config['frozen_anchors'],
                frozen_profiles=str(asset_path(config['frozen_profiles'])) if not fixture else config['frozen_profiles']))
        started=time.monotonic();paths=dispatch(pilot_job,jobs,'OPPORTUNITY_PILOT_AND_FROZEN_REGRESSION')
        stages['pilot']={'chr'+str(read(Path(p)/'summary.json')['chrom']):p for p in paths}
        times['pilot']=time.monotonic()-started;write(out/'result_locations.json',stages)
        ps=[read(Path(p)/'summary.json') for p in paths]
        actual={k:sum(p[k] for p in ps) for k in ['full_eligible_coordinates','independent_sites','frozen_sites','multiple_PS_sites','pilot_sites','pilot_eligible_rows','regression_checked','regression_errors']}
        for key,want in config['expected'].items():
            if key in actual:require(actual[key]==want,'Full aggregation mismatch: '+key)
        states=collections.Counter()
        for p in stages['phase_index'].values():states.update(read(Path(p)/'summary.json')['states'])
        for path,snap in ref_before.items():require(snapshot(path)==snap,'Reference changed during run')
        summary=dict(status='PASS_STAGE1_INDEXES_AND_PILOT',identity=identity,counts=actual,
            source_raw_rows=config['expected']['raw_rows'],source_tested_rows=config['expected']['tested_rows'],
            full_variant_states=dict(states),timing_seconds=times,resources_at_end=resources(),
            alignment_scope='Pilot coordinates only. Full genome opportunity inventory and reusable phase indexes completed.',
            statistics='No new Fisher/partial-conjunction/BH values; no change to existing tissue OR inference.',
            next_step='Return receipt for review before implementing full genome orientation and directional recurrence.',fixture=fixture)
        write(out/'summary.json',summary)
        write(out/'validation.json',dict(status=summary['status'],errors=[],identity=identity,
            sha256={n:sha(out/n) for n in ['summary.json','result_locations.json','run_identity.json','machine.json']}))
        write(out/'progress.json',dict(status='COMPLETE_STAGE1',summary=summary))
        return summary

def export(out):
    """Return compact records; full indexes remain on the server."""
    out=Path(out).resolve();tag=datetime.datetime.now().strftime('%Y%m%d_%H%M%S')+'_'+str(time.time_ns())
    archive=out/f'SEND_THIS_FILE_ASM_allele_prep_{tag}.tar.gz'
    selected=[]
    for p in out.rglob('*'):
        if not p.is_file() or p.name.startswith('SEND_THIS_FILE') or p.name=='phase_dictionary.json':continue
        if p.suffix=='.json' or p.name in ['pilot_sites.tsv.gz','pilot_profiles.tsv.gz']:
            if 'selftest_work' not in p.parts:selected.append(p)
    # Old incomplete attempts are kept server-side; their summaries/error records
    # are useful diagnostics, but no large arrays or raw profiles are packaged.
    with tarfile.open(archive,'w:gz') as tar:
        for p in sorted(selected):tar.add(p,arcname='ASM_allele_prep_receipt/'+p.relative_to(out).as_posix(),recursive=False)
    checksum=archive.with_name(archive.name+'.sha256');checksum.write_text(sha(archive)+'  '+archive.name+'\n',encoding='ascii')
    print('SEND_THIS_FILE='+str(archive),flush=True);print('CHECKSUM='+str(checksum),flush=True)
    return archive

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('command',choices=['run','status','export'])
    p.add_argument('--output',required=True);p.add_argument('--workers',type=int)
    p.add_argument('--test-record');args=p.parse_args();out=Path(args.output).resolve()
    if args.command=='status':
        print(json.dumps(read(out/'progress.json') if (out/'progress.json').exists() else {'status':'NOT_STARTED'},indent=2));return
    if args.command=='export':export(out);return
    require(args.workers is None or args.workers>0,'Workers must be positive')
    out.mkdir(parents=True,exist_ok=True)
    # Verify release before any production scientific work.
    rh=release();config=read(input_config())
    require(args.test_record is not None,'Launcher selftest record required')
    record=read(args.test_record)
    require(record['status']=='PASS' and record['release_sha256']==rh and record['failures']==0 and record['errors']==0,'Selftests do not match this release')
    write(out/('selftest_'+str(time.time_ns())+'.json'),record)
    try:execute(config,out,rh,args.workers)
    except AlreadyRunning as e:raise SystemExit(str(e))
    except Exception as e:
        failure=dict(status='FAILED',error=repr(e),traceback=traceback.format_exc(),time=time.time(),
            note='Completed checkpoints retained; do not treat partial output as a passing result.')
        write(out/('failure_'+str(time.time_ns())+'.json'),failure);write(out/'progress.json',failure)
        export(out);raise
    export(out)

if __name__=='__main__':main()
