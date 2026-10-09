"""Full aligned paired-universe tissue OR recurrence; immutable authenticated input."""
import os
for key in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS']:os.environ[key]='1'
import argparse,datetime,platform,signal,sys,tarfile,traceback,shutil
import numpy as np
import scipy
from common import *
from data import PAIR_NAMES,rows,RESULT,INDIVIDUAL,FAMILY,FILTER
from pipeline import tile_job,family_job
from reporting import select,export_tile,merge_exports


def job(name,parent,identity,memory=GIB,**kw):
    return dict(name=name,parent=str(parent),identity=digest([identity,name]),memory_estimate=int(memory),**kw)

def execute(config,out,release_hash,workers=None,fixture=False,test_record=None):
    out=Path(out).resolve();out.mkdir(parents=True,exist_ok=True)
    require(config['fixture']==fixture,'Fixture/production mismatch')
    require(config['policy']==read(ROOT/'policy.json'),'Analysis policy differs from the released contract')
    runtime=dict(python=list(sys.version_info[:3]),numpy=np.__version__,scipy=scipy.__version__)
    identity=digest(dict(release=release_hash,config=config,fixture=fixture,runtime=runtime))
    with lock(out):
        guard(out,fixture=fixture)
        if (out/'run_identity.json').exists():require(read(out/'run_identity.json')['identity']==identity,'Different inputs/release; use a separate output directory')
        else:write(out/'run_identity.json',dict(identity=identity,release_sha256=release_hash,config=config,fixture=fixture,runtime=runtime))
        # A failed restart must not retain a stale COMPLETE gate while runtime
        # metadata has changed. Preserve the previous authenticated completion.
        if (out/'validation.json').exists():
            previous=read(out/'validation.json');history=out/'completion_history'/str(time.time_ns())
            for name,h in previous['sha256'].items():verify(out/name,h)
            history.mkdir(parents=True)
            for name in [*previous['sha256'],'validation.json']:shutil.copyfile(out/name,history/name)
            (out/'validation.json').unlink()
        write(out/'progress.json',dict(status='RUNNING',stage='VERIFY_AND_RESUME'))
        if test_record is not None:write(out/('selftest_'+str(time.time_ns())+'.json'),test_record)
        if (out/'STOP_REQUESTED').exists():(out/'STOP_REQUESTED').unlink()
        started=time.monotonic();snapshots={}
        for p,h in config['root_checks'].items():
            src=Path(p).resolve();require(not out.is_relative_to(src.parent) and not src.is_relative_to(out),'Output overlaps input')
            snapshots[p]=verify(p,h)
        if not fixture:
            gate=read(Path(config['stage2_root'])/'validation.json')
            require(gate['identity']==config['stage2_identity'] and gate['status']=='PASS_GENOME_DIRECTIONAL_ASM' and not gate['errors'],'Stage2 completion changed')
        write(out/'machine.json',dict(python=sys.version,numpy=np.__version__,scipy=scipy.__version__,psutil=psutil.__version__,platform=platform.platform(),resources=resources(),requested_workers=workers,fixture=fixture))
        write(out/'analysis_policy.json',config['policy'])
        if not fixture:
            expected_pairs=sum([3398278,3440630,3460303])
            array_bound=expected_pairs*(RESULT.itemsize+INDIVIDUAL.itemsize+9*FAMILY.itemsize)+config['expected']['sites']*FILTER.itemsize
            estimate=dict(known_CpG_pair_tests=expected_pairs,maximum_individual_tests=10*expected_pairs,
                retained_array_upper_bytes=int(array_bound),disk_reserve_bytes=20*GIB,
                note='Conservative N=10 array bound. Runtime and export volume depend on actual count depths and discoveries; native pilot timing is not a full-genome time guarantee.')
            write(out/'resource_plan.json',estimate)
            require(shutil.disk_usage(out).free>=20*GIB+array_bound,'Insufficient free disk for result arrays plus reserve')
        oldsignal=None
        if os.name=='posix':
            oldsignal=signal.getsignal(signal.SIGTERM)
            signal.signal(signal.SIGTERM,lambda *_:write(out/'STOP_REQUESTED',dict(time=time.time(),reason='SIGTERM')))
        try:
            def dispatch(fn,jobs,stage):
                if fixture:
                    for j in jobs:j.update(memory_estimate=64*1024**2,fixture=True)
                return sorted(parallel(fn,jobs,out,stage,workers,fixture=fixture))
            jobs=[]
            for spec in config['tiles']:
                source=Path(spec['path']).resolve();require(not out.is_relative_to(source) and not source.is_relative_to(out),'Output overlaps source tile')
                name=f"chr{spec['chrom']:02d}_{spec['start']:010d}_{spec['end']:010d}"
                memory=max(GIB,768*1024**2+spec['sites']*4000+spec['profiles']*200)
                jobs.append(job(name,out/'tiles'/name,identity,memory,spec=spec))
            tiles=dispatch(tile_job,jobs,'ALIGNED_PAIRED_EXACT_OR')
            totals=dict(tiles=len(tiles),sites=0,profiles=0,aligned=0);regression=0
            for path in tiles:
                s=read(Path(path)/'summary.json');totals['sites']+=s['input_sites'];totals['profiles']+=s['input_profiles'];totals['aligned']+=s['aligned_profiles'];regression+=s['regression_checks']
            require(totals==config['expected'],'Full aligned input inventory mismatch')
            if not fixture:require(regression==read(asset_path('assets/real_fixture_provenance.json'))['profiles'],'Full frozen-profile regression incomplete')
            jobs=[]
            for p in range(3):
                # Bound family sorting memory using the complete eligible universe.
                maximum=sum(sum(read(Path(t)/'summary.json')['per_pair'][PAIR_NAMES[p]]['available_by_N'][2:]) for t in tiles)
                for r in range(2,11):
                    name=f'{PAIR_NAMES[p]}_r{r}'
                    jobs.append(job(name,out/'families'/name,identity,max(GIB,maximum*220+512*1024**2),pair=p,r=r,tiles=tiles))
            families=dispatch(family_job,jobs,'COMPLETE_ALIGNED_FAMILY_BH_AND_BY')
            choice=select(job('selection',out/'selection',identity,tiles=tiles,families=families))
            summary=read(Path(choice)/'summary.json')
            if not fixture:
                expected=[3398278,3440630,3460303]
                for name,n in zip(PAIR_NAMES,expected):require(sum(summary['pair_availability'][name]['tested_by_N'][2:])==n,'Existing full paired N inventory differs')
            wanted=read(Path(choice)/'return_selection.json');exports=[]
            specs={s['path']:s for s in config['tiles']}
            for number,indices in wanted.items():
                tile=tiles[int(number)];spec=specs[read(Path(tile)/'summary.json')['source']]
                exports.append(job('export_'+number,out/'profiles'/('tile_'+number),identity,max(GIB,spec['profiles']*180),tile=tile,rows=indices,spec=spec))
            exports=dispatch(export_tile,exports,'EXPORT_COUNTS_EFFECTS_AND_RECURRENCE_LINKS')
            final,reused=begin(out/'final',digest([identity,'final']))
            if not reused:
                for p in Path(choice).iterdir():
                    if p.name.endswith(('.tsv','.tsv.gz')):shutil.copyfile(p,final/p.name)
                merge_exports(exports,final)
                write(final/'policy.json',config['policy'])
                summary.update(status='PASS_FULL_ALIGNED_TISSUE_ANALYSIS',input_counts=totals,frozen_real_profile_regression=regression,
                    fixture=fixture,elapsed_seconds=time.monotonic()-started,
                    inference_limits='Full-support exact OR and Fisher partial conjunction; signed fraction effects are descriptive. BH has cross-CpG assumptions; BY uses the same complete family. No sign filtering, union/max-r FDR or mean-DeltaM test. Raw count/phase/independent-animal assumptions remain.')
                seal(out/'final',final,digest([identity,'final']),summary)
            for p,before in snapshots.items():require(snapshot(p)==before,'Root source changed during execution')
            write(out/'result_locations.json',dict(tiles=tiles,families=families,selection=choice,profiles=exports,final=str(final)))
            write(out/'summary.json',read(final/'summary.json'))
            write(out/'validation.json',dict(status='PASS_FULL_ALIGNED_TISSUE_ANALYSIS',errors=[],identity=identity,
                sha256={p:sha(out/p) for p in ['summary.json','run_identity.json','result_locations.json','machine.json','analysis_policy.json']}))
            if (out/'STOP_REQUESTED').exists():raise StopRequested('Checkpoints complete; run again to finish')
            write(out/'progress.json',dict(status='COMPLETE',summary=read(final/'summary.json')))
            return read(final/'summary.json')
        finally:
            if oldsignal is not None:signal.signal(signal.SIGTERM,oldsignal)

def export(out):
    out=Path(out);rh=release()
    if (out/'run_identity.json').exists():require(read(out/'run_identity.json')['release_sha256']==rh,'Export package differs from the executed release')
    if (out/'validation.json').exists():
        gate=read(out/'validation.json')
        for name,h in gate['sha256'].items():verify(out/name,h)
        require(current(out/'final') is not None,'Final result checkpoint missing')
    tag=datetime.datetime.now().strftime('%Y%m%d_%H%M%S')+'_'+str(time.time_ns())
    target=out/('SEND_THIS_FILE_ASM_tissue_aligned_'+tag+'.tar.gz')
    with tarfile.open(target,'w:gz') as tar:
        for p in sorted(out.rglob('*')):
            if not p.is_file() or p.name.startswith('SEND_THIS_FILE'):continue
            rel=p.relative_to(out)
            include=p.suffix=='.json' or (rel.parts[0] in ['final','selection'] and p.name.endswith(('.tsv','.tsv.gz')))
            if include:tar.add(p,arcname='results/'+rel.as_posix(),recursive=False)
        for name in [*read(ROOT/'release.json')['sha256'],'release.json']:
            tar.add(ROOT/name,arcname='package/'+name,recursive=False)
    Path(str(target)+'.sha256').write_text(sha(target)+'  '+target.name+'\n',encoding='ascii')
    print('SEND_THIS_FILE='+str(target),flush=True);return target

def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('command',choices=['run','status','stop','export'])
    parser.add_argument('--output',required=True);parser.add_argument('--workers',type=int);parser.add_argument('--test-record')
    args=parser.parse_args();out=Path(args.output).resolve()
    if args.command=='status':print(json.dumps(read(out/'progress.json') if (out/'progress.json').exists() else dict(status='NOT_STARTED'),indent=2));return
    if args.command=='stop':require(out.is_dir(),'Output missing');write(out/'STOP_REQUESTED',dict(time=time.time()));print('Active jobs will checkpoint, then stop; run resumes.');return
    if args.command=='export':
        with lock(out):export(out)
        return
    require(args.workers is None or args.workers>0,'Worker count must be positive')
    rh=release();require(args.test_record is not None,'Exact-release selftest record required; use launcher')
    test=read(args.test_record);require(test['status']=='PASS' and test['release_sha256']==rh and test['failures']==test['errors']==test['skipped']==0,'Exact-release selftests failed/mismatched')
    out.mkdir(parents=True,exist_ok=True)
    try:
        execute(read(input_config()),out,rh,args.workers,test_record=test)
    except AlreadyRunning as e:raise SystemExit(str(e))
    except StopRequested as e:write(out/'progress.json',dict(status='STOPPED_RESUMABLE',reason=str(e)))
    except Exception as e:
        failure=dict(status='FAILED_RESUMABLE',error=repr(e),traceback=traceback.format_exc());write(out/('failure_'+str(time.time_ns())+'.json'),failure);write(out/'progress.json',failure)
        with lock(out):export(out)
        raise
    with lock(out):
        export(out)
if __name__=='__main__':main()
