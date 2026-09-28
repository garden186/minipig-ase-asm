"""Stage 2: genome-wide common-anchor alignment and directional recurrence."""
import os
for name in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS']:os.environ[name]='1'
import argparse,datetime,platform,sys,tarfile,signal,traceback
import numpy as np
import scipy
from common import *
from pipeline import align_tile,family_job
from reporting import assemble,annotation_job


def job(name,parent,identity,memory,**kw):return dict(name=name,parent=str(parent),identity=digest([identity,name]),memory_estimate=int(memory),**kw)

def authenticate_job(task):
    snaps={}
    for path,h in task['files'].items():snaps[path]=verify(path,h)
    return snaps

def authenticate(config,out,workers=None,fixture=False):
    root=Path(config['stage1_root']);verify(root/'validation.json',config['stage1_gate_sha256'])
    for name,h in config['stage1_root_payloads'].items():verify(root/name,h)
    require(read(root/'run_identity.json')['identity']==config['stage1_identity'],'Stage 1 identity differs')
    require(read(root/'summary.json')['fixture']==fixture,'Stage 1 fixture/production mismatch')
    require(not out.is_relative_to(root) and not root.is_relative_to(out),'Output overlaps stage1')
    tasks=[];counter=0
    for stage in ['units','phase_index','pilot','variant_sources']:
        for name,sp in config['locations'][stage].items():
            if stage=='variant_sources' and not name.endswith('-phased'):continue
            reg=config['registry'][sp];p=Path(sp);files={str(p/'validation.json'):reg['gate_sha256']}
            for n,h in reg['sha256'].items():
                if stage!='variant_sources' or n=='phase_dictionary.json':files[str(p/n)]=h
            tasks.append(dict(name='stage1_'+name,memory_estimate=64*1024**2,files=files));counter+=len(files)
    for u in config['source_inputs']['units']:
        cp=u['cache'];files={cp['gate']:cp['gate_sha256']}
        for n,h in cp['sha256'].items():files[str(Path(cp['path'])/n)]=h
        require(not out.is_relative_to(Path(cp['path'])) and not Path(cp['path']).is_relative_to(out),'Output overlaps original cache')
        tasks.append(dict(name='cache_'+u['sample']+'_'+u['tissue'],memory_estimate=64*1024**2,files=files));counter+=len(files)
    snaps={}
    for result in parallel(authenticate_job,tasks,out,'AUTHENTICATE_STAGE1_AND_ORIGINAL_COUNTS',workers,fixture=fixture):snaps.update(result)
    write(out/'source_audit.json',dict(files=len(snaps),snapshots=snaps,sha_pinned_in_inputs=True,fixture=fixture))
    return snaps

def execute(config,out,release_hash,workers=None,fixture=False):
    out=Path(out).resolve();out.mkdir(parents=True,exist_ok=True);identity=digest(dict(config=config,release=release_hash,fixture=fixture));locations={}
    with lock(out):
        old=out/'run_identity.json'
        if old.exists():require(read(old)['identity']==identity,'Output belongs to another release or inputs')
        else:write(old,dict(identity=identity,release=release_hash,config=config,fixture=fixture))
        if (out/'STOP_REQUESTED').exists():(out/'STOP_REQUESTED').unlink()
        started=time.monotonic();guard(out,fixture=fixture);write(out/'machine.json',dict(platform=platform.platform(),python=sys.version,numpy=np.__version__,scipy=scipy.__version__,psutil=psutil.__version__,pid=os.getpid(),resources=resources(),fixture=fixture))
        previous=None
        if os.name=='posix':
            previous=signal.getsignal(signal.SIGTERM)
            signal.signal(signal.SIGTERM,lambda *_:write(out/'STOP_REQUESTED',dict(reason='SIGTERM',time=time.time())))
        try:
            snaps=authenticate(config,out,workers,fixture)
            def dispatch(fn,jobs,stage):
                if fixture:
                    for j in jobs:j['memory_estimate']=64*1024**2;j['fixture']=True
                return parallel(fn,jobs,out,stage,workers,fixture=fixture)
            jobs=[]
            for name,sp in sorted(config['locations']['pilot'].items(),key=lambda kv:int(kv[0][3:])):
                c=int(name[3:]);pos=np.load(Path(sp)/'eligible_positions.npy',mmap_mode='r')
                if not len(pos):continue
                # No tested coordinate is omitted because of old significance.
                width=config['tile_bp']
                for start in range(1,int(pos[-1])+1,width):
                    end=start+width;lo,hi=np.searchsorted(pos,[start,end]);n=int(hi-lo)
                    if n==0:continue
                    label=f'chr{c}_{start}_{end}'
                    estimate=max(GIB,768*1024**2+n*14000)
                    jobs.append(job(label,out/'tiles'/label,identity,estimate,config=config,chrom=c,start=start,end=end))
            paths=dispatch(align_tile,jobs,'GENOME_ALIGNMENT_AND_DIRECTIONAL_PC')
            paths.sort(key=lambda p:(read(Path(p)/'summary.json')['chrom'],read(Path(p)/'summary.json')['start']))
            locations['tiles']=paths;write(out/'result_locations.json',locations)
            summaries=[read(Path(p)/'summary.json') for p in paths];expected=config['source_inputs']['expected']
            require(sum(x['sites'] for x in summaries)==expected['full_eligible_coordinates'],'Full CpG total changed')
            require(sum(x['profiles'] for x in summaries)==expected['tested_rows'],'Original eligible row total changed')
            require(sum(x['pilot_sites_checked'] for x in summaries)==config['expected_pilot_sites'],'Pilot anchor regression incomplete')
            require(sum(x['pilot_profiles_checked'] for x in summaries)==config['expected_pilot_profiles'],'Pilot profile regression incomplete')
            jobs=[]
            for ti,t in enumerate(TISSUES):
                for di,d in enumerate(['ALT_gt_REF','REF_gt_ALT']):
                    for r in range(2,11):
                        name=f'{t}_{d}_r{r}';m=sum(s['family_sizes'][t][str(r)] for s in summaries)
                        jobs.append(job(name,out/'families'/name,identity,max(512*1024**2,m*180),tiles=paths,tissue=ti,direction=di,r=r))
            fam=dispatch(family_job,jobs,'FULL_BACKGROUND_DIRECTIONAL_BH')
            fam.sort(key=lambda p:Path(p).parent.name);locations['families']=fam;write(out/'result_locations.json',locations)
            # Every direction has the same eligibility universe within tissue/r.
            for ti,t in enumerate(TISSUES):
                for r in range(2,11):
                    sizes=[read(Path(p)/'summary.json')['hypotheses'] for p in fam if read(Path(p)/'summary.json')['tissue']==t and read(Path(p)/'summary.json')['r']==r]
                    require(len(sizes)==2 and sizes[0]==sizes[1],'Direction-dependent denominator')
            assembled=assemble(job('assembly',out/'assembly',identity,GIB,config=config,tiles=paths,families=fam,fixture=fixture))
            locations['assembly']=assembled;write(out/'result_locations.json',locations)
            wanted=read(Path(assembled)/'selected_source_rows.json');jobs=[]
            for spec in config['source_inputs']['units']:
                s,t=spec['sample'],spec['tissue'];bychr={k.split('|')[2]:v for k,v in wanted.items() if k.startswith(s+'|'+t+'|')}
                if bychr:jobs.append(job(s+'-'+t,out/'annotations'/(s+'-'+t),identity,512*1024**2,spec=spec,wanted=bychr))
            locations['annotations']=dispatch(annotation_job,jobs,'ORIGINAL_SELECTED_ANNOTATIONS');write(out/'result_locations.json',locations)
            if (out/'STOP_REQUESTED').exists():raise StopRequested('Stopped after final stage checkpoint; run resumes')
            for p,before in snaps.items():require(snapshot(p)==before,'Source changed during computation: '+p)
            result=read(Path(assembled)/'summary.json')
            result.update(status='PASS_GENOME_DIRECTIONAL_ASM',identity=identity,full_eligible_CpGs=sum(x['sites'] for x in summaries),original_eligible_profiles=sum(x['profiles'] for x in summaries),
                tiles=len(paths),pilot_sites_replayed=sum(x['pilot_sites_checked'] for x in summaries),pilot_profiles_replayed=sum(x['pilot_profiles_checked'] for x in summaries),
                elapsed_seconds=time.monotonic()-started,fixture=fixture,scope='54 tissue/direction/r families; two-direction correction per tissue/r. Original tissue OR inference unchanged.')
            write(out/'summary.json',result);write(out/'validation.json',dict(status=result['status'],errors=[],identity=identity,sha256={n:sha(out/n) for n in ['summary.json','result_locations.json','run_identity.json','machine.json','source_audit.json']}))
            write(out/'progress.json',dict(status='COMPLETE_STAGE2',summary=result));return result
        finally:
            if previous is not None:signal.signal(signal.SIGTERM,previous)

def export(out):
    out=Path(out);name='SEND_THIS_FILE_ASM_directional_'+datetime.datetime.now().strftime('%Y%m%d_%H%M%S')+'_'+str(time.time_ns())+'.tar.gz';path=out/name
    with tarfile.open(path,'w:gz') as tar:
        for p in sorted(out.rglob('*')):
            if p.is_file() and not p.name.startswith('SEND_THIS_FILE') and (p.suffix=='.json' or p.name.endswith('.tsv.gz') or p.name in ['hits_BH10.npy','top25.npy']):
                tar.add(p,arcname='ASM_directional_receipt/'+p.relative_to(out).as_posix(),recursive=False)
    Path(str(path)+'.sha256').write_text(sha(path)+'  '+path.name+'\n',encoding='ascii');print('SEND_THIS_FILE='+str(path),flush=True);return path

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('command',choices=['run','status','stop','export']);p.add_argument('--output',required=True);p.add_argument('--workers',type=int);p.add_argument('--test-record');a=p.parse_args()
    out=Path(a.output).resolve()
    if a.command=='status':print(json.dumps(read(out/'progress.json') if (out/'progress.json').exists() else {'status':'NOT_STARTED'},indent=2));return
    if a.command=='stop':
        require(out.is_dir(),'Output does not exist');write(out/'STOP_REQUESTED',dict(time=time.time(),reason='User requested stop'));print('Stop requested. Active jobs will checkpoint before stopping.');return
    if a.command=='export':export(out);return
    require(a.workers is None or a.workers>0,'Workers must be positive');rh=release();test=read(a.test_record)
    require(test['status']=='PASS' and test['release_sha256']==rh and test['failures']==test['errors']==test['skipped']==0,'Exact release selftests required')
    out.mkdir(parents=True,exist_ok=True);write(out/('selftest_'+str(time.time_ns())+'.json'),test)
    try:execute(read(input_config()),out,rh,a.workers)
    except AlreadyRunning as e:raise SystemExit(str(e))
    except StopRequested as e:
        write(out/'progress.json',dict(status='STOPPED_RESUMABLE',reason=str(e)));export(out);return
    except Exception as e:
        failure=dict(status='FAILED_RESUMABLE',error=repr(e),traceback=traceback.format_exc(),time=time.time());write(out/('failure_'+str(time.time_ns())+'.json'),failure);write(out/'progress.json',failure);export(out);raise
    export(out)

if __name__=='__main__':main()
