"""Read-only sources, authenticated checkpoints, bounded parallel resource use."""
from __future__ import annotations
import contextlib, hashlib, json, math, os, shutil, socket, time
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, wait, FIRST_COMPLETED
import psutil

ROOT = Path(__file__).resolve().parent
VERSION = 'asm_allele_prep_v0.1.0'
GIB = 1024**3
IDS = '0235 0242 0276 0309 0326 0349 0355 0377 0384 0385'.split()
TISSUES = ['Heart', 'Kidney', 'Liver']

class AlreadyRunning(RuntimeError):pass

def require(ok, msg):
    if not ok: raise ValueError(msg)

def read(p): return json.loads(Path(p).read_text(encoding='utf-8'))

def sha(p):
    h=hashlib.sha256()
    with open(p,'rb') as f:
        for b in iter(lambda:f.read(4*1024**2),b''): h.update(b)
    return h.hexdigest()

def digest(x): return hashlib.sha256(json.dumps(x,sort_keys=True,separators=(',',':')).encode()).hexdigest()
def snapshot(p):
    s=Path(p).stat(); return [s.st_size,s.st_mtime_ns]

def verify(p,h):
    s=snapshot(p); require(sha(p)==h and snapshot(p)==s,'Source identity mismatch: '+str(p)); return s

def write(p,x):
    p=Path(p);p.parent.mkdir(parents=True,exist_ok=True)
    tmp=p.with_name(p.name+'.writing.'+str(os.getpid()))
    tmp.write_text(json.dumps(x,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    for k in range(11):
        try: os.replace(tmp,p);return
        except PermissionError:
            if os.name!='nt' or k==10: raise
            time.sleep(.1)

def safe_child(root,name):
    root=Path(root).resolve();p=(root/name).resolve()
    require(p.is_relative_to(root),'Path escapes output: '+str(p));return p

def cgroup_value(name):
    """Read this process's unified cgroup; unavailable values remain unknown."""
    if os.name!='posix': return None
    try:
        rel=next(x.split(':',2)[2] for x in Path('/proc/self/cgroup').read_text().splitlines() if x.startswith('0::'))
        p=Path('/sys/fs/cgroup')/rel.lstrip('/')/name
        if not p.exists(): p=Path('/sys/fs/cgroup')/name
        return p.read_text().strip()
    except (OSError,StopIteration): return None

def resources():
    mem=psutil.virtual_memory(); total=mem.total;available=mem.available
    limit=cgroup_value('memory.max');used=cgroup_value('memory.current')
    if limit and limit!='max' and used:
        total=min(total,int(limit));available=min(available,max(0,int(limit)-int(used)))
    try: cpus=len(os.sched_getaffinity(0))
    except AttributeError:
        try: cpus=len(psutil.Process().cpu_affinity())
        except (AttributeError,psutil.Error): cpus=os.cpu_count() or 1
    quota=cgroup_value('cpu.max')
    if quota and not quota.startswith('max'):
        q,p=map(int,quota.split());cpus=min(cpus,max(1,math.ceil(q/p)))
    reserve=int(.02*total)+2*GIB
    return dict(cpus=max(1,cpus),total=total,available=available,reserve=reserve,
                usable=max(0,available-reserve),hostname=socket.gethostname())

def guard(out,minimum_disk=20*GIB):
    r=resources();require(r['available']>r['reserve'],'Memory reserve reached; rerun resumes checkpoints')
    require(shutil.disk_usage(out).free>=minimum_disk,'Disk reserve reached; checkpoints retained')

def release():
    m=read(ROOT/'release.json')
    for name,h in m['sha256'].items(): verify(safe_child(ROOT,name),h)
    return sha(ROOT/'release.json')

def current(parent,identity=None,verify_payload=True):
    p=Path(parent)/'current.json'
    if not p.exists(): return None
    r=read(p)
    if identity is not None: require(r['identity']==identity,'Different checkpoint inputs: '+str(parent))
    out=safe_child(parent,r['directory']);verify(out/'validation.json',r['validation_sha256'])
    v=read(out/'validation.json')
    require(v['identity']==r['identity'] and not v['errors'] and v['status'].startswith('PASS'),'Failed checkpoint')
    if verify_payload:
        for name,h in v['sha256'].items(): verify(safe_child(out,name),h)
    return out

def begin(parent,identity):
    old=current(parent,identity)
    if old: return old,True
    out=Path(parent)/('attempt_'+str(time.time_ns())+'_'+str(os.getpid()))
    out.mkdir(parents=True,exist_ok=False);return out,False

def seal(parent,out,identity,summary):
    write(out/'summary.json',summary)
    files=sorted(p for p in out.iterdir() if p.is_file() and p.name!='validation.json')
    write(out/'validation.json',dict(status='PASS_ALLELE_PREPARATION_STAGE',errors=[],identity=identity,
                                    sha256={p.name:sha(p) for p in files}))
    write(Path(parent)/'current.json',dict(identity=identity,directory=out.name,validation_sha256=sha(out/'validation.json')))
    return str(out)

@contextlib.contextmanager
def lock(out):
    p=Path(out)/'controller.lock'
    with p.open('a+b') as f:
        if os.name=='nt':
            import msvcrt
            if f.tell()==0:f.write(b' ');f.flush()
            f.seek(0)
            try:msvcrt.locking(f.fileno(),msvcrt.LK_NBLCK,1)
            except OSError as e:raise AlreadyRunning('Another controller holds the output lock') from e
        else:
            import fcntl
            try:fcntl.flock(f.fileno(),fcntl.LOCK_EX|fcntl.LOCK_NB)
            except OSError as e:raise AlreadyRunning('Another controller holds the output lock') from e
        try: yield
        finally:
            if os.name=='nt':f.seek(0);msvcrt.locking(f.fileno(),msvcrt.LK_UNLCK,1)
            else:fcntl.flock(f.fileno(),fcntl.LOCK_UN)

def parallel(fn,jobs,out,stage,workers=None):
    """Jobs carry conservative RAM estimates; capacity rechecked before launches."""
    jobs=list(jobs);r=resources();limit=min(len(jobs),r['cpus'],workers or r['cpus'])
    results=[];active={};pending=list(jobs);start=time.monotonic()
    if not jobs:return []
    require(limit>=1,'No worker capacity')
    with ProcessPoolExecutor(max_workers=limit) as pool:
        while pending or active:
            r=resources()
            # Host available already excludes resident jobs; pending reservation also
            # prevents scheduling a whole batch before workers have allocated RAM.
            committed=sum(v['memory_estimate'] for v in active.values())
            allowance=min(r['usable'],max(0,r['total']-r['reserve']-committed-psutil.Process().memory_info().rss))
            while pending and len(active)<limit:
                job=pending[0];need=job['memory_estimate']
                if need>allowance:break
                pending.pop(0);active[pool.submit(fn,job)]=job;allowance-=need
            if not active:
                raise RuntimeError('Insufficient available RAM for next job; close competing work and resume')
            done,_=wait(active,timeout=10,return_when=FIRST_COMPLETED)
            for future in done:
                job=active.pop(future);results.append(future.result())
                print('STAGE='+stage+' FINISHED='+str(len(results))+'/'+str(len(jobs))+' JOB='+job['name'],flush=True)
            write(Path(out)/'progress.json',dict(status='RUNNING',stage=stage,completed=len(results),total=len(jobs),
                    active=[j['name'] for j in active.values()],resources=resources(),elapsed_seconds=time.monotonic()-start))
    return results

# Portable release adapter. Data/fixtures remain in an external directory.
def input_config():
    value=os.environ.get('MINIPIG_ANALYSIS_CONFIG')
    require(bool(value),'Set MINIPIG_ANALYSIS_CONFIG to your private stage input JSON')
    path=Path(value).expanduser().resolve()
    require(path.is_file(),'Private stage configuration is missing')
    return path

_checked_assets={}
def asset_path(name):
    value=os.environ.get('MINIPIG_ASSETS')
    require(bool(value),'Set MINIPIG_ASSETS to the external original stage asset directory')
    path=safe_child(Path(value).expanduser().resolve(),name)
    manifest=read(ROOT/'assets_manifest.json')
    require(name in manifest,'Unregistered biological asset: '+name)
    current_snapshot=snapshot(path)
    if _checked_assets.get(str(path))!=current_snapshot:
        verify(path,manifest[name]);_checked_assets[str(path)]=current_snapshot
    return path
