"""Shared provenance, restart, locking and resource rules."""
from contextlib import contextmanager
from pathlib import Path
from datetime import datetime, timezone
import hashlib
import json
import os
import shutil
import time

IDS = '0235 0242 0276 0309 0326 0349 0355 0377 0384 0385'.split()
TISSUES = ['Heart', 'Kidney', 'Liver']
ROOT = Path(__file__).resolve().parent
COUNTS = ['H1_methylated', 'H1_unmethylated', 'H2_methylated', 'H2_unmethylated']
QCOLS = ['BH_q_within_animal', 'BH_q_tissue10', 'BH_q_all30']
POLICY = dict(GQ=20, minimum_each=3, minimum_total=10, masks='original disruption masks',
              phase='exact chromosome/PS within animal; no count pooling',
              universe='all registered tested observations, no P/q/effect/feature prefilter',
              scope='preparation and statistical pilot, no genome-wide new significance',
              missing='absence from tested cache is UNTESTED_REASON_UNRESOLVED, not zero/negative',
              effects='Zelen tests OR homogeneity, not risk-difference equality',
              implementation='0.1.0')

def read(path):
    return json.loads(Path(path).read_text(encoding='utf-8'))

def sha(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for b in iter(lambda: f.read(4*1024**2), b''):
            h.update(b)
    return h.hexdigest()

def digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(',', ':')).encode()).hexdigest()

def write(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + '.writing')
    tmp.write_text(json.dumps(value, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')
    for retry in range(11):
        try:
            os.replace(tmp, path)
            break
        except PermissionError as exc:
            # Windows sync/indexing can briefly hold a just-written progress file.
            # Persistent permissions still fail; Linux publication stays unchanged.
            if os.name!='nt' or getattr(exc,'winerror',None) not in (5,32) or retry==10:
                raise
            time.sleep(.1)

def require(ok, message):
    if not ok:
        raise ValueError(message)

def child(root, name):
    root = Path(root).resolve()
    p = (root / name).resolve()
    require(p.is_relative_to(root), 'Unsafe registered path')
    return p

def stamp():
    return datetime.now(timezone.utc).strftime('%Y%m%d_%H%M%S') + '_' + str(time.time_ns())

def snapshot(path):
    s = Path(path).stat()
    return [s.st_size, s.st_mtime_ns]

def verify_file(path, expected):
    before = snapshot(path)
    require(sha(path) == expected and snapshot(path) == before, 'Input identity changed: ' + str(path))

def release():
    manifest = read(ROOT/'release_manifest.json')
    for n, h in manifest['sha256'].items():
        verify_file(child(ROOT, n), h)
    return sha(ROOT/'release_manifest.json')

def completed(parent, identity):
    pointer = Path(parent)/'current.json'
    if not pointer.exists():
        return None
    record = read(pointer)
    require(record['identity'] == identity, 'Existing checkpoint has different inputs/settings: ' + str(parent))
    dest = child(parent, record['directory'])
    verify_file(dest/'validation.json', record['validation_sha256'])
    gate = read(dest/'validation.json')
    require(gate['identity'] == identity and gate['status'].startswith('PASS') and not gate['errors'], 'Invalid checkpoint')
    for n,h in gate['sha256'].items():
        verify_file(child(dest,n),h)
    return dest

def attempt(parent):
    p = Path(parent)/('attempt_'+stamp())
    p.mkdir(parents=True, exist_ok=False)
    return p

def seal(parent, dest, identity, files, summary):
    write(dest/'summary.json', summary)
    files = sorted(set(files) | {'summary.json'})
    write(dest/'validation.json', dict(status='PASS_PREPARATION_STAGE', errors=[], identity=identity,
          sha256={n:sha(dest/n) for n in files}))
    write(Path(parent)/'current.json', dict(directory=dest.name, identity=identity,
          validation_sha256=sha(dest/'validation.json')))
    return str(dest)

class Guard:
    def __init__(self, output, reserve_memory_gib=32, reserve_disk_gib=50, worker_gib=4):
        self.output = Path(output)
        self.mem, self.disk, self.worker = (x*1024**3 for x in (reserve_memory_gib,reserve_disk_gib,worker_gib))
    def __call__(self):
        import psutil
        require(psutil.virtual_memory().available >= self.mem, 'Host memory reserve reached; checkpoints retained')
        require(shutil.disk_usage(self.output).free >= self.disk, 'Disk reserve reached; checkpoints retained')
        require(psutil.Process().memory_info().rss <= self.worker, 'Worker RSS budget exceeded; checkpoints retained')

@contextmanager
def lock(output):
    p = Path(output)/'controller.lock'
    with p.open('a+b') as f:
        if os.name == 'nt':
            import msvcrt
            if f.tell() == 0:
                f.write(b' '); f.flush()
            f.seek(0)
            try: msvcrt.locking(f.fileno(), msvcrt.LK_NBLCK, 1)
            except OSError: raise RuntimeError('Another controller owns this workspace')
        else:
            import fcntl
            try: fcntl.flock(f.fileno(), fcntl.LOCK_EX|fcntl.LOCK_NB)
            except BlockingIOError: raise RuntimeError('Another controller owns this workspace')
        try:
            yield
        except BaseException as exc:
            write(Path(output)/'progress.json',dict(status='FAILED',error=repr(exc),
                  note='Completed checkpoints retained; no source files or unrelated jobs modified'))
            raise
        finally:
            if os.name == 'nt': f.seek(0); msvcrt.locking(f.fileno(), msvcrt.LK_UNLCK,1)
            else: fcntl.flock(f.fileno(), fcntl.LOCK_UN)
