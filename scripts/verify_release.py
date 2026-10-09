"""Verify a fixed code archive without locating any biological data."""
from pathlib import Path
import argparse,hashlib,json

def sha(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda:f.read(1024**2),b''):h.update(block)
    return h.hexdigest()

def verify(root):
    root=root.resolve();manifest=json.loads((root/'release_manifest.json').read_text(encoding='utf-8'))
    for name,expected in manifest['sha256'].items():
        p=(root/name).resolve()
        if not p.is_relative_to(root):raise ValueError('Manifest path escapes release')
        if not p.is_file() or sha(p)!=expected:raise ValueError('Changed/missing release file: '+name)
    return dict(status='PASS',version=manifest['version'],files=len(manifest['sha256']))

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,default=Path(__file__).resolve().parents[1]);a=p.parse_args();print(json.dumps(verify(a.root)))
