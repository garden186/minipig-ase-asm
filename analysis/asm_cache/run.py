"""Build the original tested-CpG caches from explicitly registered annotated tables.

Only the final-analysis cache preparation kernel is invoked. The historical
exploratory inventory and statistical pilot are not run by this interface.
"""
import argparse
import re
from pathlib import Path
from cache import process_unit
from common import child, digest, lock, read, require, sha, verify_file, write


def execute(config_path, output, guards):
    config_path = Path(config_path).resolve()
    config = read(config_path)
    units = config['units']
    require(bool(units), 'No source units')
    prepared = []
    seen = set()
    output = Path(output).resolve()
    for unit in units:
        sample, tissue = unit['sample'], unit['tissue']
        require(all(isinstance(s, str) and re.fullmatch(r'[A-Za-z0-9_]+', s)
                    for s in (sample, tissue)), 'Unsafe/invalid sample or tissue ID')
        require((sample, tissue) not in seen, 'Duplicate sample/tissue')
        seen.add((sample, tissue))
        spec = dict(unit['annotated'])
        require((spec['sample'], spec['tissue']) == (sample, tissue), 'Source role mismatch')
        for key in ('path', 'gate'):
            if key in spec:
                p = Path(spec[key])
                spec[key] = str((config_path.parent / p).resolve())
                require(not Path(spec[key]).is_relative_to(output), 'Source inside output')
        require(Path(spec['path']).is_file(), 'Missing annotated table')
        require(re.fullmatch(r'[0-9a-f]{64}', spec['sha256']) is not None, 'Invalid input SHA256')
        require(type(spec['rows']) is int and spec['rows'] > 0, 'Invalid registered row count')
        if 'gate' in spec:
            verify_file(spec['gate'], spec['gate_sha256'])
            gate = read(spec['gate'])
            require(gate['status'].startswith('PASS') and not gate.get('errors'), 'Source gate failed')
            require(gate['sha256'][Path(spec['path']).name] == spec['sha256'], 'Source gate/table mismatch')
        prepared.append(dict(sample=sample, tissue=tissue, annotated=spec))
    code = {n: sha(Path(__file__).parent / n) for n in ('cache.py', 'common.py', 'run.py')}
    identity = digest(dict(units=prepared, code=code))
    output.mkdir(parents=True, exist_ok=True)
    with lock(output):
        request = output / 'run_request.json'
        if request.exists():
            require(read(request)['identity'] == identity, 'Output belongs to different inputs/code')
        else:
            write(request, dict(identity=identity, units=prepared, code=code))
        result = []
        for unit in prepared:
            spec = unit['annotated']
            dest = Path(process_unit((spec, str(child(output, 'units/' + unit['sample'] + '-' + unit['tissue'])),
                                      digest([identity, spec]), guards)))
            gate = dest / 'validation.json'
            result.append(dict(**unit, cache=dict(path=str(dest), gate=str(gate),
                gate_sha256=sha(gate), sha256=read(gate)['sha256'])))
        registry = dict(units=result, scope='Tested-CpG cache only; raw opportunities and WGS inputs are separate')
        write(output / 'cache_registry.json', registry)
        write(output / 'validation.json', dict(status='PASS_TESTED_CPG_CACHE_EXPORT', errors=[],
            identity=identity, units=len(result), rows=sum(u['annotated']['rows'] for u in result),
            sha256={'cache_registry.json': sha(output / 'cache_registry.json')}))
        return registry


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--config', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--reserve-memory-gib', type=float, default=2)
    p.add_argument('--reserve-disk-gib', type=float, default=10)
    p.add_argument('--worker-memory-gib', type=float, default=4)
    a = p.parse_args()
    require(a.reserve_memory_gib >= 0 and a.reserve_disk_gib >= 0 and a.worker_memory_gib > 0,
            'Invalid resource limits')
    result = execute(a.config, a.output, dict(reserve_memory_gib=a.reserve_memory_gib,
        reserve_disk_gib=a.reserve_disk_gib, worker_gib=a.worker_memory_gib))
    print('PASS_TESTED_CPG_CACHE_EXPORT units=' + str(len(result['units'])))


if __name__ == '__main__':
    main()
