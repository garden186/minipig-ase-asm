"""Bounded autosomal tiles, derived indexes, and resumable SNP preparation."""
from collections import Counter, OrderedDict
from bisect import bisect_right
import math
from pathlib import Path
import subprocess
import time
from common import fingerprint, sha, read_json, write_json, stamp, capture

HALO = 1_000_000


def tiles(lengths, width):
    if width < 100_000 or width > 10_000_000:
        raise ValueError('Tile width must be between 100 kb and 10 Mb')
    return [(chrom, start, min(start + width, lengths[chrom]))
            for chrom in sorted(lengths, key=int)
            for start in range(0, lengths[chrom], width)]


def expanded_bounds(task, required):
    left = max(0, min(required[0] - 5000, task['start0'] - 5000))
    right = min(task['chrom_length'], max(required[1] + 5000, task['end0'] + 5000))
    if left < max(0, task['start0'] - HALO) or right > min(task['chrom_length'], task['end0'] + HALO):
        raise RuntimeError('Required mate context exceeds the 1 Mb operational halo; no completed tile is emitted')
    return left, right


def fresh_indexes(units, tools, output, guard):
    parent = output / 'derived_vcf_indexes'
    parent.mkdir(exist_ok=True)
    result = {}
    for unit in units:
        for role in ('phased_vcf', 'genotype_vcf'):
            key = unit['sample'], role
            if key in result:
                continue
            record = parent / (unit['sample'] + '.' + role + '.json')
            source = unit['inputs'][role]
            guard()
            if record.exists():
                saved = read_json(record)
                if saved['source'] != source or fingerprint(saved['alias']) != source or sha(saved['alias'] + '.csi') != saved['index_sha256']:
                    raise ValueError('Derived VCF index checkpoint changed')
            else:
                directory = parent / (unit['sample'] + '.' + role + '.' + stamp())
                directory.mkdir()
                alias = directory / 'input.vcf.gz'
                alias.symlink_to(source['path'])
                partial = directory / 'new.csi.partial'
                print('BUILD_DERIVED_INDEX=' + unit['sample'] + ':' + role, flush=True)
                with (directory / 'index.stderr.log').open('w') as error:
                    subprocess.run([tools['bcftools'], 'index', '-c', '-o', str(partial), str(alias)],
                                   stdout=subprocess.DEVNULL, stderr=error, check=True, timeout=7200)
                partial.rename(str(alias) + '.csi')
                count = capture([tools['bcftools'], 'index', '-n', str(alias)], directory / 'count.stderr.log', 60).strip()
                if not count.isdigit() or fingerprint(alias) != source:
                    raise ValueError('Derived VCF index verification failed')
                saved = dict(source=source, alias=str(alias), index_sha256=sha(str(alias) + '.csi'),
                             indexed_records=int(count), original_index_modified=False,
                             method='Fresh CSI built by sequentially indexing the original BGZF VCF through a derived symlink')
                write_json(record, saved)
            result[key] = saved
    return result


def prepare_tile_splits(unit, specs, output, guard):
    """Read each SNPsplit source once, writing disjoint-core padded query shards."""
    parent = output / 'prepared_tile_snps'
    parent.mkdir(exist_ok=True)
    checkpoint = parent / (unit['sample'] + '.json')
    source = unit['inputs']['snpsplit_snps']
    geometry = [list(row) for row in specs]
    if checkpoint.exists():
        saved = read_json(checkpoint)
        if saved['source'] != source or saved['geometry'] != geometry:
            raise ValueError('SNP preparation geometry or source changed')
        if any(sha(Path(saved['directory']) / name) != value for name, value in saved['sha256'].items()):
            raise ValueError('Prepared SNP shard changed')
        return Path(saved['directory']), saved['sha256']
    directory = parent / (unit['sample'] + '.' + stamp())
    directory.mkdir()
    bychrom = {}
    for chrom, start, end in specs:
        name = f'{chrom}_{start}_{end}.tsv'
        (directory / name).touch(exist_ok=False)
        bychrom.setdefault(chrom, []).append((start - HALO, end + HALO, name))
    search = {}
    for chrom, entries in bychrom.items():
        entries.sort()
        search[chrom] = ([r[0] for r in entries], entries, max(r[1]-r[0] for r in entries))
    handles = OrderedDict()
    try:
        with open(source['path'], encoding='utf-8') as source_handle:
            for i, line in enumerate(source_handle):
                if i % 100000 == 0:
                    guard()
                if not line.strip() or line.startswith('#'):
                    continue
                fields = line.rstrip('\n').split('\t')
                if len(fields) < 5:
                    raise ValueError('Malformed SNPsplit source row')
                pos = int(fields[2]) - 1
                starts, entries, span = search.get(fields[1], ([], [], 0))
                for left, right, name in entries[bisect_right(starts,pos-span):bisect_right(starts,pos)]:
                    if left <= pos < right:
                        if name not in handles:
                            if len(handles) >= 48:
                                handles.popitem(last=False)[1].close()
                            handles[name] = (directory / name).open('a', encoding='utf-8', newline='')
                        handles.move_to_end(name)
                        handles[name].write(line)
    finally:
        for handle in handles.values():
            handle.close()
    if fingerprint(source['path']) != source:
        raise ValueError('SNPsplit source changed during preparation')
    hashes = {p.name: sha(p) for p in directory.glob('*.tsv')}
    write_json(checkpoint, dict(source=source, geometry=geometry, directory=str(directory), sha256=hashes))
    return directory, hashes


def validate_geometry(tasks, lengths, units):
    expected_units = {(u['sample'], u['tissue']) for u in units}
    seen = {}
    for task in tasks:
        key = task['sample'], task['tissue'], task['chrom']
        seen.setdefault(key, []).append((task['start0'], task['end0']))
    if set(seen) != {(sample, tissue, chrom) for sample, tissue in expected_units for chrom in lengths}:
        raise ValueError('Missing or unexpected sample/tissue/chromosome scope')
    for key, spans in seen.items():
        end = 0
        for left, right in sorted(spans):
            if left != end or right <= left:
                raise ValueError('Tile overlap or gap')
            end = right
        if end != lengths[key[2]]:
            raise ValueError('Incomplete chromosome coverage')


class Guard:
    """Check live resources frequently without repeatedly walking old outputs."""
    def __init__(self, output, request, started):
        self.output, self.request, self.started = output, request, started
        self.last_size_check = float('-inf')

    def __call__(self, *ignored):
        import shutil
        from scheduler import host_resources
        now = time.monotonic()
        if now - self.started >= self.request['max_seconds']:
            raise TimeoutError('Session time budget reached; resume preserves completed tiles')
        if shutil.disk_usage(self.output).free < self.request['min_free_disk_bytes']:
            raise RuntimeError('Free disk reserve reached')
        available = host_resources()['memory_available_bytes']
        if available is not None and available < self.request['min_available_memory_bytes']:
            raise MemoryError('Available memory fell below the shared-server reserve')
        if now - self.last_size_check >= 60:
            # Do not charge source VCF symlink sizes as newly written output.
            used = sum(p.stat().st_size for p in self.output.rglob('*') if p.is_file() and not p.is_symlink())
            if used > self.request['max_output_bytes']:
                raise RuntimeError('Output size budget reached')
            self.last_size_check = now
