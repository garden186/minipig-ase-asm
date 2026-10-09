"""Provenance and I/O for a descriptive CpG assignment audit."""
from __future__ import annotations
import csv
import datetime as dt
import gzip
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import time

PACKAGE = Path(__file__).resolve().parent
ANIMALS = ('0235', '0242', '0276', '0309', '0326', '0349', '0355', '0377', '0384', '0385')
TISSUES = ('Heart', 'Kidney', 'Liver')
ROLES = ('bam', 'phased_vcf', 'genotype_vcf', 'snpsplit_snps', 'reference')
CALLS = ('METHYLATED', 'UNMETHYLATED', 'UNCALLABLE', 'MATE_CPG_CONFLICT')
STATES = ('ASSIGNED', 'NO_KNOWN_SNV_OVERLAP', 'NO_ELIGIBLE_PHASED_LINK',
          'NO_USABLE_LINK_BASE', 'NO_ALLELE_SUPPORT', 'MATE_SNP_CONFLICT',
          'CONFLICTING_HAPLOTYPES', 'MULTIPLE_PHASE_SETS', 'INCOMPLETE_FETCH_CONTEXT')
SITE_FIELDS = ['sample', 'tissue', 'chrom', 'cpg_pos1', 'genotype_disrupted'] + [
    f'{unit}_{state}_{call}' for unit in ('READ', 'FRAGMENT') for state in STATES for call in CALLS]
PHASE_FIELDS = ['sample', 'tissue', 'chrom', 'cpg_pos1', 'ps'] + [
    f'{unit}_{hap}_{call}' for unit in ('READ', 'FRAGMENT') for hap in ('H1', 'H2') for call in CALLS]
CONTRACT = {
    'version': '0.3.0', 'assembly': 'Sscrofa11.1', 'autosomes': list(map(str, range(1, 19))),
    'target': 'All observed reference CpGs, including zero-methylation sites; no CpG4 or depth/effect selection',
    'coordinates': '1-based reference-C coordinate; GA-strand observations use aligned G minus one',
    'read_unit': 'One retained primary alignment; read-self assignment only',
    'fragment_unit': 'Unmodified upstream QNAME mate grouping, minimum MAPQ 0, pair wait 2000 bp',
    'assignment': 'Unmodified September conversion-aware engine; exact PS, minimum one eligible link SNP',
    'link_selection': 'Upstream phased GT/PS and SNPsplit allele-order agreement; global C/T and A/G exclusion',
    'known_snv': 'Observed SNV positions in supplied phased and genotype VCFs, not proof of all true variants',
    'genotype_disrupted': 'Upstream DeepVariant-based mask is retained as a column, not silently discarded',
    'uncallable': 'Aligned strand-equivalent reference CpG base exists but XM has no Z/z; not unmethylated',
    'conflict': 'Contradictory mate CpG calls are a separate count, not methylated or unmethylated',
    'PS': 'H1/H2 counts are stored only with exact within-animal PS; site ASSIGNED totals pool origin availability only',
    'limits': 'Descriptive assignment opportunity; no proof of unbiased selection or true haplotype, no ASM tests',
    'evidence': 'Every retained fragment covering a reference CpG, including ambiguous/unassigned fragments; raw SNP bases and qualities on both mates, strand-aware CpG atoms, and VCF genotype/GQ/PS provenance.',
    'GQ': 'Extraction uses no GQ cutoff. Numeric phased GQ is used for later filtering; an exact-allele and genotype-matched DeepVariant GQ is a recorded fallback only when phased GQ is unavailable.',
    'replay': 'No-cutoff evidence replay must reproduce all baseline read/fragment CpG and exact-PS haplotype counts before completion.',
}

def sha(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for b in iter(lambda: f.read(1024 * 1024), b''):
            h.update(b)
    return h.hexdigest()

def signature(obj):
    return hashlib.sha256(json.dumps(obj, sort_keys=True, ensure_ascii=True).encode()).hexdigest()

def write_json(path, obj):
    path = Path(path)
    tmp = path.with_name(path.name + '.partial')
    tmp.write_text(json.dumps(obj, ensure_ascii=False, indent=2, sort_keys=True) + '\n', encoding='utf-8')
    os.replace(tmp, path)

def read_json(path):
    return json.loads(Path(path).read_text(encoding='utf-8'))

def fingerprint(path):
    p = Path(path).resolve(strict=True)
    s = p.stat()
    if not p.is_file() or s.st_size == 0:
        raise ValueError('Missing or empty input: ' + str(p))
    return {'path': str(p), 'size_bytes': s.st_size, 'mtime_ns': s.st_mtime_ns}

def write_table(path, fields, rows):
    path = Path(path)
    tmp = path.with_name(path.name + '.partial')
    opener = gzip.open if path.name.endswith('.gz') else open
    with opener(tmp, 'wt', encoding='utf-8', newline='') as f:
        w = csv.writer(f, delimiter='\t', lineterminator='\n')
        w.writerow(fields)
        w.writerows(rows)
    os.replace(tmp, path)

def rows(path):
    opener = gzip.open if str(path).endswith('.gz') else open
    with opener(path, 'rt', encoding='utf-8', newline='') as f:
        yield from csv.DictReader(f, delimiter='\t')

def stream(command, error_path, deadline=None):
    with open(error_path, 'w', encoding='utf-8') as err:
        p = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=err, text=True,
                             encoding='utf-8', bufsize=1024 * 1024)
        exhausted = False
        try:
            for line in p.stdout:
                if deadline and time.monotonic() > deadline:
                    raise TimeoutError('Execution time budget exceeded')
                yield line
            exhausted = True
        finally:
            p.stdout.close()
            if not exhausted and p.poll() is None:
                p.terminate()
            try:
                rc = p.wait(timeout=10)
            except subprocess.TimeoutExpired:
                p.kill()
                rc = p.wait()
            if exhausted and rc:
                raise RuntimeError('External command failed; inspect ' + str(error_path))

def capture(command, error_path, timeout=120):
    with open(error_path, 'w', encoding='utf-8') as err:
        p = subprocess.run(command, stdout=subprocess.PIPE, stderr=err, text=True,
                           encoding='utf-8', timeout=timeout, check=False)
    if p.returncode:
        raise RuntimeError('External command failed; inspect ' + str(error_path))
    return p.stdout

def rss_bytes():
    try:
        import resource
        return int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss) * (1 if sys.platform == 'darwin' else 1024)
    except ImportError:
        return None

def stamp():
    return dt.datetime.now(dt.timezone(dt.timedelta(hours=9))).strftime('%Y%m%d_%H%M%S_%f')

def verify_package():
    m = read_json(PACKAGE / 'package_manifest.json')
    for name, expected in m['sha256'].items():
        if sha(PACKAGE / name) != expected:
            raise ValueError('Package checksum mismatch: ' + name)
    return sha(PACKAGE / 'package_manifest.json')
