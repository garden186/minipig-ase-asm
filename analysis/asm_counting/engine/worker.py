"""Run one bounded region or one chromosome in an isolated process."""
from __future__ import annotations
import csv
import gzip
from pathlib import Path
import sys
import time
from common import read_json, write_json, fingerprint, stream, sha, rss_bytes
from profile import CORE, engine, Accumulator, DETAIL_FIELDS, run_records
from validate import validate

def snv_context(task, fetch, links, variants=None):
    contexts = {}
    split = CORE.load_snpsplit_alleles(task['split_subset'], fetch)
    for role in ('genotype_vcf', 'phased_vcf'):
        sample_col = 9 + CORE.read_vcf_samples(query_path(task, role)).index(task['sample'])
        for fields in CORE.iter_vcf_records(task['tools']['bcftools'], query_path(task, role), fetch):
            if variants is not None:
                from evidence import variant_row
                variants.append(variant_row(role,fields,sample_col,links))
            alleles = [fields[3].upper()] + fields[4].upper().split(',')
            if any(a not in ('A', 'C', 'G', 'T') for a in alleles):
                continue
            p = int(fields[1]) - 1
            flags = contexts.setdefault(p, set())
            fmt = dict(zip(fields[8].split(':'), fields[sample_col].split(':')))
            gt = fmt.get('GT', '.')
            tokens = gt.replace('|', '/').split('/')
            try:
                called = [alleles[int(x)] for x in tokens]
            except (ValueError, IndexError):
                flags.add('MISSING_OR_INVALID_GENOTYPE')
                continue
            if len(called) != 2:
                flags.add('NON_DIPLOID_GENOTYPE')
            elif called[0] == called[1]:
                flags.add('HOMOZYGOUS_SNV')
            elif role == 'genotype_vcf':
                flags.add('GENOTYPE_HETEROZYGOUS')
            elif p in links:
                flags.add('ELIGIBLE_PHASED_LINK')
            elif '|' not in gt:
                flags.add('UNPHASED_HETEROZYGOUS')
            elif fmt.get('PS', '.') in ('', '.'):
                flags.add('PHASED_WITHOUT_PS')
            elif p not in split:
                flags.add('PHASED_NOT_IN_SNPSPLIT')
            elif tuple(called) != tuple(split[p]):
                flags.add('SNPSPLIT_ALLELE_ORDER_MISMATCH')
            elif frozenset(called) in CORE.CONVERSION_AMBIGUOUS:
                flags.add('CONVERSION_AMBIGUOUS_LINK_EXCLUDED')
            else:
                flags.add('PHASED_OTHER_NOT_SELECTED')
            if fields[6] not in ('PASS', '.'):
                flags.add('VCF_NONPASS_RECORD_PRESENT')
    for p in links:
        contexts.setdefault(p, set()).add('ELIGIBLE_PHASED_LINK')
    for flags in contexts.values():
        if flags == {'GENOTYPE_HETEROZYGOUS'}:
            flags.add('NO_USABLE_PHASE_RECORD')
    return contexts

def run_once(task, output):
    output = Path(output)
    output.mkdir(parents=True, exist_ok=False)
    started = time.monotonic()
    if any(fingerprint(v['path']) != v for v in task['inputs'].values()):
        raise ValueError('Input fingerprint changed before the region scan')
    if sha(task['split_subset']) != task['split_subset_sha256']:
        raise ValueError('Prepared SNP subset changed')
    for role, path in task.get('query_paths', {}).items():
        if fingerprint(path) != task['inputs'][role]: raise ValueError('VCF alias changed')
        if sha(path + '.csi') != task['query_index_sha256'][role]: raise ValueError('Derived VCF index changed')
    write_json(output / 'task.json', task)
    serial = [0]
    def logged(command):
        serial[0] += 1
        return stream(list(command), output / f'tool_{serial[0]:03d}.stderr.log', started + task['max_seconds'])
    CORE.stream_command = logged
    tools = task['tools']
    region = CORE.Region(task['chrom'], task['start0'], task['end0'])
    fetch = CORE.Region(region.contig, task.get('fetch_start0', max(0, region.start0 - 5000)), task.get('fetch_end0', min(task['chrom_length'], region.end0 + 5000)))
    reference_region = CORE.Region(region.contig, region.start0, min(task['chrom_length'], region.end0 + 1))
    positions, reference_qc = CORE.enumerate_shared_reference_cpgs(
        samtools=tools['samtools'], fasta=task['inputs']['reference']['path'],
        region=reference_region, excluded_het_positions0=[])
    positions = [p for p in positions if region.start0 <= p < region.end0]
    disrupted, genotype_qc = engine.helpers.disrupted_reference_cpgs(core=CORE, bcftools=tools['bcftools'],
        genotype_vcf=query_path(task, 'genotype_vcf'), sample=task['sample'], region=fetch, cpg_positions0=positions)
    _, links, phase_qc = CORE.load_variant_context(bcftools=tools['bcftools'],
        vcf_path=query_path(task, 'phased_vcf'), sample=task['sample'], snpsplit_path=task['split_subset'],
        extended_region=fetch, include_conversion_snps=False)
    variants = []
    contexts = snv_context(task, fetch, links, variants)
    from evidence import save_catalog
    save_catalog(output,variants,links,contexts)
    handle = None
    evidence_handle = None
    try:
        detail = None
        if task['write_detail']:
            handle = gzip.open(output / 'observation_detail.tsv.gz.partial', 'wt', encoding='utf-8', newline='')
            detail = csv.writer(handle, delimiter='\t', lineterminator='\n')
            detail.writerow(DETAIL_FIELDS)
        evidence_handle = gzip.open(output/'fragment_evidence.jsonl.gz.partial','wt',encoding='utf-8',compresslevel=1)
        acc = Accumulator(task['sample'], task['tissue'], task['chrom'], positions, disrupted, links,
                          contexts, fetch=fetch, detail=detail,evidence=evidence_handle)
        lines = logged([tools['samtools'], 'view', '-@', str(task['samtools_threads']), task['inputs']['bam']['path'], fetch.samtools])
        counting_seconds = run_records(lines, acc, task['max_records'], task['max_seconds'], task['max_memory_gib'])
        if task.get('adaptive_context') and acc.required_fetch_bounds != [fetch.start0, fetch.end0]:
            write_json(output/'context_expansion.json', {'required_bounds':acc.required_fetch_bounds})
            raise ContextExpansion(acc.required_fetch_bounds)
        evidence_handle.close()
        evidence_handle = None
        (output/'fragment_evidence.jsonl.gz.partial').rename(output/'fragment_evidence.jsonl.gz')
        if handle:
            handle.close()
            handle = None
            (output / 'observation_detail.tsv.gz.partial').rename(output / 'observation_detail.tsv.gz')
        acc.write(output)
        from reassign import write_asm_inputs
        write_asm_inputs(output)
        write_json(output / 'context_qc.json', {'reference': reference_qc, 'genotype': genotype_qc,
            'phase': phase_qc, 'known_snv_positions': len(contexts), 'eligible_link_snps': len(links),
            'core_region': region.samtools, 'fetch_region': fetch.samtools})
        del acc
        if any(fingerprint(v['path']) != v for v in task['inputs'].values()):
            raise ValueError('Input fingerprint changed during the region scan')
        gate = validate(output)
        write_json(output / 'timing.json', {'total_seconds': time.monotonic()-started, 'counting_seconds': counting_seconds,
            'peak_rss_bytes': rss_bytes(), 'output_bytes': sum(p.stat().st_size for p in output.iterdir() if p.is_file()),
            'status': gate['status']})
    finally:
        if handle:
            handle.close()
        if evidence_handle:
            evidence_handle.close()


class ContextExpansion(Exception):
    def __init__(self, bounds):
        self.bounds = bounds
        super().__init__('Additional mate or alignment context is required')


def query_path(task, role):
    return task.get('query_paths', {}).get(role, task['inputs'][role]['path'])


def run(task, output):
    if not task.get('adaptive_context'):
        return run_once(task, output)
    from genome_support import expanded_bounds
    output = Path(output)
    if output.exists(): raise FileExistsError(str(output))
    current = dict(task)
    for attempt in range(1, 9):
        candidate = output.with_name(output.name + '.context_' + str(attempt))
        # Existing failed attempts are never overwritten on a retry.
        while candidate.exists():
            attempt += 1
            candidate = output.with_name(output.name + '.context_' + str(attempt))
        try:
            run_once(current, candidate)
        except ContextExpansion as exc:
            current['fetch_start0'], current['fetch_end0'] = expanded_bounds(task, exc.bounds)
            continue
        # Move only our derived sibling directory after successful validation.
        if candidate.parent.resolve() != output.parent.resolve(): raise ValueError('Output parent mismatch')
        candidate.rename(output)
        return
    raise RuntimeError('Adaptive context retry limit reached; partial evidence is preserved')


if __name__ == '__main__':
    destination = Path(sys.argv[2])
    try:
        run(read_json(sys.argv[1]), destination)
    except BaseException as exc:
        if destination.exists():
            write_json(destination / 'failure.json', {'status': 'FAILED_OR_RESOURCE_STOPPED', 'error': str(exc)})
        raise
