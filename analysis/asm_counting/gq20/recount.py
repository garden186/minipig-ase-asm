"""Reassign GQ20 fragment evidence without reading BAMs or changing CpG calls."""
from collections import Counter, defaultdict
import gzip
import json
from pathlib import Path
import sys
import time
from bootstrap import read_json, write_json, rows, write_table, sha, CALLS, STATES, POLICY
from evidence import selected_links, classify, methylation
from profile import merge_cpg
from depth import COUNT_FIELDS, FIELDS, summarize, dump

INPUTS = ('task.json', 'evidence_validation.json', 'link_catalog.json', 'fragment_evidence.jsonl.gz',
          'cpg_assignment.tsv.gz', 'cpg_asm_inputs.tsv.gz')


def source_check(source, expected_hash=None):
    gate_hash = sha(source / 'validation.json')
    if expected_hash and gate_hash != expected_hash:
        raise ValueError('Source tile registration checksum differs')
    gate = read_json(source / 'validation.json')
    if gate['status'] != 'PASS_DESCRIPTIVE_CPG_ASSIGNMENT' or gate['errors']:
        raise ValueError('Incomplete source tile')
    for name in INPUTS:
        if gate['sha256'].get(name) != sha(source / name):
            raise ValueError('Source checksum differs: ' + name)
    evidence_gate = read_json(source / 'evidence_validation.json')
    if evidence_gate['status'] != 'PASS_UNFILTERED_EVIDENCE_REPLAY' or evidence_gate['errors']:
        raise ValueError('Source evidence replay did not pass')
    return gate_hash


def count(source, output, expected_hash=None):
    started = time.monotonic()
    gate_hash = source_check(source, expected_hash)
    output.mkdir(parents=True, exist_ok=False)
    task = read_json(source / 'task.json')
    catalog = read_json(source / 'link_catalog.json')
    links = selected_links(catalog, 20)
    dropped = set(map(int, catalog['links'])) - links.keys()
    metadata, expected, base_assigned = {}, {}, {}
    for r in rows(source / 'cpg_assignment.tsv.gz'):
        pos = int(r['cpg_pos1']) - 1
        if pos in metadata or not task['start0'] <= pos < task['end0']:
            raise ValueError('Duplicate or out-of-core source CpG')
        if tuple(r[k] for k in ('sample', 'tissue', 'chrom')) != tuple(task[k] for k in ('sample', 'tissue', 'chrom')):
            raise ValueError('Source CpG identity differs')
        metadata[pos] = r['genotype_disrupted']
        expected[pos] = [sum(int(r[f'FRAGMENT_{state}_{call}']) for state in STATES) for call in CALLS]
        base_assigned[pos] = [int(r[f'FRAGMENT_ASSIGNED_{call}']) for call in CALLS[:2]]
    baseline = {}
    for r in rows(source / 'cpg_asm_inputs.tsv.gz'):
        pos, ps = int(r['cpg_pos1']) - 1, r['ps']
        if (pos, ps) in baseline or metadata.get(pos) != r['genotype_disrupted']:
            raise ValueError('Duplicate or mismatched baseline phase row')
        if tuple(r[k] for k in ('sample', 'tissue', 'chrom')) != tuple(task[k] for k in ('sample', 'tissue', 'chrom')):
            raise ValueError('Source phase identity differs')
        if r['assignment_GQ_cutoff'] != 'NONE':
            raise ValueError('Source extraction already has a GQ cutoff')
        baseline[pos, ps] = [int(r['H1_depth']), int(r['H2_depth'])]
    totals = defaultdict(lambda: [0, 0, 0, 0])
    after = defaultdict(lambda: [0, 0, 0, 0])
    transition = Counter()
    counts = Counter({k: 0 for k in ('fragment_rows', 'fragment_rows_reclassified',
        'fragment_rows_with_unchanged_link_evidence', 'covered_physical_CpGs',
        'masked_covered_physical_CpGs', 'unmasked_valid_fragment_CpG_observations',
        'unmasked_physical_CpGs_with_valid_calls')})
    with gzip.open(source / 'fragment_evidence.jsonl.gz', 'rt', encoding='utf-8') as f:
        for serial, line in enumerate(f, 1):
            e = json.loads(line)
            if e['schema_version'] != '0.2.0' or e['fragment_id'] != serial:
                raise ValueError('Invalid evidence identity')
            counts['fragment_rows'] += 1
            obs = merge_cpg([methylation(r) for r in e['reads']])
            old = e['baseline']
            # Identical linking inputs preserve the validated baseline assignment.
            if any(atom[0] in dropped for r in e['reads'] for atom in r['snps']):
                state, ps, hap = classify(e['reads'], links, catalog['contexts'], e['incomplete_fetch_context'])
                counts['fragment_rows_reclassified'] += 1
            else:
                state, ps, hap = old['state'], old['ps'], old['hap']
                counts['fragment_rows_with_unchanged_link_evidence'] += 1
            if state not in STATES or old['state'] not in STATES:
                raise ValueError('Unknown assignment state')
            if obs:
                transition[old['state'], state] += 1
            for pos, call in obs.items():
                if pos not in metadata:
                    raise ValueError('Evidence includes an unregistered CpG')
                index = CALLS.index(call)
                totals[pos][index] += 1
                if state == 'ASSIGNED' and index < 2:
                    if hap not in ('H1', 'H2') or ps in ('', '.'):
                        raise ValueError('Invalid assigned phase identity')
                    after[pos, ps][index + (2 if hap == 'H2' else 0)] += 1
    for pos, old in expected.items():
        if totals[pos] != old:
            raise ValueError('Fragment methylation totals changed at CpG ' + str(pos + 1))
    for pos, mask in metadata.items():
        counts['covered_physical_CpGs'] += 1
        counts['masked_covered_physical_CpGs'] += mask == '1'
        if mask == '0':
            valid = sum(expected[pos][:2])
            counts['unmasked_valid_fragment_CpG_observations'] += valid
            counts['unmasked_physical_CpGs_with_valid_calls'] += valid > 0
    per_site_after, per_site_before = Counter(), Counter()
    for (pos, ps), values in after.items():
        per_site_after[pos] += sum(values)
    for (pos, ps), values in baseline.items():
        per_site_before[pos] += sum(values)
    for pos in metadata:
        if per_site_before[pos] != sum(base_assigned[pos]) or per_site_after[pos] > sum(expected[pos][:2]):
            raise ValueError('Haplotype counts and CpG observation denominator disagree')

    def values():
        for pos, ps in sorted(baseline.keys() | after.keys()):
            h = after.get((pos, ps), [0] * 4)
            d = [sum(h[:2]), sum(h[2:])]
            b = baseline.get((pos, ps), [0, 0])
            status = ('BOTH_HAPLOTYPES_OBSERVED' if min(d) else
                      'ONE_HAPLOTYPE_UNOBSERVED' if sum(d) else 'NO_VALID_HAPLOTYPE_AFTER_GQ20')
            yield [task[k] for k in ('sample', 'tissue', 'chrom')] + [pos + 1, ps, metadata[pos]] + h + d + b + [20, status]

    write_table(output / 'cpg_depth_by_phase_set.tsv.gz', FIELDS, values())
    metrics = dump(output, summarize(rows(output / 'cpg_depth_by_phase_set.tsv.gz')))
    counts.update(metrics)
    counts['unmasked_valid_CpGs_without_GQ20_assignment'] = (
        counts['unmasked_physical_CpGs_with_valid_calls'] - counts['GQ20_physical_CpGs_with_any_valid_haplotype'])
    if counts['unmasked_valid_CpGs_without_GQ20_assignment'] < 0:
        raise ValueError('Invalid physical CpG denominator')
    write_table(output / 'fragment_transitions.tsv', ['baseline_state', 'GQ20_state', 'fragment_rows'],
                (list(k) + [v] for k, v in sorted(transition.items())))
    quality = Counter()
    for p, v in catalog['links'].items():
        if task['start0'] <= int(p) < task['end0']:
            quality[v['gq_source'], '.' if v['gq'] is None else str(v['gq']), int(int(p) in links)] += 1
    write_table(output / 'core_link_GQ_distribution.tsv', ['gq_source', 'GQ', 'retained_at_GQ20', 'SNP_records'],
                (list(k) + [v] for k, v in sorted(quality.items())))
    counts['core_eligible_SNP_records_before_GQ'] = sum(quality.values())
    counts['core_eligible_SNP_records_GQ20'] = sum(v for k, v in quality.items() if k[2])
    gate = dict(status='PASS_GQ20_FRAGMENT_DEPTH', errors=[], policy=POLICY,
                task={k: task[k] for k in ('sample', 'tissue', 'chrom', 'start0', 'end0', 'task_id')},
                source=str(source.resolve()), source_gate_sha256=gate_hash, summary=dict(counts),
                elapsed_seconds=time.monotonic() - started,
                sha256={p.name: sha(p) for p in output.iterdir() if p.is_file()},
                limits='GQ20 only; no ASM test, phase accuracy guarantee, or new base/mapping quality filter.')
    write_json(output / 'validation.json', gate)
    return gate


if __name__ == '__main__':
    job = read_json(Path(sys.argv[1]))
    count(Path(job['source']), Path(job['output']), job['validation_sha256'])
