"""Independently reconcile CpG, exact-PS, entity, and optional membership counts."""
from __future__ import annotations
from collections import Counter
from pathlib import Path
import sys
from common import CALLS, STATES, SITE_FIELDS as BASE_FIELDS, PHASE_FIELDS, rows, read_json, sha, write_json, write_table

def require(condition, message):
    if not condition:
        raise ValueError(message)

def validate(directory):
    directory = Path(directory)
    task = read_json(directory / 'task.json')
    fields = BASE_FIELDS + ['fragment_assigned_using_mate_only']
    expected_phase, site_counts, summaries, histogram = {}, {}, Counter(), Counter()
    last = None
    nsites = 0
    for r in rows(directory / 'cpg_assignment.tsv.gz'):
        require(list(r) == fields, 'CpG table schema mismatch')
        key = (r['sample'], r['tissue'], r['chrom'], int(r['cpg_pos1']))
        require(last is None or key > last, 'CpG rows must be strictly sorted and unique')
        require(key[:3] == (task['sample'], task['tissue'], task['chrom']) and task['start0'] <= key[3]-1 < task['end0'], 'CpG row outside its disjoint core')
        last = key
        require(r['genotype_disrupted'] in ('0', '1'), 'Invalid genotype mask')
        values = {f: int(r[f]) for f in fields[5:]}
        require(all(v >= 0 for v in values.values()), 'Negative count')
        if task.get('adaptive_context'):
            require(not any(values['FRAGMENT_INCOMPLETE_FETCH_CONTEXT_' + c] for c in CALLS), 'Incomplete fragment context cannot complete a genome tile')
        cov = {}
        for unit in ('READ', 'FRAGMENT'):
            totals = {call: sum(values[f'{unit}_{s}_{call}'] for s in STATES) for call in CALLS}
            assigned = {call: values[f'{unit}_ASSIGNED_{call}'] for call in CALLS}
            cov[unit] = sum(totals.values())
            for scope in ('ALL_REFERENCE', 'GENOTYPE_RETAINED'):
                if scope == 'GENOTYPE_RETAINED' and r['genotype_disrupted'] == '1':
                    continue
                summaries[(scope, unit, 'covered_cpg_sites')] += int(cov[unit] > 0)
                summaries[(scope, unit, 'sites_with_any_methylated_observation')] += int(totals['METHYLATED'] > 0)
                for call in CALLS:
                    summaries[(scope, unit, 'total_' + call)] += totals[call]
                    summaries[(scope, unit, 'assigned_' + call)] += assigned[call]
                valid = totals['METHYLATED'] + totals['UNMETHYLATED']
                avail = assigned['METHYLATED'] + assigned['UNMETHYLATED']
                if valid:
                    bucket = '0' if avail == 0 else '1' if avail == valid else '(0,0.1]' if avail * 10 <= valid else '(0.1,0.25]' if avail * 4 <= valid else '(0.25,0.5]' if avail * 2 <= valid else '(0.5,1)'
                    histogram[(scope, unit, bucket)] += 1
        require(cov['FRAGMENT'] <= cov['READ'], 'Fragment coverage exceeds read coverage')
        require(all(values[f'READ_{s}_MATE_CPG_CONFLICT'] == 0 for s in STATES), 'A single read cannot have a mate CpG conflict')
        require(all(values[f'READ_INCOMPLETE_FETCH_CONTEXT_{c}'] == 0 for c in CALLS), 'Read-self context must not depend on mate fetch')
        require(values['fragment_assigned_using_mate_only'] <= sum(values['FRAGMENT_ASSIGNED_' + c] for c in CALLS), 'Invalid mate rescue count')
        expected_phase[key] = [values[f'{unit}_ASSIGNED_{call}'] for unit in ('READ', 'FRAGMENT') for call in CALLS]
        site_counts[key] = values if (directory / 'observation_detail.tsv.gz').exists() else None
        nsites += 1
    last = None
    nphase = 0
    for r in rows(directory / 'cpg_haplotype_by_ps.tsv.gz'):
        require(list(r) == PHASE_FIELDS, 'Phase table schema mismatch')
        key = (r['sample'], r['tissue'], r['chrom'], int(r['cpg_pos1']))
        ordered = key + (r['ps'],)
        require(last is None or ordered > last, 'Phase rows must be strictly sorted and unique')
        require(key in expected_phase and r['ps'] not in ('', '.', 'NA'), 'Unknown site or missing PS')
        last = ordered
        for ui, unit in enumerate(('READ', 'FRAGMENT')):
            for ci, call in enumerate(CALLS):
                a, b = int(r[f'{unit}_H1_{call}']), int(r[f'{unit}_H2_{call}'])
                require(a >= 0 and b >= 0, 'Negative haplotype count')
                expected_phase[key][ui * 4 + ci] -= a + b
        nphase += 1
    require(all(all(n == 0 for n in v) for v in expected_phase.values()), 'Exact-PS counts do not reconcile with assigned totals')
    source = read_json(directory / 'count_summary.json')
    require(source['covered_reference_cpgs'] == nsites and source['phase_rows'] == nphase, 'Saved row count mismatch')
    for unit in ('READ', 'FRAGMENT'):
        total = sum(summaries[('ALL_REFERENCE', unit, 'total_' + call)] for call in CALLS)
        require(total == source['qc'].get(unit + '_cpg_observations', 0), 'Observation count does not match source QC')
    entity_totals = Counter()
    for r in rows(directory / 'entity_assignment_context.tsv'):
        require(r['unit'] in ('READ', 'FRAGMENT') and r['state'] in STATES, 'Invalid entity state')
        require(int(r['n_entities']) > 0, 'Invalid entity count')
        entity_totals[r['unit']] += int(r['n_entities'])
    for unit in ('READ', 'FRAGMENT'):
        require(entity_totals[unit] == source['qc'].get(unit + '_entities_covering_reference_cpg', 0), 'Entity count mismatch')
    detail_counts = Counter()
    detail_rows = 0
    if (directory / 'observation_detail.tsv.gz').exists():
        identity = read_json(directory / 'task.json')
        for r in rows(directory / 'observation_detail.tsv.gz'):
            require(r['state'] in STATES and r['unit'] in ('READ', 'FRAGMENT'), 'Invalid membership state')
            require((r['state'] == 'ASSIGNED') == (r['haplotype'] in ('H1', 'H2') and r['ps'] != '.'), 'Membership haplotype/PS mismatch')
            for item in r['cpg_positions1_and_calls'].split(','):
                pos, call = item.split(':')
                require(call in CALLS, 'Invalid membership CpG call')
                key = (identity['sample'], identity['tissue'], identity['chrom'], int(pos))
                require(key in site_counts, 'Membership refers to an absent CpG')
                name = r['unit'] + '_' + r['state'] + '_' + call
                site_counts[key][name] -= 1
            detail_counts[r['unit']] += 1
            detail_rows += 1
        require(all(all(v == 0 for k, v in values.items() if k != 'fragment_assigned_using_mate_only') for values in site_counts.values()), 'Membership does not reconstruct CpG counts')
        require(detail_counts == entity_totals, 'Membership does not reconstruct entity counts')
    write_table(directory / 'measurement_summary.tsv', ['reference_scope', 'unit', 'metric', 'count'],
                (list(k) + [n] for k, n in sorted(summaries.items())))
    write_table(directory / 'assignment_fraction_histogram.tsv', ['reference_scope', 'unit', 'assigned_fraction_bin', 'n_cpg_sites'],
                (list(k) + [n] for k, n in sorted(histogram.items())))
    files = ['cpg_assignment.tsv.gz', 'cpg_haplotype_by_ps.tsv.gz', 'entity_assignment_context.tsv',
             'count_summary.json', 'task.json', 'measurement_summary.tsv', 'assignment_fraction_histogram.tsv']
    if detail_rows or (directory / 'observation_detail.tsv.gz').exists():
        files.append('observation_detail.tsv.gz')
    if (directory / 'context_qc.json').exists():
        files.append('context_qc.json')
    if (directory/'link_catalog.json').exists():
        from evidence import verify_evidence
        verify_evidence(directory)
        from check_asm_inputs import check
        check(directory)
        files.extend(['link_catalog.json','variant_records.tsv.gz','fragment_evidence.jsonl.gz','evidence_validation.json','cpg_asm_inputs.tsv.gz'])
    result = {'status': 'PASS_DESCRIPTIVE_CPG_ASSIGNMENT', 'errors': [], 'covered_cpg_rows': nsites,
              'exact_ps_rows': nphase, 'membership_rows_reconciled': detail_rows,
              'sha256': {name: sha(directory / name) for name in files},
              'limits': 'Arithmetic and membership reconciliation; no independent biological haplotype truth or sampling-bias validation'}
    write_json(directory / 'validation.json', result)
    return result

if __name__ == '__main__':
    try:
        result = validate(sys.argv[1])
        print(result['status'])
    except Exception as exc:
        write_json(Path(sys.argv[1]) / 'validation_failure.json', {'status': 'FAIL', 'errors': [str(exc)]})
        raise
