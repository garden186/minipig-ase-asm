"""Validate registered depth outputs without relying on directory names."""
from pathlib import Path
import sys
from bootstrap import read_json, rows, sha, POLICY
from depth import summarize, HIST_FIELDS, OP_FIELDS


def tile(directory, source_hash, release_hash, audit_rows=False):
    gate = read_json(directory / 'validation.json')
    if (gate['status'] != 'PASS_GQ20_FRAGMENT_DEPTH' or gate['errors'] or
            gate['policy'] != POLICY or gate['source_gate_sha256'] != source_hash or
            gate.get('release_sha256') != release_hash):
        raise ValueError('Invalid depth checkpoint: ' + str(directory))
    for name, expected in gate['sha256'].items():
        if sha(directory / name) != expected:
            raise ValueError('Depth output checksum differs: ' + name)
    if audit_rows:
        hist, metrics, opp = summarize(rows(directory / 'cpg_depth_by_phase_set.tsv.gz'))
        observed_hist = {(r['population'], r['metric'], int(r['depth'])): int(r['observations'])
                         for r in rows(directory / 'depth_histogram.tsv')}
        observed_opp = [[r['unit']] + [int(r[k]) for k in OP_FIELDS[1:]]
                        for r in rows(directory / 'depth_opportunity.tsv')]
        if dict(hist) != observed_hist or opp != observed_opp:
            raise ValueError('Depth summaries do not reproduce serialized counts')
        if any(gate['summary'].get(k) != v for k, v in metrics.items()):
            raise ValueError('Depth metrics differ from serialized counts')
    return gate


def run(directory, audit_rows=False):
    gate = read_json(directory / 'validation.json')
    if gate['status'] not in ('PASS_GENOME_GQ20_DEPTH', 'PASS_BOUNDED_GQ20_DEPTH') or gate['errors']:
        raise ValueError('The depth run is not complete')
    request = read_json(directory / 'run_request.json')
    if len(gate['regions']) != len(request['selected_tasks']):
        raise ValueError('Depth run has incomplete geometry')
    expected_tasks = {r['task_id']: r for r in request['selected_tasks']}
    if len(expected_tasks) != len(request['selected_tasks']):
        raise ValueError('Duplicate requested tile')
    seen = set()
    for item in gate['regions']:
        task_id = item['task_id']
        if task_id in seen or task_id not in expected_tasks:
            raise ValueError('Duplicate or unexpected completed tile')
        seen.add(task_id)
        src = expected_tasks[task_id]
        destination = directory / item['directory']
        if sha(destination / 'validation.json') != item['validation_sha256']:
            raise ValueError('Depth tile registration differs')
        found = tile(destination, src['validation_sha256'], gate['release_sha256'], audit_rows)
        if found['task']['task_id'] != task_id:
            raise ValueError('Depth tile identity differs')
    for name, digest in gate['sha256'].items():
        if sha(directory / name) != digest:
            raise ValueError('Cohort output checksum differs: ' + name)
    return {'status': 'PASS_REGISTERED_GQ20_DEPTH_REVIEW', 'regions': len(seen),
            'row_audit': audit_rows, 'errors': []}


if __name__ == '__main__':
    print(run(Path(sys.argv[1]), '--row-audit' in sys.argv[2:]))
