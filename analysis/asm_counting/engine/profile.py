"""Count CpG coverage before assignment, using the pinned assignment engine."""
from __future__ import annotations
from array import array
import bisect
from collections import Counter
import csv
import gzip
from pathlib import Path
import sys
import time

from common import CALLS, STATES, SITE_FIELDS, PHASE_FIELDS, write_table, write_json, rss_bytes
sys.path.insert(0, str(Path(__file__).resolve().parent / 'upstream'))
import engine

CORE = engine.CORE
NCOUNT = len(STATES) * len(CALLS) * 2
SITE_FIELDS = SITE_FIELDS + ['fragment_assigned_using_mate_only']
DETAIL_FIELDS = ['unit', 'qname', 'alignment_flags', 'alignment_starts0', 'state', 'ps',
                 'haplotype', 'known_snv_positions1', 'eligible_link_positions1',
                 'variant_context', 'cpg_positions1_and_calls', 'minimum_mapq', 'link_base_quality_bin']

def overlap(record, positions):
    found = []
    for _q, start, length in CORE.aligned_blocks(record):
        found.extend(positions[bisect.bisect_left(positions, start):bisect.bisect_left(positions, start + length)])
    return found

def cpg_observations(record, positions):
    """Retain uncallable covered CpG bases as well as Z/z calls."""
    xg, xm = record.tags.get('XG'), record.tags.get('XM', '')
    if xg not in ('CT', 'GA') or record.seq == '*' or len(xm) != len(record.seq):
        raise ValueError('A usable alignment has missing or invalid XG/XM/SEQ tags')
    shift = int(xg == 'GA')
    out = {}
    for qstart, start, length in CORE.aligned_blocks(record):
        for cpos in positions[bisect.bisect_left(positions, start - shift):
                              bisect.bisect_left(positions, start + length - shift)]:
            qpos = qstart + cpos + shift - start
            symbol = xm[qpos]
            out[cpos] = 'METHYLATED' if symbol == 'Z' else 'UNMETHYLATED' if symbol == 'z' else 'UNCALLABLE'
    return out

def assign(records, links, positions, proxy, known_positions, contexts):
    known = sorted(set(p for record in records for p in overlap(record, known_positions)))
    eligible = sorted(set(p for record in records for p in overlap(record, positions)))
    calls, conflicts = {}, set()
    for record in records:
        for pos, base in proxy.extract_snp_calls(record, positions).items():
            CORE.merge_call(calls, conflicts, pos, base)
    states, votes = CORE.infer_ps_states(CORE.Fragment('audit', {}, calls), links, 1)
    ps, hap = '.', '.'
    if len(votes) > 1:
        state = 'MULTIPLE_PHASE_SETS'
    elif len(states) == 1:
        ps, hap = next(iter(states.items()))
        state = 'ASSIGNED'
    elif votes:
        state = 'CONFLICTING_HAPLOTYPES'
    elif conflicts:
        state = 'MATE_SNP_CONFLICT'
    elif not known:
        state = 'NO_KNOWN_SNV_OVERLAP'
    elif not eligible:
        state = 'NO_ELIGIBLE_PHASED_LINK'
    elif not calls:
        state = 'NO_USABLE_LINK_BASE'
    else:
        state = 'NO_ALLELE_SUPPORT'
    flags = sorted({flag for p in known for flag in contexts[p]})
    qualities = []
    for record in records:
        quality = getattr(record, 'audit_quality', '*')
        for qstart, start, length in CORE.aligned_blocks(record):
            for p in positions[bisect.bisect_left(positions, start):bisect.bisect_left(positions, start + length)]:
                q = qstart + p - start
                if quality != '*' and q < len(quality):
                    qualities.append(ord(quality[q]) - 33)
    minimum = min(qualities) if qualities else -1
    bq = 'MISSING' if minimum < 0 else 'BELOW_20' if minimum < 20 else '20_TO_29' if minimum < 30 else 'AT_LEAST_30'
    return {'state': state, 'ps': ps, 'hap': hap, 'known': known, 'eligible': eligible,
            'flags': flags, 'bq': bq, 'conflicting_snp_positions': len(conflicts)}

def merge_cpg(observations):
    result = {}
    for p in set().union(*(set(o) for o in observations)):
        called = {o[p] for o in observations if p in o and o[p] != 'UNCALLABLE'}
        result[p] = ('MATE_CPG_CONFLICT' if len(called) > 1 else next(iter(called)) if called else 'UNCALLABLE')
    return result

def count_bin(n):
    return '0' if n == 0 else '1' if n == 1 else '2_OR_MORE'

class Accumulator:
    def __init__(self, sample, tissue, chrom, positions, disrupted, links, contexts, fetch=None, detail=None, evidence=None):
        self.sample, self.tissue, self.chrom = sample, tissue, chrom
        self.positions, self.disrupted = positions, disrupted
        self.links, self.contexts = links, contexts
        self.link_positions, self.known_positions = sorted(links), sorted(contexts)
        self.qc, self.entities = Counter(), Counter()
        self.normalization_qc = Counter()
        self.proxy = engine.conversion_core(links, self.normalization_qc)
        self.sites, self.phases = {}, {}
        self.fetch, self.detail = fetch, detail
        self.evidence, self.evidence_rows = evidence, 0
        self.required_fetch_bounds = [fetch.start0, fetch.end0] if fetch else None

    def add(self, unit, records, result, observations):
        if not observations:
            return
        state = result['state']
        u = int(unit == 'FRAGMENT')
        self.qc[unit + '_entities_covering_reference_cpg'] += 1
        self.qc[unit + '_cpg_observations'] += len(observations)
        mapq = min(r.mapq for r in records)
        self.entities[(unit, state, count_bin(len(result['known'])), count_bin(len(result['eligible'])),
                       ','.join(result['flags']) or 'NONE', result['bq'], str(mapq),
                       ','.join(sorted({r.tags['XG'] for r in records})))] += 1
        for p, call in observations.items():
            v = self.sites.get(p)
            if v is None:
                v = self.sites[p] = array('Q', [0]) * (NCOUNT + 1)
            v[(u * len(STATES) + STATES.index(state)) * len(CALLS) + CALLS.index(call)] += 1
            if state == 'ASSIGNED':
                values = self.phases.get((p, result['ps']))
                if values is None:
                    values = self.phases[(p, result['ps'])] = array('Q', [0]) * 16
                offset = (u * 2 + int(result['hap'] == 'H2')) * 4 + CALLS.index(call)
                values[offset] += 1
        if self.detail:
            self.detail.writerow([unit, records[0].qname, ','.join(str(r.flag) for r in records),
                ','.join(str(r.pos0) for r in records), state, result['ps'], result['hap'],
                ','.join(str(p + 1) for p in result['known']) or '.',
                ','.join(str(p + 1) for p in result['eligible']) or '.', ','.join(result['flags']) or '.',
                ','.join(str(p + 1) + ':' + observations[p] for p in sorted(observations)), mapq, result['bq']])

    def consume(self, records):
        observations = [cpg_observations(r, self.positions) for r in records]
        if not any(observations):
            return
        if self.required_fetch_bounds is not None:
            for record in records:
                self.required_fetch_bounds[0] = min(self.required_fetch_bounds[0], record.pos0)
                self.required_fetch_bounds[1] = max(self.required_fetch_bounds[1], record.reference_end0)
                if (len(records) == 1 and record.flag & 1 and not record.flag & 8
                        and record.mate_rname in ('=', record.rname) and record.mate_pos0 >= 0):
                    self.required_fetch_bounds[0] = min(self.required_fetch_bounds[0], record.mate_pos0)
                    self.required_fetch_bounds[1] = max(self.required_fetch_bounds[1], record.mate_pos0 + 1)
        read_results = [assign((r,), self.links, self.link_positions, self.proxy,
                               self.known_positions, self.contexts) for r in records]
        for record, result, obs in zip(records, read_results, observations):
            self.add('READ', (record,), result, obs)
        result = assign(records, self.links, self.link_positions, self.proxy, self.known_positions, self.contexts)
        if self.fetch is not None and len(records) == 1:
            r = records[0]
            if (r.flag & 1 and not r.flag & 8 and r.mate_rname in ('=', r.rname)
                    and not self.fetch.start0 <= r.mate_pos0 < self.fetch.end0):
                result.update(state='INCOMPLETE_FETCH_CONTEXT', ps='.', hap='.')
        merged = merge_cpg(observations)
        if self.evidence is not None:
            from evidence import write_fragment
            self.evidence_rows += 1
            write_fragment(self.evidence,self.evidence_rows,records,self.positions,self.known_positions,result)
        self.add('FRAGMENT', records, result, merged)
        if result['state'] == 'ASSIGNED':
            for p in merged:
                covering = [rr for rr, obs in zip(read_results, observations) if p in obs]
                if covering and all(rr['state'] != 'ASSIGNED' for rr in covering):
                    self.sites[p][-1] += 1

    def write(self, directory):
        directory = Path(directory)
        prefix = [self.sample, self.tissue, self.chrom]
        mask = {p: int(self.disrupted[i]) for i, p in enumerate(self.positions)}
        write_table(directory / 'cpg_assignment.tsv.gz', SITE_FIELDS,
                    (prefix + [p + 1, mask[p]] + list(v) for p, v in sorted(self.sites.items())))
        write_table(directory / 'cpg_haplotype_by_ps.tsv.gz', PHASE_FIELDS,
                    (prefix + [p + 1, ps] + list(v) for (p, ps), v in sorted(self.phases.items())))
        write_table(directory / 'entity_assignment_context.tsv',
                    ['unit', 'state', 'known_snv_count_bin', 'eligible_link_count_bin', 'variant_context',
                     'minimum_link_base_quality_bin', 'minimum_mapq', 'XG', 'n_entities'],
                    (list(k) + [n] for k, n in sorted(self.entities.items())))
        write_json(directory / 'count_summary.json', {
            'qc': dict(self.qc), 'covered_reference_cpgs': len(self.sites),
            'conversion_normalizations_across_read_and_fragment_evaluations': dict(self.normalization_qc),
            'reference_cpgs_in_scope': len(self.positions), 'reference_cpgs_with_no_retained_alignment': len(self.positions)-len(self.sites),
            'phase_rows': len(self.phases), 'peak_rss_bytes': rss_bytes(),
            'scope': 'READ and FRAGMENT denominators are separate; entity counts differ from CpG observation counts'})

def run_records(records, accumulator, max_records, max_seconds, max_memory_gib):
    started = time.monotonic()
    qc, previous = Counter(), [-1]
    def checked():
        for line in records:
            if not line.strip() or line.startswith('@'):
                continue
            qc['raw_alignments_seen'] += 1
            if qc['raw_alignments_seen'] > max_records:
                raise RuntimeError('Alignment record budget exceeded')
            if time.monotonic() - started > max_seconds:
                raise TimeoutError('Counting time budget exceeded')
            if qc['raw_alignments_seen'] % 5000 == 0:
                memory = rss_bytes()
                if memory and memory > max_memory_gib * 1024**3:
                    raise MemoryError('Counting memory budget exceeded')
            record = CORE.parse_sam_line(line)
            record.audit_quality = line.rstrip('\n').split('\t')[10]
            if record.rname != accumulator.chrom or record.pos0 < previous[0]:
                raise ValueError('Unsorted stream or unexpected chromosome')
            previous[0] = record.pos0
            if CORE.record_is_usable(record, 0):
                if record.tags.get('XG') not in ('CT', 'GA') or len(record.tags.get('XM', '')) != len(record.seq):
                    raise ValueError('Missing or invalid bisulfite tags on a usable alignment')
            yield record
    for group in engine.helpers.stream_record_groups(core=CORE, records=checked(), min_mapq=0, pair_wait_bp=2000, qc=qc):
        accumulator.consume(group)
    accumulator.qc.update(qc)
    return time.monotonic() - started
