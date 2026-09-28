"""Exact depth distributions; physical CpG opportunity never pools phase sets."""
from collections import Counter, defaultdict
import math
from bootstrap import rows, write_table

COUNT_FIELDS = ['H1_methylated', 'H1_unmethylated', 'H2_methylated', 'H2_unmethylated']
FIELDS = ['sample', 'tissue', 'chrom', 'cpg_pos1', 'ps', 'genotype_disrupted'] + COUNT_FIELDS + [
    'H1_depth', 'H2_depth', 'baseline_H1_depth', 'baseline_H2_depth',
    'assignment_GQ_cutoff', 'observation_status']
HIST_FIELDS = ['population', 'metric', 'depth', 'observations']
OP_FIELDS = ['unit', 'minimum_depth_each_haplotype', 'baseline_observations', 'GQ20_observations']
DEPTHS = tuple(range(1, 101)) + (150, 200, 300, 500, 1000)


def summarize(table):
    """Build counters from serialized rows independently of fragment reassignment."""
    hist = Counter()
    metrics = Counter({k: 0 for k in (
        'all_phase_union_rows', 'masked_phase_union_rows', 'unmasked_phase_union_rows',
        'unmasked_physical_CpGs_in_phase_union', 'phase_rows_losing_all_valid_assignment',
        'phase_rows_gaining_valid_assignment')})
    for label in ('BASELINE', 'GQ20'):
        for suffix in ('physical_CpGs_with_any_valid_haplotype', 'valid_assigned_fragment_CpG_observations',
                       'phase_rows_any_valid', 'phase_rows_both_valid'):
            metrics[label + '_' + suffix] = 0
    joint = {'BASELINE': Counter(), 'GQ20': Counter()}
    physical = {'BASELINE': Counter(), 'GQ20': Counter()}
    last = None
    current = None
    best = [0, 0]
    any_valid = [False, False]

    def finish():
        if current is not None:
            for i, label in enumerate(('BASELINE', 'GQ20')):
                physical[label][best[i]] += 1
                metrics[label + '_physical_CpGs_with_any_valid_haplotype'] += any_valid[i]
            metrics['unmasked_physical_CpGs_in_phase_union'] += 1

    for r in table:
        if list(r) != FIELDS:
            raise ValueError('Unexpected depth row schema')
        key = (int(r['cpg_pos1']), r['ps'])
        if last is not None and key <= last:
            raise ValueError('Unsorted or duplicate CpG/phase-set row')
        last = key
        if r['ps'] in ('', '.') or r['genotype_disrupted'] not in ('0', '1'):
            raise ValueError('Missing phase set or invalid genotype mask')
        h = [int(r[k]) for k in COUNT_FIELDS]
        d = [int(r['H1_depth']), int(r['H2_depth'])]
        b = [int(r['baseline_H1_depth']), int(r['baseline_H2_depth'])]
        if min(h + d + b) < 0 or d != [sum(h[:2]), sum(h[2:])] or not sum(d + b):
            raise ValueError('Invalid fragment depths')
        if r['assignment_GQ_cutoff'] != '20':
            raise ValueError('Only GQ20 is supported')
        status = ('BOTH_HAPLOTYPES_OBSERVED' if min(d) > 0 else
                  'ONE_HAPLOTYPE_UNOBSERVED' if sum(d) else 'NO_VALID_HAPLOTYPE_AFTER_GQ20')
        if r['observation_status'] != status:
            raise ValueError('Depth observation status differs')
        metrics['all_phase_union_rows'] += 1
        if r['genotype_disrupted'] == '1':
            metrics['masked_phase_union_rows'] += 1
            continue
        if current != key[0]:
            finish()
            current, best, any_valid = key[0], [0, 0], [False, False]
        metrics['unmasked_phase_union_rows'] += 1
        for i, (label, values) in enumerate((('BASELINE', b), ('GQ20', d))):
            smaller = min(values)
            best[i] = max(best[i], smaller)
            any_valid[i] = any_valid[i] or bool(sum(values))
            joint[label][smaller] += 1
            metrics[label + '_valid_assigned_fragment_CpG_observations'] += sum(values)
            if sum(values):
                metrics[label + '_phase_rows_any_valid'] += 1
            if smaller:
                metrics[label + '_phase_rows_both_valid'] += 1
        for label, include in (('PHASE_UNION', True), ('GQ20_ANY_VALID', bool(sum(d))),
                               ('GQ20_BOTH_VALID', bool(min(d)))):
            if include:
                for name, value in zip(('H1', 'H2', 'SMALLER', 'LARGER', 'SUM'),
                                       (d[0], d[1], min(d), max(d), sum(d))):
                    hist[label, name, value] += 1
        if sum(b) and not sum(d):
            metrics['phase_rows_losing_all_valid_assignment'] += 1
        if not sum(b) and sum(d):
            metrics['phase_rows_gaining_valid_assignment'] += 1
    finish()
    opportunity = []
    for label, counts in (('CpG_PHASE_SET', joint), ('PHYSICAL_CpG', physical)):
        for threshold in DEPTHS:
            opportunity.append([label, threshold] + [sum(n for d, n in counts[g].items() if d >= threshold)
                                                     for g in ('BASELINE', 'GQ20')])
    return hist, metrics, opportunity


def dump(directory, result):
    hist, metrics, opportunity = result
    write_table(directory / 'depth_histogram.tsv', HIST_FIELDS,
                (list(k) + [v] for k, v in sorted(hist.items())))
    write_table(directory / 'depth_opportunity.tsv', OP_FIELDS, opportunity)
    return dict(metrics)


def quantiles(hist):
    """Nearest-rank empirical quantiles, including zero depth where appropriate."""
    n = sum(hist.values())
    if not n:
        return [0, '.', '.', '.', '.', '.', '.']
    ranks = [max(1, math.ceil(q * n)) for q in (.10, .25, .50, .75, .90, .95)]
    values, cumulative = [], 0
    for depth, count in sorted(hist.items()):
        cumulative += count
        while len(values) < len(ranks) and cumulative >= ranks[len(values)]:
            values.append(depth)
    return [n] + values


def aggregate(registry, output):
    """Sum disjoint tiles; preserve sample, tissue and chromosome provenance."""
    hist, opportunity, metrics, chromosomes, transitions, gqs = Counter(), Counter(), Counter(), Counter(), Counter(), Counter()
    for item in registry:
        directory = output / item['directory']
        from bootstrap import read_json
        gate = read_json(directory / 'validation.json')
        unit = (item['sample'], item['tissue'])
        for r in rows(directory / 'depth_histogram.tsv'):
            hist[unit + (r['population'], r['metric'], int(r['depth']))] += int(r['observations'])
        for r in rows(directory / 'depth_opportunity.tsv'):
            for label, field in (('BASELINE', 'baseline_observations'), ('GQ20', 'GQ20_observations')):
                opportunity[unit + (r['unit'], int(r['minimum_depth_each_haplotype']), label)] += int(r[field])
        for k, v in gate['summary'].items():
            metrics[unit + (k,)] += v
            chromosomes[unit + (item['chrom'], k)] += v
        for r in rows(directory / 'fragment_transitions.tsv'):
            transitions[unit + (r['baseline_state'], r['GQ20_state'])] += int(r['fragment_rows'])
        for r in rows(directory / 'core_link_GQ_distribution.tsv'):
            gqs[unit + (r['gq_source'], r['GQ'], r['retained_at_GQ20'])] += int(r['SNP_records'])
    write_table(output / 'sample_tissue_depth_histogram.tsv', ['sample', 'tissue'] + HIST_FIELDS,
                (list(k) + [v] for k, v in sorted(hist.items())))
    grouped = defaultdict(Counter)
    for item in registry:
        for population in ('PHASE_UNION', 'GQ20_ANY_VALID', 'GQ20_BOTH_VALID'):
            for metric in ('H1', 'H2', 'SMALLER', 'LARGER', 'SUM'):
                grouped[item['sample'], item['tissue'], population, metric]
    for k, v in hist.items():
        grouped[k[:-1]][k[-1]] += v
    write_table(output / 'sample_tissue_depth_quantiles.tsv',
                ['sample', 'tissue', 'population', 'metric', 'observations', 'p10', 'p25', 'median', 'p75', 'p90', 'p95'],
                (list(k) + quantiles(h) for k, h in sorted(grouped.items())))
    write_table(output / 'sample_tissue_depth_opportunity.tsv',
                ['sample', 'tissue', 'unit', 'minimum_depth_each_haplotype', 'GQ_setting', 'observations'],
                (list(k) + [v] for k, v in sorted(opportunity.items())))
    write_table(output / 'sample_tissue_summary.tsv', ['sample', 'tissue', 'metric', 'value'],
                (list(k) + [v] for k, v in sorted(metrics.items())))
    write_table(output / 'chromosome_summary.tsv', ['sample', 'tissue', 'chrom', 'metric', 'value'],
                (list(k) + [v] for k, v in sorted(chromosomes.items())))
    write_table(output / 'sample_tissue_fragment_transitions.tsv',
                ['sample', 'tissue', 'baseline_state', 'GQ20_state', 'fragment_rows'],
                (list(k) + [v] for k, v in sorted(transitions.items())))
    write_table(output / 'sample_tissue_link_GQ_distribution.tsv',
                ['sample', 'tissue', 'gq_source', 'GQ', 'retained_at_GQ20', 'SNP_records'],
                (list(k) + [v] for k, v in sorted(gqs.items())))
