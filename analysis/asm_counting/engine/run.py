"""Preflight and run a checkpointed 30-BAM descriptive assignment profile."""
from __future__ import annotations
import argparse
from collections import Counter
import datetime as dt
import os
from pathlib import Path
import re
import shutil
import sys
import tarfile
import time
from common import (ANIMALS, TISSUES, ROLES, CONTRACT, PACKAGE, fingerprint, sha, signature,
                    read_json, write_json, rows, write_table, capture, stamp, verify_package)

from scheduler import host_resources, concurrency_plan, schedule, resource_guard

def parser():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--mode', choices=('preflight', 'pilot'), default='pilot')
    p.add_argument('--manifest', type=Path, default=PACKAGE / 'sejong_input_manifest.tsv')
    p.add_argument('--units', help='Optional exact comma-separated sample-tissue units; default all 30')
    p.add_argument('--output-root', type=Path, required=True)
    p.add_argument('--resume', type=Path)
    p.add_argument('--tmp-root', type=Path, required=True)
    p.add_argument('--max-readers-per-device', type=int, choices=(1,2), default=1)
    p.add_argument('--min-available-memory-gib', type=float, default=32.0)
    p.add_argument('--min-free-disk-gib', type=float, default=50.0)
    p.add_argument('--pilot-result', type=Path, help='Completed identical-code and identical-input pilot required for full mode')
    p.add_argument('--max-hours', type=float, default=4.0)
    p.add_argument('--max-memory-gib', type=float, default=8.0)
    p.add_argument('--max-output-gib', type=float, default=20.0)
    p.add_argument('--jobs', type=int, choices=range(1,17), default=4, help='Requested independent BAM workers; capped by current CPU and available RAM')
    p.add_argument('--samtools-threads', type=int, choices=(0, 1), default=1)
    p.add_argument('--samtools', default='samtools')
    p.add_argument('--bcftools', default='bcftools')
    return p

def select(args):
    data = list(rows(args.manifest))
    keys = [r['sample'] + '-' + r['tissue'] for r in data]
    if len(data) != 30 or set(keys) != {a + '-' + t for a in ANIMALS for t in TISSUES}:
        raise ValueError('The manifest must contain the verified complete 10-by-3 roster exactly once')
    chosen = args.units.split(',') if args.units else keys
    if len(chosen) != len(set(chosen)) or set(chosen) - set(keys):
        raise ValueError('Unknown or duplicated requested unit')
    return [r for r in data if r['sample'] + '-' + r['tissue'] in chosen]

def index_for(path, bam=False):
    candidates = [path + s for s in (('.bai', '.csi') if bam else ('.tbi', '.csi'))]
    if bam:
        candidates += [str(Path(path).with_suffix('.bai'))]
    for p in candidates:
        if Path(p).is_file():
            return p
    raise ValueError('No existing index found for ' + path)

def preflight(data, args, output):
    missing = []
    for row in data:
        for role in ROLES:
            if not Path(row[role]).is_file():
                missing.append(row['sample']+'-'+row['tissue']+' '+role+': '+row[role])
    if missing:
        write_json(output/'preflight.json', {'status':'FAIL_PREFLIGHT', 'errors':missing, 'bam_contents_read':False})
        raise ValueError('Required inputs are missing; inspect preflight.json')
    tools, versions = {}, {}
    for name in ('samtools', 'bcftools'):
        executable = shutil.which(getattr(args, name))
        if not executable:
            raise ValueError('Required tool is unavailable: ' + name)
        tools[name] = str(Path(executable).resolve())
        versions[name] = capture([executable, '--version'], output / (name + '.stderr.log'), 30).splitlines()[:4]
    versions.update(python=sys.version, executable=sys.executable)
    reference_paths = {r['reference'] for r in data}
    if len(reference_paths) != 1:
        raise ValueError('A single explicit matching reference FASTA is required')
    reference = next(iter(reference_paths))
    lengths = {}
    with open(reference + '.fai', encoding='utf-8') as f:
        for line in f:
            a = line.split('\t')
            lengths[a[0]] = int(a[1])
    if any(str(c) not in lengths for c in range(1, 19)):
        raise ValueError('Reference must contain the recorded Sscrofa11.1 autosomes named 1 through 18')
    units, seen_headers, errors = [], set(), []
    for row in data:
        unit = row['sample'] + '-' + row['tissue']
        try:
            paths = {role: row[role] for role in ROLES}
            for role in ('bam', 'phased_vcf', 'genotype_vcf'):
                paths[role + '_index'] = index_for(row[role], role == 'bam')
            paths['reference_fai'] = reference + '.fai'
            inputs = {role: fingerprint(path) for role, path in paths.items()}
            if inputs['bam']['size_bytes'] != int(row['expected_bam_bytes']):
                raise ValueError('BAM size differs from the recorded source')
            capture([tools['samtools'], 'quickcheck', '-v', row['bam']], output / (unit + '.quickcheck.log'))
            header = capture([tools['samtools'], 'view', '-H', row['bam']], output / (unit + '.header.log'))
            (output / (unit + '.bam_header.sam')).write_text(header, encoding='utf-8')
            dictionary, sorted_bam = {}, False
            for line in header.splitlines():
                fields = line.split('\t')
                tags = dict(x.split(':', 1) for x in fields[1:] if ':' in x)
                if fields[0] == '@HD':
                    sorted_bam = tags.get('SO') == 'coordinate'
                if fields[0] == '@SQ':
                    dictionary[tags['SN']] = int(tags['LN'])
            if not sorted_bam or any(dictionary.get(str(c)) != lengths[str(c)] for c in range(1, 19)):
                raise ValueError('BAM sorting or reference chromosome lengths disagree')
            for role in ('phased_vcf', 'genotype_vcf'):
                key = (role, row[role], row['sample'])
                if key in seen_headers:
                    continue
                header = capture([tools['bcftools'], 'view', '-h', row[role]], output / (row['sample'] + '.' + role + '.log'))
                sample_names = next(line.split('\t')[9:] for line in header.splitlines() if line.startswith('#CHROM'))
                vcf_lengths = {c: int(n) for c, n in re.findall(r'##contig=<ID=([^,>]+),length=(\d+)', header)}
                if sample_names != [row['sample']] or any(vcf_lengths.get(str(c)) != lengths[str(c)] for c in range(1, 19)):
                    raise ValueError('VCF sample or reference chromosome lengths disagree: ' + role)
                seen_headers.add(key)
            units.append({'sample': row['sample'], 'tissue': row['tissue'], 'inputs': inputs,
                          'manifest_paths': paths, 'bam_device': Path(inputs['bam']['path']).stat().st_dev})
        except Exception as exc:
            errors.append(unit + ': ' + str(exc))
    result = {'status': 'PASS_PREFLIGHT' if not errors else 'FAIL_PREFLIGHT', 'errors': errors, 'units': units,
              'tools': tools, 'versions': versions, 'lengths': {str(c): lengths[str(c)] for c in range(1, 19)},
              'declared_assembly': CONTRACT['assembly'], 'genetic_identity_verified': False,
              'reference_note': 'Recorded reference path and matching chromosome dictionaries; not a whole-FASTA checksum or sample identity test'}
    write_json(output / 'preflight.json', result)
    if errors:
        raise ValueError('Input preflight failed; inspect preflight.json')
    return result

def prepare_split(unit, chroms, output, intervals=None):
    directory = output / 'prepared_snps' / unit['sample']
    source = unit['inputs']['snpsplit_snps']
    for candidate in sorted((output/'prepared_snps').glob(unit['sample']+'*')):
        record=candidate/'manifest.json'
        if record.exists():
            m=read_json(record)
            if (m['source'] != source or m.get('intervals') != intervals or
                    set(m['sha256']) != {'chr'+c+'.tsv' for c in chroms} or
                    any(sha(candidate/n) != s for n,s in m['sha256'].items())):
                raise ValueError('Prepared SNP source or cache changed')
            return candidate
    if directory.exists():
        directory=directory.with_name(directory.name+'_attempt_'+stamp())
    record=directory/'manifest.json'
    directory.mkdir(parents=True, exist_ok=False)
    handles = {c: (directory / ('chr' + c + '.tsv.partial')).open('w', encoding='utf-8', newline='') for c in chroms}
    counts = Counter()
    try:
        with open(source['path'], encoding='utf-8') as f:
            for line in f:
                if not line.strip() or line.startswith('#'):
                    continue
                fields = line.rstrip('\n').split('\t')
                if len(fields) < 5:
                    raise ValueError('Malformed source SNPsplit row')
                if (fields[1] in handles and (intervals is None or
                        intervals[fields[1]][0] <= int(fields[2])-1 < intervals[fields[1]][1])):
                    handles[fields[1]].write(line)
                    counts[fields[1]] += 1
    finally:
        for h in handles.values():
            h.close()
    for c in chroms:
        (directory / ('chr' + c + '.tsv.partial')).rename(directory / ('chr' + c + '.tsv'))
    if fingerprint(source['path']) != source:
        raise ValueError('SNPsplit source changed during preparation')
    write_json(record, {'source': source, 'counts': dict(counts), 'intervals': intervals,
                       'sha256': {p.name: sha(p) for p in directory.glob('*.tsv')}})
    return directory

def complete_region(path, task_id):
    p = Path(path)
    if not (p / 'validation.json').exists() or not (p / 'timing.json').exists():
        return False
    v = read_json(p / 'validation.json')
    if v['status'] != 'PASS_DESCRIPTIVE_CPG_ASSIGNMENT' or v['errors']:
        return False
    if read_json(p / 'task.json')['task_id'] != task_id:
        raise ValueError('Checkpoint task identity differs')
    if any(sha(p / n) != s for n, s in v['sha256'].items()):
        raise ValueError('A validated checkpoint was modified')
    return True

def forecast(pilot_path, request):
    gate = read_json(pilot_path / 'validation.json')
    old = read_json(pilot_path / 'run_request.json')
    if gate['status'] != 'PASS_COHORT_CPG_ASSIGNMENT_PILOT' or gate['errors']:
        raise ValueError('A completed pilot is required')
    if old['mode'] != 'pilot' or len(gate['regions']) != len(old['units'])*3:
        raise ValueError('Unexpected pilot scope or region count')
    if old['package_sha256'] != request['package_sha256'] or old['units'] != request['units'] or old['tools'] != request['tools'] or old['versions'] != request['versions'] or old['samtools_threads'] != request['samtools_threads']:
        raise ValueError('Pilot code, unit inputs, or tool paths differ')
    if any(sha(pilot_path/n) != s for n,s in gate['sha256'].items()):
        raise ValueError('Pilot aggregate changed')
    timings, sizes, densities = [], [], []
    for item in gate['regions']:
        directory = pilot_path / item['directory']
        if sha(directory/'validation.json') != item['validation_sha256'] or sha(directory/'timing.json') != item['timing_sha256']:
            raise ValueError('Pilot validation or timing record changed')
        task = read_json(directory / 'task.json')
        if not complete_region(directory, task['task_id']):
            raise ValueError('Pilot checkpoint is incomplete')
        t = read_json(directory / 'timing.json')
        timings.append(t['total_seconds'])
        sizes.append(t['output_bytes'])
        densities.append(read_json(directory / 'count_summary.json')['covered_reference_cpgs'] / (task['end0'] - task['start0']))
    unit_bp = len(request['units']) * sum(request['lengths'].values())
    pilot_bp = len(request['units']) * 300000
    scale = unit_bp / pilot_bp
    return {'status': 'FORECAST_ONLY', 'pilot_regions': len(timings),
            'full_autosome_unit_bp': unit_bp, 'pilot_unit_bp': pilot_bp,
            'forecast_serial_worker_seconds_with_factor_2': 2 * sum(timings) * scale,
            'forecast_seconds_with_factor_2': 2 * sum(timings) * scale / (min(request.get('jobs', 1), len(request['units'])) * 0.8),
            'forecast_jobs': min(request.get('jobs', 1), len(request['units'])),
            'assumed_parallel_efficiency': 0.8,
            'forecast_output_bytes_with_factor_2': 2 * sum(sizes) * scale,
            'forecast_largest_chromosome_memory_bytes': 512 * 1024**2 + max(densities, default=0) * max(request['lengths'].values()) * 2400,
            'limits': 'Linear opportunity extrapolation with factor 2 and assumed 80% parallel efficiency, not measured full-genome time/storage/RAM. Fixed pilot coordinates may misrepresent genome-wide density. The evidence-preserving pilot includes repeated indexed-query and validation overhead. Full execution is not enabled in this release. Preparation/preflight/final receipt time is not extrapolated.'}

def aggregate(output, request, completed):
    counts, hist, totals = Counter(), Counter(), Counter()
    registry = []
    for task, directory in completed:
        if not complete_region(directory, task['task_id']):
            raise ValueError('Final checkpoint reconciliation failed')
        unit = (task['sample'], task['tissue'])
        for r in rows(directory / 'measurement_summary.tsv'):
            counts[unit + (r['reference_scope'], r['unit'], r['metric'])] += int(r['count'])
        for r in rows(directory / 'assignment_fraction_histogram.tsv'):
            hist[unit + (r['reference_scope'], r['unit'], r['assigned_fraction_bin'])] += int(r['n_cpg_sites'])
        for k, n in read_json(directory / 'count_summary.json')['qc'].items():
            totals[unit + (k,)] += n
        registry.append({'directory': str(directory.relative_to(output)), 'validation_sha256': sha(directory / 'validation.json'),
                         'timing_sha256': sha(directory / 'timing.json')})
    if (output/'validation.json').exists():
        prior=read_json(output/'validation.json')
        if prior['regions'] != registry or any(sha(output/n) != s for n,s in prior['sha256'].items()):
            raise ValueError('A previously validated cohort result changed')
        return
    write_table(output / 'cohort_measurement_summary.tsv', ['sample','tissue','reference_scope','unit','metric','count'],
                (list(k)+[n] for k,n in sorted(counts.items())))
    write_table(output / 'cohort_assignment_fraction_histogram.tsv', ['sample','tissue','reference_scope','unit','assigned_fraction_bin','n_cpg_sites'],
                (list(k)+[n] for k,n in sorted(hist.items())))
    write_table(output / 'cohort_processing_qc.tsv', ['sample','tissue','metric','count'],
                (list(k)+[n] for k,n in sorted(totals.items())))
    write_json(output / 'validation.json', {'status': 'PASS_COHORT_CPG_ASSIGNMENT_' + request['mode'].upper(),
        'errors': [], 'units': len(request['units']), 'regions': registry,
        'raw_genome_scope': 'All 18 autosomes' if request['mode']=='full' else 'Fixed central 100 kb on chromosomes 1, 10, 18 per unit',
        'sha256': {p.name: sha(p) for p in output.glob('cohort_*.tsv')},
        'limits': 'Descriptive counts only. No selection-bias correction, random-sampling proof, ASM significance or whole-genome pilot extrapolation as an observed result.'})

def check_full_budget(plan,request):
    if (plan['forecast_seconds_with_factor_2'] > request['max_seconds'] or
            plan['forecast_output_bytes_with_factor_2'] > request['max_output_bytes'] or
            plan['forecast_largest_chromosome_memory_bytes'] > request['max_memory_gib']*1024**3):
        raise RuntimeError('Full forecast exceeds a declared resource budget; no full BAM scan started')

def main():
    args = parser().parse_args()
    if sys.flags.optimize or sys.version_info < (3,10):
        raise ValueError('Use Python >=3.10 without optimization flags')
    if args.max_hours <= 0 or args.max_memory_gib <= 0 or args.max_output_gib <= 0:
        raise ValueError('Resource budgets must be positive')
    from server_resources import configure_scratch
    scratch = configure_scratch(args.tmp_root)
    if args.min_available_memory_gib < 2 or args.min_free_disk_gib < 2:
        raise ValueError('Memory and disk reserves must be at least 2 GiB')
    package_sha = verify_package()
    if args.resume:
        output = args.resume.resolve(strict=True)
        request = read_json(output / 'run_request.json')
        if request['mode'] != 'pilot':
            raise ValueError('Only a pilot or full execution can be resumed')
        if request['mode']=='full':
            plan_path=output/'measured_full_forecast.json'
            if not plan_path.exists():
                raise ValueError('A full run without a completed resource forecast cannot be resumed')
            plan=read_json(plan_path)
            check_full_budget(plan,request)
        if package_sha != request['package_sha256']:
            raise ValueError('Resume requires identical package code')
        for u in request['units']:
            if any(fingerprint(v['path']) != v for v in u['inputs'].values()):
                raise ValueError('Resume input fingerprint differs')
    else:
        output = args.output_root.resolve() / ('wgbs_cpg_assignment_profile_v0.3.0_' + stamp())
        output.mkdir(parents=True, exist_ok=False)
        print('OUTPUT_DIRECTORY=' + str(output), flush=True)
        data = select(args)
        resource = concurrency_plan(args.jobs, args.max_memory_gib, args.samtools_threads, host_resources())
        resource.update(input_bam_bytes_recorded=sum(int(r['expected_bam_bytes']) for r in data),
                        max_hours=args.max_hours, max_memory_gib_per_worker=args.max_memory_gib,
                        max_output_gib=args.max_output_gib, raw_inputs_modified=False,
                        initial_runtime_expectation='Pilot: provisionally 5 to 30 minutes including preparation, not a server benchmark. Full mode requires measured pilot planning.')
        write_json(output / 'resource_plan.json', resource)
        print(f"JOBS_REQUESTED={args.jobs} JOBS_PLANNED={resource['effective_jobs']}", flush=True)
        preflight_started = time.monotonic()
        try:
            pf = preflight(data, args, output)
        except BaseException as exc:
            write_json(output/'run_status.json', {'status':'FAILED_PREFLIGHT', 'error':str(exc)})
            archive=Path(str(output)+'_receipt_'+stamp()+'.tar.gz')
            with tarfile.open(archive,'w:gz') as tar:
                for p in sorted(output.rglob('*')):
                    if p.is_file():
                        tar.add(p,arcname=output.name+'/'+str(p.relative_to(output)),recursive=False)
            print('SEND_THIS_FILE='+str(archive),flush=True)
            raise
        device_count = len({u['bam_device'] for u in pf['units']})
        resource['effective_jobs'] = min(resource['effective_jobs'], device_count * args.max_readers_per_device)
        resource['max_readers_per_device'] = args.max_readers_per_device
        resource['input_filesystem_devices'] = device_count
        resource['temporary_directory'] = str(scratch)
        resource['maximum_active_compute_threads'] = resource['effective_jobs'] * (2 + args.samtools_threads)
        resource['max_worker_memory_gib_total'] = resource['effective_jobs'] * args.max_memory_gib
        write_json(output / 'resource_plan.json', resource)
        request = {'mode': args.mode, 'units': pf['units'], 'tools': pf['tools'], 'versions': pf['versions'],
            'lengths': pf['lengths'], 'contract': CONTRACT, 'package_sha256': package_sha,
            'samtools_threads': args.samtools_threads, 'jobs': resource['effective_jobs'], 'max_seconds': args.max_hours*3600,
            'max_memory_gib': args.max_memory_gib, 'max_readers_per_device': args.max_readers_per_device,
            'min_available_memory_bytes': int(args.min_available_memory_gib*1024**3),
            'min_free_disk_bytes': int(args.min_free_disk_gib*1024**3),
            'temporary_root': str(args.tmp_root.resolve()), 'max_output_bytes': int(args.max_output_gib*1024**3)}
        write_json(output / 'run_request.json', request)
        write_json(output / 'preflight_timing.json', {'seconds': time.monotonic()-preflight_started})
    write_json(output / ('scratch_' + stamp() + '.json'), {'directory':str(scratch)})
    try:
        resource_guard(output, request, time.monotonic())
        if request['mode'] == 'preflight':
            print('PASS_PREFLIGHT=' + str(output), flush=True)
            return
        if request['mode'] in ('plan','full') and not args.resume:
            if not args.pilot_result:
                raise ValueError('Provide --pilot-result from the completed identical-input pilot')
            plan = forecast(args.pilot_result, request)
            write_json(output / 'measured_full_forecast.json', plan)
            if request['mode'] == 'plan':
                return
            check_full_budget(plan,request)
        resource = concurrency_plan(request['jobs'], request['max_memory_gib'], request['samtools_threads'], host_resources())
        # A resumed run may use fewer workers if the host is busier, never more than recorded.
        execution_request = dict(request, jobs=resource['effective_jobs'])
        write_json(output / ('execution_resources_' + stamp() + '.json'), resource)
        if request['mode']=='full':
            fresh_plan = dict(read_json(output/'measured_full_forecast.json'))
            fresh_plan['forecast_seconds_with_factor_2'] *= fresh_plan['forecast_jobs'] / min(execution_request['jobs'], len(request['units']))
            check_full_budget(fresh_plan, request)
        started = time.monotonic()
        chroms = ['1','10','18'] if request['mode']=='pilot' else list(map(str,range(1,19)))
        intervals = None
        if request['mode']=='pilot':
            intervals = {c:[max(0,(request['lengths'][c]-100000)//2-5000),
                            min(request['lengths'][c],(request['lengths'][c]-100000)//2+105000)] for c in chroms}
        split_by_sample = {}
        for unit in request['units']:
            resource_guard(output, request, started)
            if unit['sample'] not in split_by_sample:
                print('PREPARE_SNPS=' + unit['sample'], flush=True)
                split_by_sample[unit['sample']] = prepare_split(unit, chroms, output, intervals)
        preparation_seconds = time.monotonic()-started
        subset_hashes = {(sample,c):sha(path/('chr'+c+'.tsv'))
                         for sample,path in split_by_sample.items() for c in chroms}
        tasks = []
        for chrom in chroms:
            for unit in request['units']:
                split = split_by_sample[unit['sample']]
                length = request['lengths'][chrom]
                start = (length-100000)//2 if request['mode']=='pilot' else 0
                end = start+100000 if request['mode']=='pilot' else length
                task = dict(unit, chrom=chrom, start0=start, end0=end, chrom_length=length,
                    tools=request['tools'], split_subset=str(split / ('chr'+chrom+'.tsv')),
                    split_subset_sha256=subset_hashes[unit['sample'],chrom], write_detail=request['mode']=='pilot',
                    max_records=1000000 if request['mode']=='pilot' else 200000000,
                    max_seconds=request['max_seconds'], max_memory_gib=request['max_memory_gib'],
                    samtools_threads=request['samtools_threads'], package_sha256=package_sha)
                task['task_id'] = signature(task)
                tasks.append(task)
        completed, execution_timing = schedule(tasks, execution_request, output, complete_region, started)
        execution_timing['snp_preparation_seconds'] = preparation_seconds
        write_json(output/('execution_timing_' + stamp() + '.json'), execution_timing)
        aggregate(output, request, completed)
        write_json(output/'progress.json', {'status':'COMPLETE', 'completed_regions':len(completed)})
        if request['mode']=='pilot':
            write_json(output/'measured_full_forecast.json', forecast(output, request))
        if request['mode']=='pilot':
            plan=read_json(output/'measured_full_forecast.json')
            print(f"FULL_FORECAST_HOURS={plan['forecast_seconds_with_factor_2']/3600:.2f} FORECAST_JOBS={plan['forecast_jobs']} FORECAST_ONLY=YES", flush=True)
        print('COMPLETE=' + str(output), flush=True)
    except BaseException as exc:
        write_json(output/'run_status.json', {'status':'INCOMPLETE_OR_RESOURCE_STOPPED', 'error':str(exc),
                                            'completed_checkpoints_preserved':True})
        raise
    finally:
        if request['mode'] == 'pilot':
            archive = Path(str(output) + '_receipt_' + stamp() + '.tar.gz')
            with tarfile.open(archive, 'w:gz') as tar:
                for p in sorted(output.rglob('*')):
                    if p.is_file():
                        tar.add(p, arcname=output.name+'/'+str(p.relative_to(output)), recursive=False)
            print('SEND_THIS_FILE=' + str(archive), flush=True)

if __name__ == '__main__':
    main()
