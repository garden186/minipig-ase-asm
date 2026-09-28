"""Bounded independent region processes with deterministic checkpoint reuse."""
from __future__ import annotations
from collections import Counter
import math
import os
from pathlib import Path
import platform
import shutil
import signal
import subprocess
import sys
import time
from common import PACKAGE, read_json, write_json

GIB = 1024**3

def host_resources():
    """Read current limits; do not infer hardware from historical host names."""
    cpus = len(os.sched_getaffinity(0)) if hasattr(os, 'sched_getaffinity') else (os.cpu_count() or 1)
    available = None
    limits = {}
    if Path('/proc/meminfo').is_file():
        mem = dict(line.split(':', 1) for line in Path('/proc/meminfo').read_text().splitlines())
        available = int(mem['MemAvailable'].split()[0]) * 1024
        limits['physical_memory_bytes'] = int(mem['MemTotal'].split()[0]) * 1024
    # Common cgroup v2 and v1 locations; affinity remains an independent limit.
    for maximum, current in [('/sys/fs/cgroup/memory.max', '/sys/fs/cgroup/memory.current'),
                             ('/sys/fs/cgroup/memory/memory.limit_in_bytes', '/sys/fs/cgroup/memory/memory.usage_in_bytes')]:
        try:
            cap = int(Path(maximum).read_text().strip())
            remaining = max(0, cap - int(Path(current).read_text().strip()))
            available = remaining if available is None else min(available, remaining)
            limits['cgroup_memory_limit_bytes'] = cap
        except (OSError, ValueError):
            pass
    try:
        quota, period = Path('/sys/fs/cgroup/cpu.max').read_text().split()
        cpus = min(cpus, max(1, int(quota) // int(period)))
    except (OSError, ValueError):
        try:
            quota = int(Path('/sys/fs/cgroup/cpu/cpu.cfs_quota_us').read_text())
            period = int(Path('/sys/fs/cgroup/cpu/cpu.cfs_period_us').read_text())
            if quota > 0:
                cpus = min(cpus, max(1, quota // period))
        except (OSError, ValueError):
            pass
    from server_resources import current_cgroup_limits
    extra = current_cgroup_limits()
    if extra['cpu_capacity'] is not None:
        cpus = min(cpus, max(1, int(extra['cpu_capacity'])))
    if extra['memory_remaining_bytes'] is not None:
        available = extra['memory_remaining_bytes'] if available is None else min(available, extra['memory_remaining_bytes'])
    limits['current_cgroup_ancestors'] = extra
    return dict(hostname=platform.node(), platform=platform.platform(), cpu_capacity=cpus,
                memory_available_bytes=available, load_average=os.getloadavg() if hasattr(os, 'getloadavg') else None,
                detected_limits=limits, measurement='Current snapshot; CPU capacity is not a count of idle CPUs')

def concurrency_plan(jobs, memory_gib, samtools_threads, host):
    if not 1 <= jobs <= 16 or memory_gib <= 0:
        raise ValueError('Use 1 to 16 jobs and a positive per-worker memory budget')
    per_job_threads = 2 + samtools_threads
    cpu_reserve = max(1, math.ceil(host['cpu_capacity'] * 0.125))
    cpu_slots = max(1, (host['cpu_capacity'] - cpu_reserve) // per_job_threads)
    memory = host['memory_available_bytes']
    memory_slots = jobs if memory is None else int((memory * 0.75) // ((memory_gib + 0.25) * GIB))
    if memory_slots < 1:
        raise RuntimeError('Available memory cannot accommodate one worker budget and the memory reserve')
    effective = min(jobs, cpu_slots, memory_slots)
    return dict(requested_jobs=jobs, effective_jobs=effective, cpu_capacity=host['cpu_capacity'],
                cpu_slots=cpu_slots, memory_slots=memory_slots if memory is not None else None,
                maximum_active_compute_threads=effective * per_job_threads,
                max_worker_memory_gib_total=effective * memory_gib,
                memory_available_bytes=memory, host=host,
                policy='Reserve 12.5% CPU capacity and 25% available RAM; allow 0.25 GiB extra per job; at most one active job per BAM; spread readers over resolved filesystem devices',
                memory_limit_note='Worker RSS guard plus controller available-memory guard, not an operating-system memory reservation')

def terminate_child(child):
    """Stop only a process tree launched by this controller."""
    if child.poll() is not None:
        return
    try:
        if os.name == 'posix':
            os.killpg(child.pid, signal.SIGTERM)
        else:
            child.terminate()
        child.wait(timeout=10)
    except ProcessLookupError:
        pass
    except subprocess.TimeoutExpired:
        if os.name == 'posix':
            os.killpg(child.pid, signal.SIGKILL)
        else:
            child.kill()
        child.wait()

def resource_guard(output, request, started):
    if time.monotonic() - started >= request['max_seconds']:
        raise TimeoutError('Run time budget reached; completed checkpoints are preserved')
    used = sum(p.stat().st_size for p in output.rglob('*') if p.is_file())
    if used >= request['max_output_bytes'] or shutil.disk_usage(output).free < request.get('min_free_disk_bytes', 2 * GIB):
        raise RuntimeError('Output storage budget or free-space reserve reached')
    memory = host_resources()['memory_available_bytes']
    if memory is not None and memory < request.get('min_available_memory_bytes', 2 * GIB):
        raise RuntimeError('Available host or cgroup memory fell below the declared reserve')

def task_stem(task):
    stem = task['sample'] + '-' + task['tissue'] + '_chr' + task['chrom']
    return stem + ('_' + str(task['start0']) + '_' + str(task['end0']) if task.get('tile') else '')

def schedule(tasks, request, output, complete_region, started=None, command_factory=None,
             guard=resource_guard, poll_seconds=0.25):
    """Execute processes concurrently; results retain the original task ordering."""
    started = time.monotonic() if started is None else started
    command_factory = command_factory or (lambda task, destination: [sys.executable, str(PACKAGE/'worker.py'), str(task), str(destination)])
    parent = output / 'regions'
    parent.mkdir(exist_ok=True)
    pending, completed, active = [], {}, []
    existing, attempts = {}, {}
    for path in parent.iterdir():
        if path.is_dir() and '_attempt_' in path.name and '.context_' not in path.name:
            existing.setdefault(path.name.rsplit('_attempt_', 1)[0], []).append(path)
        elif path.name.endswith('.json') and '_task_' in path.name:
            stem, number = path.stem.rsplit('_task_', 1)
            attempts[stem] = max(attempts.get(stem, 0), int(number))
    for task in tasks:
        stem = task_stem(task)
        reusable = [p for p in sorted(existing.get(stem, [])) if p.is_dir() and '.context_' not in p.name and complete_region(p, task['task_id'])]
        if reusable:
            completed[task['task_id']] = (task, reusable[-1])
            print('REUSE=' + stem, flush=True)
        else:
            pending.append(task)
    seen = [t['task_id'] for t in tasks]
    if len(set(seen)) != len(seen):
        raise ValueError('Duplicate region task identity')
    events = []
    reused = len(completed)
    last_guard, last_progress, peak_active = float('-inf'), float('-inf'), 0
    try:
        while pending or active:
            now = time.monotonic()
            if now - last_guard >= 5:
                guard(output, request, started)
                last_guard = now
            # Check every active process before launching replacements.
            for item in list(active):
                rc = item['child'].poll()
                if rc is None:
                    if request.get('enforce_process_memory'):
                        import psutil
                        try:
                            proc = psutil.Process(item['child'].pid)
                            used = sum(p.memory_info().rss for p in [proc] + proc.children(recursive=True))
                        except (psutil.NoSuchProcess, psutil.AccessDenied):
                            used = 0
                        if used > request['max_memory_gib'] * GIB:
                            raise MemoryError('Worker process tree exceeded its memory budget')
                    continue
                item['log'].close()
                active.remove(item)
                events.append({'event':'EXIT','task_id':item['task']['task_id'],'time_monotonic':time.monotonic()})
                if rc or not complete_region(item['destination'], item['task']['task_id']):
                    raise RuntimeError('Region did not pass; inspect retained log: ' + task_stem(item['task']))
                completed[item['task']['task_id']] = (item['task'], item['destination'])
                print(f"DONE={task_stem(item['task'])} COMPLETED={len(completed)}/{len(tasks)}", flush=True)
            while pending and len(active) < request['jobs']:
                busy = {a['task']['inputs']['bam']['path'] for a in active}
                devices = Counter(a['task'].get('bam_device') for a in active)
                candidates = [i for i, t in enumerate(pending) if t['inputs']['bam']['path'] not in busy
                              and devices[t.get('bam_device')] < request.get('max_readers_per_device', 16)]
                if not candidates:
                    break
                index = min(candidates, key=lambda i: (devices[pending[i].get('bam_device')], i))
                task = pending.pop(index)
                stem = task_stem(task)
                number = attempts.get(stem, 0) + 1
                attempts[stem] = number
                destination = parent / (stem+f'_attempt_{number:03d}')
                task_path = parent / (stem+f'_task_{number:03d}.json')
                write_json(task_path, task)
                log = open(parent/(stem+f'_attempt_{number:03d}.log'), 'w', encoding='utf-8')
                try:
                    child = subprocess.Popen(command_factory(task_path, destination), stdout=log, stderr=subprocess.STDOUT,
                                             start_new_session=(os.name == 'posix'))
                except BaseException:
                    log.close()
                    raise
                active.append(dict(task=task, destination=destination, child=child, log=log))
                events.append({'event':'START','task_id':task['task_id'],'device':task.get('bam_device'),
                               'time_monotonic':time.monotonic(),'active_on_device':devices[task.get('bam_device')]+1})
                peak_active = max(peak_active, len(active))
                print(f'SCAN={stem} ACTIVE={len(active)} JOBS={request["jobs"]}', flush=True)
            if now - last_progress >= 10:
                write_json(output/'progress.json', dict(status='RUNNING', completed_regions=len(completed),
                           expected_regions=len(tasks), active_regions=len(active), pending_regions=len(pending),
                           elapsed_seconds=now-started, effective_jobs=request['jobs']))
                last_progress = now
            if pending or active:
                time.sleep(poll_seconds)
    finally:
        for item in active:
            terminate_child(item['child'])
            item['log'].close()
        write_json(output/'scheduler_events.json', {'events':events, 'max_readers_per_device':request.get('max_readers_per_device',16)})
    timing = dict(elapsed_seconds=time.monotonic()-started, effective_jobs=request['jobs'],
                  maximum_observed_concurrent_jobs=peak_active, reused_regions=reused, completed_regions=len(completed))
    return [completed[t['task_id']] for t in tasks], timing
