"""Read CPU and memory limits for the fragment-counting controller."""

import os
import platform
from pathlib import Path


def host_resources():
    """Read current CPU and memory limits."""
    cpus = (
        len(os.sched_getaffinity(0)) if hasattr(os, "sched_getaffinity") else (os.cpu_count() or 1)
    )
    available = None
    limits = {}
    if Path("/proc/meminfo").is_file():
        mem = dict(line.split(":", 1) for line in Path("/proc/meminfo").read_text().splitlines())
        available = int(mem["MemAvailable"].split()[0]) * 1024
        limits["physical_memory_bytes"] = int(mem["MemTotal"].split()[0]) * 1024
    # Common cgroup v2 and v1 locations; affinity remains an independent limit.
    for maximum, current in [
        ("/sys/fs/cgroup/memory.max", "/sys/fs/cgroup/memory.current"),
        (
            "/sys/fs/cgroup/memory/memory.limit_in_bytes",
            "/sys/fs/cgroup/memory/memory.usage_in_bytes",
        ),
    ]:
        try:
            cap = int(Path(maximum).read_text().strip())
            remaining = max(0, cap - int(Path(current).read_text().strip()))
            available = remaining if available is None else min(available, remaining)
            limits["cgroup_memory_limit_bytes"] = cap
        except (OSError, ValueError):
            pass
    try:
        quota, period = Path("/sys/fs/cgroup/cpu.max").read_text().split()
        cpus = min(cpus, max(1, int(quota) // int(period)))
    except (OSError, ValueError):
        try:
            quota = int(Path("/sys/fs/cgroup/cpu/cpu.cfs_quota_us").read_text())
            period = int(Path("/sys/fs/cgroup/cpu/cpu.cfs_period_us").read_text())
            if quota > 0:
                cpus = min(cpus, max(1, quota // period))
        except (OSError, ValueError):
            pass
    return dict(
        hostname=platform.node(),
        platform=platform.platform(),
        cpu_capacity=cpus,
        memory_available_bytes=available,
        load_average=os.getloadavg() if hasattr(os, "getloadavg") else None,
        detected_limits=limits,
        measurement="Current snapshot; CPU capacity is not a count of idle CPUs",
    )
