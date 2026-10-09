#!/usr/bin/env python3
"""BAM, VCF, and fragment helpers for exact CpG-K counting."""

from __future__ import annotations

import bisect
import collections
import heapq
import importlib.util
import subprocess
import sys
from pathlib import Path
from typing import Iterable, Iterator, MutableMapping, Sequence

from cpgk_common import CpGKError


def load_core(path: Path):
    spec = importlib.util.spec_from_file_location("reference_fixed_cpgk_core", path)
    if spec is None or spec.loader is None:
        raise CpGKError(f"could not import CpG core: {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)

    def managed_stream(command: Sequence[str]) -> Iterator[str]:
        process = subprocess.Popen(
            list(command), text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
            bufsize=1024 * 1024,
        )
        assert process.stdout is not None
        try:
            yield from process.stdout
        finally:
            process.stdout.close()
            stderr = ""
            if process.stderr is not None:
                stderr = process.stderr.read()
                process.stderr.close()
            returncode = process.wait()
            if returncode != 0:
                raise CpGKError(
                    f"command failed ({returncode}): {' '.join(command)}\n{stderr.strip()}"
                )

    module.stream_command = managed_stream
    return module


def called_nonreference_allele_indices(gt_text: str) -> list[int]:
    separator = "|" if "|" in gt_text else "/" if "/" in gt_text else None
    if separator is None:
        return []
    fields = gt_text.split(separator)
    if len(fields) != 2 or "." in fields:
        return []
    try:
        return sorted({int(value) for value in fields if int(value) > 0})
    except ValueError:
        return []


def disrupted_reference_cpgs(
    *, core, bcftools: str, genotype_vcf: str, sample: str, region,
    cpg_positions0: Sequence[int],
) -> tuple[bytearray, dict[str, int]]:
    samples = core.read_vcf_samples(genotype_vcf)
    if sample not in samples:
        raise CpGKError(f"sample {sample!r} is absent from {genotype_vcf}")
    sample_column = 9 + samples.index(sample)
    disrupted = bytearray(len(cpg_positions0))
    qc: MutableMapping[str, int] = collections.Counter()
    for fields in core.iter_vcf_records(bcftools, genotype_vcf, region):
        qc["genotype_vcf_records_seen"] += 1
        if fields[0] != region.contig or fields[6] not in {"PASS", "."}:
            continue
        format_values = dict(zip(fields[8].split(":"), fields[sample_column].split(":")))
        indices = called_nonreference_allele_indices(format_values.get("GT", "."))
        if not indices:
            continue
        alternate = fields[4].split(",")
        try:
            called = [alternate[index - 1] for index in indices]
        except IndexError:
            qc["genotype_allele_index_out_of_range"] += 1
            continue
        pos0 = int(fields[1]) - 1
        ref = fields[3]
        event_start = pos0
        event_end = pos0 + max(1, len(ref)) - 1
        if any(len(allele) != len(ref) for allele in called if allele != "*"):
            event_start = max(region.start0, event_start - 1)
            event_end = min(region.end0 - 1, event_end + 1)
        left = bisect.bisect_left(cpg_positions0, event_start - 1)
        right = bisect.bisect_right(cpg_positions0, event_end)
        for index in range(left, right):
            cpos0 = cpg_positions0[index]
            if cpos0 <= event_end and cpos0 + 1 >= event_start:
                disrupted[index] = 1
    qc["reference_cpgs_disrupted"] = sum(disrupted)
    return disrupted, dict(qc)


def sam_record_lines(
    *, core, samtools: str, bam: str, region, threads: int,
) -> Iterator[str]:
    command = [samtools, "view", "-@", str(threads), bam, region.samtools]
    for line in core.stream_command(command):
        if line and not line.startswith("@"):
            yield line


def stream_record_groups(
    *, core, records: Iterable[object], min_mapq: int, pair_wait_bp: int,
    qc: MutableMapping[str, int],
) -> Iterator[tuple[object, ...]]:
    active: dict[str, tuple[object, int, int]] = {}
    deadlines: list[tuple[int, int, str]] = []
    serial = 0

    def finalize(qname: str) -> tuple[object, ...] | None:
        item = active.pop(qname, None)
        if item is None:
            return None
        first, _deadline, _serial = item
        qc["fragments_emitted"] += 1
        qc["singleton_fragments"] += 1
        return (first,)

    for record in records:
        qc["sam_records_seen"] += 1
        if not core.record_is_usable(record, min_mapq):
            qc["sam_records_filtered"] += 1
            continue
        while deadlines and deadlines[0][0] < record.pos0:
            deadline, old_serial, qname = heapq.heappop(deadlines)
            current = active.get(qname)
            if current is None or current[2] != old_serial or current[1] != deadline:
                continue
            group = finalize(qname)
            if group is not None:
                qc["pair_cache_expired"] += 1
                yield group
        if record.qname in active:
            first, _deadline, _old_serial = active.pop(record.qname)
            qc["fragments_emitted"] += 1
            qc["paired_fragments"] += 1
            yield first, record
            continue
        is_paired = bool(record.flag & 0x1)
        mate_mapped = not bool(record.flag & 0x8)
        same_contig = record.mate_rname in {"=", record.rname}
        if is_paired and mate_mapped and same_contig and record.mate_pos0 >= 0:
            deadline = max(record.reference_end0, record.mate_pos0 + pair_wait_bp)
            serial += 1
            active[record.qname] = (record, deadline, serial)
            heapq.heappush(deadlines, (deadline, serial, record.qname))
        else:
            qc["fragments_emitted"] += 1
            qc["singleton_fragments"] += 1
            yield (record,)
    for qname in list(active):
        group = finalize(qname)
        if group is not None:
            yield group


def merged_calls(core, records: Sequence[object], extractor, positions=()) -> tuple[dict, int]:
    calls: dict = {}
    conflicts: set = set()
    for record in records:
        extracted = extractor(record, positions) if positions else extractor(record)
        for key, value in extracted.items():
            core.merge_call(calls, conflicts, key, value)
    return calls, len(conflicts)
