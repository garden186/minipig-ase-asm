#!/usr/bin/env python3
"""Strict three-CpG fragment bimodality discovery and PS-linked ASM testing.

The implementation intentionally uses only the Python standard library plus
samtools/bcftools.  It consumes coordinate-sorted Bismark BAMs as SAM streams,
normalizes XM CpG calls to the reference-C coordinate, merges paired mates by
QNAME, and never treats the two mates as independent observations.

This version intentionally supports exactly three consecutive CpGs. A fully
observed fragment-window is U only for 000, M only for 111, and X for each of
the other six patterns. No continuous methylation-fraction threshold is used
to define the three fragment states.
"""

from __future__ import annotations

import argparse
import bisect
import collections
import dataclasses
import datetime as dt
import gzip
import heapq
import json
import math
import multiprocessing as mp
import os
import queue
import re
import shutil
import subprocess
import sys
import tempfile
import traceback
import zlib
from array import array
from typing import Dict, Iterable, Iterator, List, Mapping, MutableMapping, Optional, Sequence, Tuple


VERSION = "0.7.0-cpg3-strict"
DEFAULT_FRAGMENT_CPGS = 3
ALLOWED_FRAGMENT_CPGS = (3,)
CPG3_PATTERN_LABELS = ("000", "001", "010", "011", "100", "101", "110", "111")
REQUIRED_HAPLOTYPE_PURE_FRAGMENTS = 5
ALLOWED_MAX_MIXED_FRACTIONS = (0.10, 0.20, 0.30)
DISCOVERY_WORKER_BATCH_RECORDS = 4096
CONVERSION_AMBIGUOUS = {frozenset(("C", "T")), frozenset(("G", "A"))}
CIGAR_RE = re.compile(r"(\d+)([MIDNSHP=X])")
QUERY_AND_REF = {"M", "=", "X"}
QUERY_ONLY = {"I", "S"}
REF_ONLY = {"D", "N"}
ACTIVE_TMP_ROOT: Optional[str] = None


def log(message: str) -> None:
    stamp = dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds")
    print(f"[{stamp}] {message}", file=sys.stderr, flush=True)


def die(message: str) -> "NoReturn":
    raise SystemExit(f"ERROR: {message}")


@dataclasses.dataclass(frozen=True)
class Region:
    contig: str
    start0: int
    end0: int

    @property
    def samtools(self) -> str:
        if self.start0 == 0:
            return f"{self.contig}:1-{self.end0}"
        return f"{self.contig}:{self.start0 + 1}-{self.end0}"

    @property
    def label(self) -> str:
        return f"{self.contig}_{self.start0 + 1}_{self.end0}"


@dataclasses.dataclass(frozen=True)
class SNPInfo:
    pos0: int
    ps: str
    h1: str
    h2: str


@dataclasses.dataclass
class SamRecord:
    qname: str
    flag: int
    rname: str
    pos0: int
    mapq: int
    cigar: str
    mate_rname: str
    mate_pos0: int
    tlen: int
    seq: str
    tags: Dict[str, str]

    @property
    def reference_end0(self) -> int:
        rpos = self.pos0
        for length, op in parse_cigar(self.cigar):
            if op in QUERY_AND_REF or op in REF_ONLY:
                rpos += length
        return rpos


@dataclasses.dataclass
class Fragment:
    qname: str
    cpg_calls: Dict[int, int]
    snp_calls: Dict[int, str]
    records: int = 0
    cpg_conflicts: int = 0
    snp_conflicts: int = 0


@dataclasses.dataclass
class CandidateWindow:
    index: int
    region_id: str
    cpg_positions0: Tuple[int, ...]
    u: int
    x: int
    m: int
    pattern_counts: Tuple[int, int, int, int, int, int, int, int]


@dataclasses.dataclass
class CandidateRegion:
    region_id: str
    start_index: int
    end_index: int
    window_indices: Tuple[int, ...]
    start0: int
    end0: int


@dataclasses.dataclass(frozen=True)
class FetchBatch:
    region: Region
    candidate_regions: Tuple[CandidateRegion, ...]


def parse_cigar(cigar: str) -> List[Tuple[int, str]]:
    if cigar == "*":
        return []
    parsed = [(int(n), op) for n, op in CIGAR_RE.findall(cigar)]
    if "".join(f"{n}{op}" for n, op in parsed) != cigar:
        raise ValueError(f"Unsupported or malformed CIGAR: {cigar}")
    return parsed


def parse_region(text: str, contig_lengths: Mapping[str, int]) -> Region:
    text = text.replace(",", "")
    if ":" not in text:
        contig = text
        if contig not in contig_lengths:
            die(f"contig {contig!r} is absent from the reference .fai")
        return Region(contig, 0, contig_lengths[contig])
    contig, coords = text.rsplit(":", 1)
    if "-" not in coords:
        die("region must be CONTIG or CONTIG:START-END")
    start_s, end_s = coords.split("-", 1)
    start1, end1 = int(start_s), int(end_s)
    if contig not in contig_lengths:
        die(f"contig {contig!r} is absent from the reference .fai")
    if start1 < 1 or end1 < start1 or end1 > contig_lengths[contig]:
        die(f"invalid region {text!r} for contig length {contig_lengths[contig]}")
    return Region(contig, start1 - 1, end1)


def read_fai(path: str) -> Dict[str, int]:
    result: Dict[str, int] = {}
    with open(path, "rt", encoding="utf-8") as handle:
        for line in handle:
            if not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            result[fields[0]] = int(fields[1])
    return result


def run_checked(command: Sequence[str], *, capture: bool = True) -> str:
    proc = subprocess.run(
        list(command),
        text=True,
        stdout=subprocess.PIPE if capture else None,
        stderr=subprocess.PIPE,
        check=False,
    )
    if proc.returncode != 0:
        stderr = (proc.stderr or "").strip()
        die(f"command failed ({proc.returncode}): {' '.join(command)}\n{stderr}")
    return proc.stdout or ""


def stream_command(command: Sequence[str]) -> Iterator[str]:
    proc = subprocess.Popen(
        list(command),
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        bufsize=1024 * 1024,
    )
    assert proc.stdout is not None
    try:
        for line in proc.stdout:
            yield line
    finally:
        proc.stdout.close()
        stderr = proc.stderr.read() if proc.stderr is not None else ""
        rc = proc.wait()
        if rc != 0:
            die(f"command failed ({rc}): {' '.join(command)}\n{stderr.strip()}")


def parse_sam_line(line: str) -> SamRecord:
    fields = line.rstrip("\n").split("\t")
    if len(fields) < 11:
        raise ValueError("SAM record has fewer than 11 fields")
    tags: Dict[str, str] = {}
    for field in fields[11:]:
        parts = field.split(":", 2)
        if len(parts) == 3 and parts[0] in {"XM", "XR", "XG", "XX"}:
            tags[parts[0]] = parts[2]
    pos1 = int(fields[3])
    mate_pos1 = int(fields[7]) if fields[7].lstrip("-").isdigit() else 0
    return SamRecord(
        qname=fields[0],
        flag=int(fields[1]),
        rname=fields[2],
        pos0=pos1 - 1,
        mapq=int(fields[4]),
        cigar=fields[5],
        mate_rname=fields[6],
        mate_pos0=mate_pos1 - 1 if mate_pos1 > 0 else -1,
        tlen=int(fields[8]),
        seq=fields[9],
        tags=tags,
    )


def record_is_usable(record: SamRecord, min_mapq: int) -> bool:
    excluded_flags = 0x4 | 0x100 | 0x200 | 0x400 | 0x800
    return (
        not (record.flag & excluded_flags)
        and record.rname != "*"
        and record.cigar != "*"
        and record.mapq >= min_mapq
    )


def aligned_blocks(record: SamRecord) -> Iterator[Tuple[int, int, int]]:
    """Yield (query_start, reference_start0, length) for M/= /X blocks."""
    qpos = 0
    rpos = record.pos0
    for length, op in parse_cigar(record.cigar):
        if op in QUERY_AND_REF:
            yield qpos, rpos, length
            qpos += length
            rpos += length
        elif op in QUERY_ONLY:
            qpos += length
        elif op in REF_ONLY:
            rpos += length
        elif op in {"H", "P"}:
            continue
        else:  # pragma: no cover - parse_cigar already constrains operations
            raise ValueError(f"unsupported CIGAR operation {op}")


def extract_cpg_calls(record: SamRecord) -> Dict[int, int]:
    """Return reference-C positions (0-based) to 0/1 methylation calls."""
    xm = record.tags.get("XM")
    xg = record.tags.get("XG")
    if xm is None or xg not in {"CT", "GA"}:
        return {}
    if len(xm) != len(record.seq):
        raise ValueError(
            f"{record.qname}: XM length {len(xm)} != SEQ length {len(record.seq)}"
        )
    calls: Dict[int, int] = {}
    for qstart, rstart, length in aligned_blocks(record):
        for offset, symbol in enumerate(xm[qstart : qstart + length]):
            if symbol not in {"Z", "z"}:
                continue
            aligned_pos0 = rstart + offset
            cpg_c_pos0 = aligned_pos0 if xg == "CT" else aligned_pos0 - 1
            if cpg_c_pos0 >= 0:
                calls[cpg_c_pos0] = 1 if symbol == "Z" else 0
    return calls


def extract_snp_calls(record: SamRecord, sorted_snp_positions0: Sequence[int]) -> Dict[int, str]:
    if not sorted_snp_positions0 or record.seq == "*":
        return {}
    calls: Dict[int, str] = {}
    sequence = record.seq.upper()
    for qstart, rstart, length in aligned_blocks(record):
        rend = rstart + length
        left = bisect.bisect_left(sorted_snp_positions0, rstart)
        right = bisect.bisect_left(sorted_snp_positions0, rend, lo=left)
        for pos0 in sorted_snp_positions0[left:right]:
            qpos = qstart + (pos0 - rstart)
            if 0 <= qpos < len(sequence):
                base = sequence[qpos]
                if base in {"A", "C", "G", "T"}:
                    calls[pos0] = base
    return calls


def merge_call(target: MutableMapping, conflicts: set, key, value) -> None:
    if key in conflicts:
        return
    if key not in target:
        target[key] = value
    elif target[key] != value:
        target.pop(key, None)
        conflicts.add(key)


def fragment_from_record(
    record: SamRecord,
    snp_positions0: Sequence[int],
) -> Tuple[Fragment, set, set]:
    cpg_calls = extract_cpg_calls(record)
    snp_calls = extract_snp_calls(record, snp_positions0)
    return (
        Fragment(record.qname, cpg_calls, snp_calls, records=1),
        set(),
        set(),
    )


def merge_record_into_fragment(
    state: Tuple[Fragment, set, set],
    record: SamRecord,
    snp_positions0: Sequence[int],
) -> Tuple[Fragment, set, set]:
    fragment, cpg_conflicts, snp_conflicts = state
    for key, value in extract_cpg_calls(record).items():
        merge_call(fragment.cpg_calls, cpg_conflicts, key, value)
    for key, value in extract_snp_calls(record, snp_positions0).items():
        merge_call(fragment.snp_calls, snp_conflicts, key, value)
    fragment.records += 1
    fragment.cpg_conflicts = len(cpg_conflicts)
    fragment.snp_conflicts = len(snp_conflicts)
    return fragment, cpg_conflicts, snp_conflicts


def stream_fragments_from_records(
    records: Iterable[SamRecord],
    snp_positions0: Sequence[int],
    *,
    min_mapq: int,
    pair_wait_bp: int,
    qc: MutableMapping[str, int],
) -> Iterator[Fragment]:
    """Pair coordinate-sorted SAM records with a bounded active cache."""
    active: Dict[str, Tuple[Tuple[Fragment, set, set], int, int]] = {}
    deadlines: List[Tuple[int, int, str]] = []
    serial = 0

    def finalize(qname: str) -> Optional[Fragment]:
        item = active.pop(qname, None)
        if item is None:
            return None
        state, _deadline, _serial = item
        fragment = state[0]
        qc["fragments_emitted"] += 1
        if fragment.records == 1:
            qc["singleton_fragments"] += 1
        qc["overlap_cpg_conflicts"] += fragment.cpg_conflicts
        qc["overlap_snp_conflicts"] += fragment.snp_conflicts
        return fragment

    for record in records:
        qc["sam_records_seen"] += 1
        if not record_is_usable(record, min_mapq):
            qc["sam_records_filtered"] += 1
            continue

        while deadlines and deadlines[0][0] < record.pos0:
            deadline, old_serial, qname = heapq.heappop(deadlines)
            current = active.get(qname)
            if current is None or current[2] != old_serial or current[1] != deadline:
                continue
            fragment = finalize(qname)
            if fragment is not None:
                qc["pair_cache_expired"] += 1
                yield fragment

        if record.qname in active:
            state, _deadline, _old_serial = active.pop(record.qname)
            merged = merge_record_into_fragment(state, record, snp_positions0)
            fragment = merged[0]
            qc["fragments_emitted"] += 1
            qc["paired_fragments"] += 1
            qc["overlap_cpg_conflicts"] += fragment.cpg_conflicts
            qc["overlap_snp_conflicts"] += fragment.snp_conflicts
            yield fragment
            continue

        is_paired = bool(record.flag & 0x1)
        mate_mapped = not bool(record.flag & 0x8)
        same_contig = record.mate_rname in {"=", record.rname}
        if is_paired and mate_mapped and same_contig and record.mate_pos0 >= 0:
            state = fragment_from_record(record, snp_positions0)
            deadline = max(record.reference_end0, record.mate_pos0 + pair_wait_bp)
            serial += 1
            active[record.qname] = (state, deadline, serial)
            heapq.heappush(deadlines, (deadline, serial, record.qname))
        else:
            fragment = fragment_from_record(record, snp_positions0)[0]
            qc["fragments_emitted"] += 1
            qc["singleton_fragments"] += 1
            yield fragment

    for qname in list(active):
        fragment = finalize(qname)
        if fragment is not None:
            yield fragment


def sam_records(
    samtools: str,
    bam: str,
    region: Region,
    threads: int,
) -> Iterator[SamRecord]:
    for line in sam_record_lines(samtools, bam, region, threads):
        yield parse_sam_line(line)


def sam_record_lines(
    samtools: str,
    bam: str,
    region: Region,
    threads: int,
) -> Iterator[str]:
    command = [samtools, "view", "-@", str(threads), bam, region.samtools]
    for line in stream_command(command):
        if line and not line.startswith("@"):
            yield line


def read_vcf_samples(vcf_path: str) -> List[str]:
    opener = gzip.open if vcf_path.endswith(".gz") else open
    with opener(vcf_path, "rt", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("#CHROM"):
                return line.rstrip("\n").split("\t")[9:]
            if not line.startswith("#"):
                break
    die(f"could not find #CHROM header in {vcf_path}")


def load_snpsplit_alleles(path: str, region: Region) -> Dict[int, Tuple[str, str]]:
    result: Dict[int, Tuple[str, str]] = {}
    with open(path, "rt", encoding="utf-8") as handle:
        for line_no, line in enumerate(handle, 1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 5:
                die(f"{path}:{line_no}: expected at least 5 tab-separated fields")
            chrom, pos_s, allele_field = fields[1], fields[2], fields[-1]
            if chrom != region.contig:
                continue
            pos0 = int(pos_s) - 1
            if pos0 < region.start0 or pos0 >= region.end0:
                continue
            allele_parts = allele_field.upper().split("/")
            if len(allele_parts) != 2 or any(len(a) != 1 for a in allele_parts):
                continue
            result[pos0] = (allele_parts[0], allele_parts[1])
    return result


def iter_vcf_records(
    bcftools: str,
    vcf_path: str,
    region: Region,
) -> Iterator[List[str]]:
    command = [bcftools, "view", "-H", "-r", region.samtools, vcf_path]
    for line in stream_command(command):
        fields = line.rstrip("\n").split("\t")
        if len(fields) >= 10:
            yield fields


def load_variant_context(
    *,
    bcftools: str,
    vcf_path: str,
    sample: str,
    snpsplit_path: str,
    extended_region: Region,
    include_conversion_snps: bool,
) -> Tuple[array, Dict[int, SNPInfo], Dict[str, int]]:
    samples = read_vcf_samples(vcf_path)
    if sample not in samples:
        die(f"VCF sample {sample!r} not found; available: {', '.join(samples)}")
    sample_col = 9 + samples.index(sample)
    split_alleles = load_snpsplit_alleles(snpsplit_path, extended_region)
    excluded: List[int] = []
    link_snps: Dict[int, SNPInfo] = {}
    qc = collections.Counter()

    for fields in iter_vcf_records(bcftools, vcf_path, extended_region):
        qc["vcf_records_seen"] += 1
        chrom, pos_s, _id, ref, alt = fields[:5]
        if chrom != extended_region.contig or len(ref) != 1:
            continue
        alleles = [ref.upper()] + [a.upper() for a in alt.split(",")]
        if any(len(a) != 1 for a in alleles):
            continue
        fmt_keys = fields[8].split(":")
        sample_values = fields[sample_col].split(":")
        fmt = dict(zip(fmt_keys, sample_values))
        gt_text = fmt.get("GT", ".")
        separator = "|" if "|" in gt_text else "/" if "/" in gt_text else None
        if separator is None:
            continue
        gt_parts = gt_text.split(separator)
        if len(gt_parts) != 2 or "." in gt_parts:
            continue
        try:
            gt = (int(gt_parts[0]), int(gt_parts[1]))
            h1, h2 = alleles[gt[0]], alleles[gt[1]]
        except (ValueError, IndexError):
            continue
        if h1 == h2:
            continue
        pos0 = int(pos_s) - 1
        if not excluded or excluded[-1] != pos0:
            excluded.append(pos0)
        qc["heterozygous_snv"] += 1

        if separator != "|":
            qc["unphased_het_snv"] += 1
            continue
        ps = fmt.get("PS", ".")
        if ps in {None, "", "."}:
            qc["phased_without_ps"] += 1
            continue
        if pos0 not in split_alleles:
            qc["phased_not_in_snpsplit"] += 1
            continue
        g1, g2 = split_alleles[pos0]
        if (h1, h2) != (g1, g2):
            qc["vcf_snpsplit_allele_mismatch"] += 1
            continue
        if not include_conversion_snps and frozenset((h1, h2)) in CONVERSION_AMBIGUOUS:
            qc["conversion_ambiguous_excluded"] += 1
            continue
        link_snps[pos0] = SNPInfo(pos0, str(ps), h1, h2)
        qc["primary_safe_link_snp"] += 1

    return array("I", excluded), link_snps, dict(qc)


def sorted_contains(values: Sequence[int], target: int) -> bool:
    index = bisect.bisect_left(values, target)
    return index < len(values) and values[index] == target


def enumerate_shared_reference_cpgs(
    *,
    samtools: str,
    fasta: str,
    region: Region,
    excluded_het_positions0: Sequence[int],
) -> Tuple[array, Dict[str, int]]:
    command = [samtools, "faidx", fasta, region.samtools]
    proc = subprocess.Popen(
        command,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        bufsize=1024 * 1024,
    )
    assert proc.stdout is not None
    positions = array("I")
    previous = ""
    pos0 = region.start0
    raw_cpg = 0
    excluded_cpg = 0
    for line in proc.stdout:
        if line.startswith(">"):
            continue
        for base in line.strip().upper():
            if previous == "C" and base == "G":
                cpos0 = pos0 - 1
                if cpos0 >= region.start0 and pos0 < region.end0:
                    raw_cpg += 1
                    if sorted_contains(excluded_het_positions0, cpos0) or sorted_contains(
                        excluded_het_positions0, cpos0 + 1
                    ):
                        excluded_cpg += 1
                    else:
                        positions.append(cpos0)
            previous = base
            pos0 += 1
    proc.stdout.close()
    stderr = proc.stderr.read() if proc.stderr is not None else ""
    rc = proc.wait()
    if rc != 0:
        die(f"samtools faidx failed ({rc}): {stderr.strip()}")
    if pos0 != region.end0:
        die(
            f"reference length mismatch while reading {region.samtools}: "
            f"observed {pos0 - region.start0}, expected {region.end0 - region.start0}"
        )
    qc = {
        "reference_cpg_raw": raw_cpg,
        "reference_cpg_excluded_het_c_or_g": excluded_cpg,
        "reference_cpg_shared": len(positions),
    }
    return positions, qc


def position_to_cpg_index(cpg_positions0: Sequence[int], pos0: int) -> Optional[int]:
    index = bisect.bisect_left(cpg_positions0, pos0)
    if index < len(cpg_positions0) and cpg_positions0[index] == pos0:
        return index
    return None


def complete_window_patterns(
    fragment_calls: Mapping[int, int],
    cpg_positions0: Sequence[int],
    fragment_cpgs: int = DEFAULT_FRAGMENT_CPGS,
) -> Iterator[Tuple[int, str, int, int]]:
    if fragment_cpgs != 3:
        die("strict CpG3 analysis requires fragment_cpgs=3")
    indexed: List[Tuple[int, int]] = []
    for pos0, call in fragment_calls.items():
        index = position_to_cpg_index(cpg_positions0, pos0)
        if index is not None:
            if call not in (0, 1):
                die(f"invalid CpG methylation call {call!r}; expected 0 or 1")
            indexed.append((index, call))
    indexed.sort()
    run_start = 0
    while run_start < len(indexed):
        run_end = run_start + 1
        while run_end < len(indexed) and indexed[run_end][0] == indexed[run_end - 1][0] + 1:
            run_end += 1
        run = indexed[run_start:run_end]
        if len(run) >= 3:
            for offset in range(0, len(run) - 2):
                window = run[offset : offset + 3]
                calls = tuple(call for _index, call in window)
                pattern_code = (calls[0] << 2) | (calls[1] << 1) | calls[2]
                methylated = calls[0] + calls[1] + calls[2]
                state = "U" if pattern_code == 0 else "M" if pattern_code == 7 else "X"
                yield window[0][0], state, methylated, pattern_code
        run_start = run_end


def count_discovery_fragments(
    fragments: Iterable[Fragment],
    *,
    cpg_positions0: Sequence[int],
    fragment_cpgs: int,
    qc: MutableMapping[str, int],
) -> Tuple[array, array, array, array, array, array, array, array]:
    if fragment_cpgs != 3:
        die("strict CpG3 analysis requires fragment_cpgs=3")
    window_count = max(0, len(cpg_positions0) - 2)
    pattern_counts = tuple(array("I", [0]) * window_count for _ in range(8))
    for fragment in fragments:
        any_window = False
        for window_index, _state, _methylated, pattern_code in complete_window_patterns(
            fragment.cpg_calls, cpg_positions0
        ):
            any_window = True
            pattern_counts[pattern_code][window_index] += 1
        if any_window:
            qc["fragments_with_complete_3cpg_window"] += 1
    return pattern_counts


def discovery_worker(
    worker_id: int,
    input_queue,
    result_queue,
    cpg_positions0: Sequence[int],
    fragment_cpgs: int,
    min_mapq: int,
    pair_wait_bp: int,
) -> None:
    """Consume one deterministic QNAME shard and return compressed count arrays."""
    try:
        saw_end_marker = False

        def records() -> Iterator[SamRecord]:
            nonlocal saw_end_marker
            while True:
                item = input_queue.get()
                if isinstance(item, tuple) and len(item) == 2 and item[0] == "END":
                    saw_end_marker = True
                    # Advance every shard to the final coordinate observed in the
                    # original stream so pair-cache expiration matches one-worker mode.
                    yield SamRecord(
                        qname=f"__discovery_worker_end_{worker_id}",
                        flag=0x4,
                        rname="*",
                        pos0=int(item[1]),
                        mapq=0,
                        cigar="*",
                        mate_rname="*",
                        mate_pos0=-1,
                        tlen=0,
                        seq="*",
                        tags={},
                    )
                    return
                for line in item:
                    yield parse_sam_line(line)

        qc: MutableMapping[str, int] = collections.Counter()
        fragments = stream_fragments_from_records(
            records(),
            (),
            min_mapq=min_mapq,
            pair_wait_bp=pair_wait_bp,
            qc=qc,
        )
        pattern_counts = count_discovery_fragments(
            fragments,
            cpg_positions0=cpg_positions0,
            fragment_cpgs=fragment_cpgs,
            qc=qc,
        )
        if saw_end_marker:
            qc["sam_records_seen"] -= 1
            qc["sam_records_filtered"] -= 1
        result_queue.put(
            (
                "OK",
                worker_id,
                tuple(zlib.compress(counts.tobytes(), level=1) for counts in pattern_counts),
                dict(qc),
            )
        )
    except Exception:
        result_queue.put(("ERROR", worker_id, traceback.format_exc()))


def put_worker_item(input_queue, item, process: mp.Process) -> None:
    while True:
        try:
            input_queue.put(item, timeout=1.0)
            return
        except queue.Full:
            if not process.is_alive():
                raise RuntimeError(
                    f"discovery worker {process.name} exited before consuming its input"
                )


def decode_count_array(payload: bytes, expected_length: int) -> array:
    values = array("I")
    values.frombytes(zlib.decompress(payload))
    if len(values) != expected_length:
        raise RuntimeError(
            f"parallel discovery returned {len(values)} counts; expected {expected_length}"
        )
    return values


def parallel_discovery_pass(
    *,
    samtools: str,
    unsplit_bam: str,
    region: Region,
    cpg_positions0: Sequence[int],
    fragment_cpgs: int,
    threads: int,
    workers: int,
    min_mapq: int,
    pair_wait_bp: int,
) -> Tuple[Tuple[array, ...], Dict[str, int]]:
    if "fork" not in mp.get_all_start_methods():
        die("--workers > 1 requires a platform with multiprocessing fork support")
    context = mp.get_context("fork")
    input_queues = [context.Queue(maxsize=4) for _ in range(workers)]
    result_queue = context.Queue()
    processes = [
        context.Process(
            target=discovery_worker,
            name=f"discovery-{worker_id}",
            args=(
                worker_id,
                input_queues[worker_id],
                result_queue,
                cpg_positions0,
                fragment_cpgs,
                min_mapq,
                pair_wait_bp,
            ),
        )
        for worker_id in range(workers)
    ]
    for process in processes:
        process.start()

    batches: List[List[str]] = [[] for _ in range(workers)]
    last_pos0 = region.start0
    records_dispatched = 0
    try:
        for line in sam_record_lines(samtools, unsplit_bam, region, threads):
            fields = line.split("\t", 4)
            if len(fields) < 4:
                raise ValueError("SAM record has fewer than four fields")
            worker_id = zlib.crc32(fields[0].encode("utf-8")) % workers
            last_pos0 = int(fields[3]) - 1
            batches[worker_id].append(line)
            records_dispatched += 1
            if len(batches[worker_id]) >= DISCOVERY_WORKER_BATCH_RECORDS:
                put_worker_item(
                    input_queues[worker_id], batches[worker_id], processes[worker_id]
                )
                batches[worker_id] = []
            if records_dispatched % 5_000_000 == 0:
                log(f"parallel discovery dispatch: SAM records={records_dispatched:,}")
        for worker_id in range(workers):
            if batches[worker_id]:
                put_worker_item(
                    input_queues[worker_id], batches[worker_id], processes[worker_id]
                )
            put_worker_item(
                input_queues[worker_id], ("END", last_pos0), processes[worker_id]
            )
    except Exception:
        for process in processes:
            if process.is_alive():
                process.terminate()
        for process in processes:
            process.join()
        raise

    worker_results = []
    while len(worker_results) < workers:
        try:
            result = result_queue.get(timeout=1.0)
        except queue.Empty:
            failed = [
                process
                for process in processes
                if process.exitcode is not None and process.exitcode != 0
            ]
            if failed:
                die(
                    "parallel discovery worker failed: "
                    + ", ".join(f"{p.name} exit={p.exitcode}" for p in failed)
                )
            continue
        if result[0] == "ERROR":
            for process in processes:
                if process.is_alive():
                    process.terminate()
            for process in processes:
                process.join()
            die(f"parallel discovery worker {result[1]} failed:\n{result[2]}")
        worker_results.append(result)

    for process in processes:
        process.join()
        if process.exitcode != 0:
            die(f"parallel discovery worker {process.name} exited with {process.exitcode}")

    window_count = max(0, len(cpg_positions0) - 2)
    pattern_counts = tuple(array("I", [0]) * window_count for _ in range(8))
    qc: MutableMapping[str, int] = collections.Counter()
    for _status, _worker_id, pattern_payloads, local_qc in worker_results:
        if len(pattern_payloads) != 8:
            die(f"parallel discovery returned {len(pattern_payloads)} pattern arrays; expected 8")
        local_arrays = tuple(
            decode_count_array(payload, window_count) for payload in pattern_payloads
        )
        for destination, source in zip(pattern_counts, local_arrays):
            for index, value in enumerate(source):
                destination[index] += value
        qc.update(local_qc)
    qc["python_worker_processes"] = workers
    return pattern_counts, dict(qc)


def discovery_pass(
    *,
    samtools: str,
    unsplit_bam: str,
    region: Region,
    cpg_positions0: Sequence[int],
    fragment_cpgs: int,
    threads: int,
    workers: int,
    min_mapq: int,
    pair_wait_bp: int,
) -> Tuple[Tuple[array, ...], Dict[str, int]]:
    if workers > 1:
        return parallel_discovery_pass(
            samtools=samtools,
            unsplit_bam=unsplit_bam,
            region=region,
            cpg_positions0=cpg_positions0,
            fragment_cpgs=fragment_cpgs,
            threads=threads,
            workers=workers,
            min_mapq=min_mapq,
            pair_wait_bp=pair_wait_bp,
        )
    qc: MutableMapping[str, int] = collections.Counter()
    records = sam_records(samtools, unsplit_bam, region, threads)
    fragments = stream_fragments_from_records(
        records,
        (),
        min_mapq=min_mapq,
        pair_wait_bp=pair_wait_bp,
        qc=qc,
    )
    pattern_counts = count_discovery_fragments(
        fragments,
        cpg_positions0=cpg_positions0,
        fragment_cpgs=fragment_cpgs,
        qc=qc,
    )
    qc["python_worker_processes"] = 1
    return pattern_counts, dict(qc)


def is_bimodal_candidate(
    u: int,
    x: int,
    m: int,
    *,
    min_fragments: int,
    min_mode_count: int,
    min_minor_mode_fraction: float,
    max_mixed_fraction: float,
) -> bool:
    total = u + x + m
    pure = u + m
    if total < min_fragments or u < min_mode_count or m < min_mode_count or pure == 0:
        return False
    return (
        min(u, m) / pure >= min_minor_mode_fraction
        and x / total <= max_mixed_fraction
    )


def build_candidates(
    cpg_positions0: Sequence[int],
    pattern_counts: Sequence[Sequence[int]],
    args: argparse.Namespace,
) -> Tuple[List[CandidateWindow], List[CandidateRegion]]:
    if len(pattern_counts) != 8:
        die(f"strict CpG3 discovery requires 8 pattern arrays; observed {len(pattern_counts)}")
    window_count = max(0, len(cpg_positions0) - 2)
    if any(len(counts) != window_count for counts in pattern_counts):
        die("CpG3 pattern count arrays have inconsistent lengths")

    def counts_for(index: int) -> Tuple[Tuple[int, ...], int, int, int]:
        counts = tuple(int(pattern_counts[code][index]) for code in range(8))
        return counts, counts[0], sum(counts[1:7]), counts[7]

    candidate_indices = []
    for index in range(window_count):
        _counts, u, x, m = counts_for(index)
        if is_bimodal_candidate(
            u,
            x,
            m,
            min_fragments=args.min_fragments,
            min_mode_count=args.min_mode_count,
            min_minor_mode_fraction=args.min_minor_mode_fraction,
            max_mixed_fraction=args.max_mixed_fraction,
        ):
            candidate_indices.append(index)

    region_groups: List[List[int]] = []
    for index in candidate_indices:
        if not region_groups or index > region_groups[-1][-1] + 2:
            region_groups.append([index])
        else:
            region_groups[-1].append(index)

    regions: List[CandidateRegion] = []
    index_to_region: Dict[int, str] = {}
    for number, group in enumerate(region_groups, 1):
        region_id = f"BR{number:07d}"
        start_index = group[0]
        end_index = group[-1] + 2
        region = CandidateRegion(
            region_id=region_id,
            start_index=start_index,
            end_index=end_index,
            window_indices=tuple(group),
            start0=cpg_positions0[start_index],
            end0=cpg_positions0[end_index] + 2,
        )
        regions.append(region)
        for index in group:
            index_to_region[index] = region_id

    windows = []
    for index in candidate_indices:
        counts, u, x, m = counts_for(index)
        windows.append(
            CandidateWindow(
                index=index,
                region_id=index_to_region[index],
                cpg_positions0=tuple(cpg_positions0[index : index + 3]),
                u=u,
                x=x,
                m=m,
                pattern_counts=counts,
            )
        )
    return windows, regions


def infer_ps_states(
    fragment: Fragment,
    snp_by_pos0: Mapping[int, SNPInfo],
    min_link_snps: int,
) -> Tuple[Dict[str, str], Dict[str, Tuple[int, int]]]:
    votes: Dict[str, List[int]] = collections.defaultdict(lambda: [0, 0])
    for pos0, base in fragment.snp_calls.items():
        info = snp_by_pos0.get(pos0)
        if info is None:
            continue
        if base == info.h1:
            votes[info.ps][0] += 1
        elif base == info.h2:
            votes[info.ps][1] += 1
    states: Dict[str, str] = {}
    vote_summary: Dict[str, Tuple[int, int]] = {}
    for ps, (h1_votes, h2_votes) in votes.items():
        vote_summary[ps] = (h1_votes, h2_votes)
        if h1_votes >= min_link_snps and h2_votes == 0:
            states[ps] = "H1"
        elif h2_votes >= min_link_snps and h1_votes == 0:
            states[ps] = "H2"
    return states, vote_summary


def fisher_exact_two_sided(a: int, b: int, c: int, d: int) -> float:
    """Two-sided Fisher exact p-value for [[a,b],[c,d]]."""
    row1, row2 = a + b, c + d
    col1 = a + c
    total = row1 + row2
    if total == 0:
        return math.nan

    def log_choose(n: int, k: int) -> float:
        if k < 0 or k > n:
            return -math.inf
        return math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1)

    def probability(x: int) -> float:
        return math.exp(log_choose(col1, x) + log_choose(total - col1, row1 - x) - log_choose(total, row1))

    lo = max(0, row1 - (total - col1))
    hi = min(row1, col1)
    observed = probability(a)
    return min(1.0, sum(probability(x) for x in range(lo, hi + 1) if probability(x) <= observed * (1 + 1e-12)))


def bh_adjust(pvalues: Sequence[float]) -> List[float]:
    result = [math.nan] * len(pvalues)
    valid = [(p, index) for index, p in enumerate(pvalues) if not math.isnan(p)]
    valid.sort(reverse=True)
    running = 1.0
    m = len(valid)
    for reverse_rank, (pvalue, index) in enumerate(valid, 1):
        rank = m - reverse_rank + 1
        running = min(running, pvalue * m / rank)
        result[index] = min(1.0, running)
    return result


def new_haplotype_counter() -> List[int]:
    # H1_U,H1_X,H1_M,H2_U,H2_X,H2_M,H1_meth,H1_calls,H2_meth,H2_calls
    return [0] * 10


def make_fetch_batches(
    candidate_regions: Sequence[CandidateRegion],
    fetch_bounds: Region,
    fetch_padding_bp: int,
    fetch_batch_gap_bp: int,
) -> List[FetchBatch]:
    batches: List[FetchBatch] = []
    current_regions: List[CandidateRegion] = []
    current_start = -1
    current_end = -1
    for candidate_region in sorted(candidate_regions, key=lambda item: item.start0):
        start0 = max(fetch_bounds.start0, candidate_region.start0 - fetch_padding_bp)
        end0 = min(fetch_bounds.end0, candidate_region.end0 + fetch_padding_bp)
        if current_regions and start0 <= current_end + fetch_batch_gap_bp:
            current_regions.append(candidate_region)
            current_end = max(current_end, end0)
            continue
        if current_regions:
            batches.append(
                FetchBatch(
                    Region(fetch_bounds.contig, current_start, current_end),
                    tuple(current_regions),
                )
            )
        current_regions = [candidate_region]
        current_start, current_end = start0, end0
    if current_regions:
        batches.append(
            FetchBatch(
                Region(fetch_bounds.contig, current_start, current_end),
                tuple(current_regions),
            )
        )
    return batches


def collect_haplotype_counts_for_bam(
    *,
    source_haplotype: str,
    samtools: str,
    bam: str,
    analysis_region: Region,
    fetch_bounds: Region,
    candidate_regions: Sequence[CandidateRegion],
    cpg_positions0: Sequence[int],
    fragment_cpgs: int,
    snp_by_pos0: Mapping[int, SNPInfo],
    stats: MutableMapping[Tuple[int, str], List[int]],
    threads: int,
    min_mapq: int,
    pair_wait_bp: int,
    fetch_padding_bp: int,
    fetch_batch_gap_bp: int,
    min_link_snps: int,
    qc: MutableMapping[str, int],
) -> None:
    sorted_snp_positions0 = sorted(snp_by_pos0)
    batches = make_fetch_batches(
        candidate_regions, fetch_bounds, fetch_padding_bp, fetch_batch_gap_bp
    )
    qc[f"{source_haplotype}_fetch_batches"] = len(batches)
    progress_every = max(1, len(batches) // 20) if batches else 1
    for batch_number, batch in enumerate(batches, 1):
        region_windows = {
            index
            for candidate_region in batch.candidate_regions
            for index in candidate_region.window_indices
        }
        local_qc: MutableMapping[str, int] = collections.Counter()
        records = sam_records(samtools, bam, batch.region, threads)
        fragments = stream_fragments_from_records(
            records,
            sorted_snp_positions0,
            min_mapq=min_mapq,
            pair_wait_bp=pair_wait_bp,
            qc=local_qc,
        )
        for fragment in fragments:
            patterns = [
                pattern
                for pattern in complete_window_patterns(
                    fragment.cpg_calls, cpg_positions0, fragment_cpgs
                )
                if pattern[0] in region_windows
            ]
            if not patterns:
                continue
            qc[f"{source_haplotype}_fragments_with_candidate_window"] += 1
            states, _votes = infer_ps_states(fragment, snp_by_pos0, min_link_snps)
            supported_ps = [ps for ps, state in states.items() if state == source_haplotype]
            discordant_ps = [ps for ps, state in states.items() if state != source_haplotype]
            qc[f"{source_haplotype}_discordant_ps_links"] += len(discordant_ps)
            if not supported_ps:
                qc[f"{source_haplotype}_candidate_fragments_without_safe_ps_link"] += 1
                continue
            for window_index, pattern, methylated, _pattern_code in patterns:
                for ps in supported_ps:
                    counts = stats.setdefault((window_index, ps), new_haplotype_counter())
                    base_index = 0 if source_haplotype == "H1" else 3
                    pattern_offset = {"U": 0, "X": 1, "M": 2}[pattern]
                    counts[base_index + pattern_offset] += 1
                    if source_haplotype == "H1":
                        counts[6] += methylated
                        counts[7] += fragment_cpgs
                    else:
                        counts[8] += methylated
                        counts[9] += fragment_cpgs
        for key, value in local_qc.items():
            qc[f"{source_haplotype}_{key}"] += value
        if batch_number % progress_every == 0 or batch_number == len(batches):
            log(
                f"{source_haplotype} progress: fetch batches={batch_number:,}/{len(batches):,}; "
                f"window×PS rows={len(stats):,}"
            )



def configure_tmp_root(path: str) -> str:
    global ACTIVE_TMP_ROOT
    resolved = os.path.realpath(os.path.abspath(path))
    if resolved == "/tmp" or resolved.startswith("/tmp/"):
        die("--tmp-root must not be /tmp or a path below /tmp")
    os.makedirs(resolved, exist_ok=True)
    if not os.path.isdir(resolved) or not os.access(resolved, os.W_OK):
        die(f"tmp root is not a writable directory: {resolved}")
    for name in ("TMPDIR", "TMP", "TEMP"):
        os.environ[name] = resolved
    tempfile.tempdir = resolved
    ACTIVE_TMP_ROOT = resolved
    return resolved


def atomic_text_writer(path: str, gzipped: bool = False):
    if ACTIVE_TMP_ROOT is None:
        die("temporary root is not configured")
    parent = os.path.dirname(os.path.abspath(path))
    os.makedirs(parent, exist_ok=True)
    fd, temp_path = tempfile.mkstemp(
        prefix=f".{os.path.basename(path)}.", dir=ACTIVE_TMP_ROOT
    )
    os.close(fd)
    if gzipped:
        handle = gzip.open(temp_path, "wt", encoding="utf-8", newline="")
    else:
        handle = open(temp_path, "wt", encoding="utf-8", newline="")

    class Context:
        def __enter__(self):
            return handle

        def __exit__(self, exc_type, exc, tb):
            handle.close()
            if exc_type is None:
                os.replace(temp_path, path)
            else:
                try:
                    os.unlink(temp_path)
                except FileNotFoundError:
                    pass
            return False

    return Context()


def write_candidates(
    prefix: str,
    sample: str,
    tissue: str,
    contig: str,
    windows: Sequence[CandidateWindow],
    regions: Sequence[CandidateRegion],
) -> Tuple[str, str]:
    windows_path = f"{prefix}.candidate_windows.tsv.gz"
    regions_path = f"{prefix}.candidate_regions.bed"
    pattern_headers = [f"pattern_{label}_fragments" for label in CPG3_PATTERN_LABELS]
    header = [
        "sample",
        "tissue",
        "window_index",
        "region_id",
        "chrom",
        "start0",
        "end0",
        "cpg_positions_1based",
        *pattern_headers,
        "U_fragments",
        "X_fragments",
        "M_fragments",
        "total_fragments",
        "minor_mode_fraction",
        "mixed_fraction",
    ]
    with atomic_text_writer(windows_path, gzipped=True) as out:
        out.write("\t".join(header) + "\n")
        for window in windows:
            u, x, m = window.u, window.x, window.m
            pure = u + m
            total = pure + x
            row = [
                sample,
                tissue,
                str(window.index),
                window.region_id,
                contig,
                str(window.cpg_positions0[0]),
                str(window.cpg_positions0[-1] + 2),
                ",".join(str(p + 1) for p in window.cpg_positions0),
                *[str(value) for value in window.pattern_counts],
                str(u),
                str(x),
                str(m),
                str(total),
                f"{min(u, m) / pure:.8g}" if pure else "NA",
                f"{x / total:.8g}" if total else "NA",
            ]
            out.write("\t".join(row) + "\n")
    with atomic_text_writer(regions_path) as out:
        for region in regions:
            out.write(
                f"{contig}\t{region.start0}\t{region.end0}\t{region.region_id}"
                f"\t{len(region.window_indices)}\n"
            )
    return windows_path, regions_path


def load_candidates(
    path: str,
    cpg_positions0: Sequence[int],
    fragment_cpgs: int,
) -> Tuple[List[CandidateWindow], List[CandidateRegion]]:
    if fragment_cpgs != 3:
        die("strict CpG3 analysis requires fragment_cpgs=3")
    windows: List[CandidateWindow] = []
    opener = gzip.open if path.endswith(".gz") else open
    pattern_columns = [f"pattern_{label}_fragments" for label in CPG3_PATTERN_LABELS]
    with opener(path, "rt", encoding="utf-8") as handle:
        header = handle.readline().rstrip("\n").split("\t")
        columns = {name: index for index, name in enumerate(header)}
        required = {
            "window_index", "region_id", "cpg_positions_1based",
            "U_fragments", "X_fragments", "M_fragments", *pattern_columns,
        }
        missing = required - columns.keys()
        if missing:
            die(f"candidate file is missing columns: {', '.join(sorted(missing))}")
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            index = int(fields[columns["window_index"]])
            positions0 = tuple(
                int(value) - 1
                for value in fields[columns["cpg_positions_1based"]].split(",")
            )
            expected = tuple(cpg_positions0[index : index + 3])
            if positions0 != expected:
                die(
                    f"candidate CpG skeleton mismatch at window_index={index}; "
                    "confirm that discover and haplotype used the same reference and strict CpG3 code"
                )
            pattern_counts = tuple(int(fields[columns[name]]) for name in pattern_columns)
            u, x, m = pattern_counts[0], sum(pattern_counts[1:7]), pattern_counts[7]
            reported = (
                int(fields[columns["U_fragments"]]),
                int(fields[columns["X_fragments"]]),
                int(fields[columns["M_fragments"]]),
            )
            if reported != (u, x, m):
                die(f"candidate pattern/state count mismatch at window_index={index}")
            windows.append(
                CandidateWindow(
                    index=index,
                    region_id=fields[columns["region_id"]],
                    cpg_positions0=positions0,
                    u=u,
                    x=x,
                    m=m,
                    pattern_counts=pattern_counts,
                )
            )
    by_region: Dict[str, List[int]] = collections.defaultdict(list)
    for window in windows:
        by_region[window.region_id].append(window.index)
    regions: List[CandidateRegion] = []
    for region_id, indices in by_region.items():
        indices.sort()
        start_index = indices[0]
        end_index = indices[-1] + 2
        regions.append(
            CandidateRegion(
                region_id,
                start_index,
                end_index,
                tuple(indices),
                cpg_positions0[start_index],
                cpg_positions0[end_index] + 2,
            )
        )
    regions.sort(key=lambda item: item.start0)
    return windows, regions


def format_float(value: float) -> str:
    return "NA" if math.isnan(value) else f"{value:.10g}"


def write_haplotype_stats(
    *,
    prefix: str,
    sample: str,
    tissue: str,
    contig: str,
    windows_by_index: Mapping[int, CandidateWindow],
    stats: Mapping[Tuple[int, str], List[int]],
    min_haplotype_fragments: int,
) -> str:
    if min_haplotype_fragments != REQUIRED_HAPLOTYPE_PURE_FRAGMENTS:
        die(
            "strict CpG3 primary association requires exactly "
            f"{REQUIRED_HAPLOTYPE_PURE_FRAGMENTS} U+M fragments per haplotype"
        )
    path = f"{prefix}.haplotype_uxm.tsv.gz"
    rows = []
    for (window_index, ps), counts in sorted(
        stats.items(),
        key=lambda item: (
            item[0][0],
            int(item[0][1]) if item[0][1].isdigit() else math.inf,
            item[0][1],
        ),
    ):
        h1_u, h1_x, h1_m, h2_u, h2_x, h2_m, h1_meth, h1_calls, h2_meth, h2_calls = counts
        h1_pure, h2_pure = h1_u + h1_m, h2_u + h2_m
        tested = h1_pure >= min_haplotype_fragments and h2_pure >= min_haplotype_fragments
        if tested:
            pvalue = fisher_exact_two_sided(h1_m, h1_u, h2_m, h2_u)
            reason = "TESTED"
        else:
            pvalue = math.nan
            insufficient = []
            if h1_pure < min_haplotype_fragments:
                insufficient.append("H1")
            if h2_pure < min_haplotype_fragments:
                insufficient.append("H2")
            reason = "INSUFFICIENT_" + "_AND_".join(insufficient) + "_PURE_FRAGMENTS"
        h1_mode = h1_m / h1_pure if h1_pure else math.nan
        h2_mode = h2_m / h2_pure if h2_pure else math.nan
        h1_mean = h1_meth / h1_calls if h1_calls else math.nan
        h2_mean = h2_meth / h2_calls if h2_calls else math.nan
        rows.append(
            {
                "window_index": window_index,
                "ps": ps,
                "counts": counts,
                "h1_pure": h1_pure,
                "h2_pure": h2_pure,
                "pvalue": pvalue,
                "reason": reason,
                "h1_mode": h1_mode,
                "h2_mode": h2_mode,
                "h1_mean": h1_mean,
                "h2_mean": h2_mean,
            }
        )
    qvalues = bh_adjust([row["pvalue"] for row in rows])
    header = [
        "sample", "tissue", "window_index", "region_id", "chrom", "start0", "end0",
        "cpg_positions_1based", "PS",
        "H1_U", "H1_X", "H1_M", "H2_U", "H2_X", "H2_M",
        "H1_pure_U_plus_M", "H2_pure_U_plus_M", "required_pure_per_haplotype",
        "H1_M_mode_fraction", "H2_M_mode_fraction", "delta_M_mode_fraction_H1_minus_H2",
        "H1_mean_methylation", "H2_mean_methylation", "delta_mean_methylation_H1_minus_H2",
        "fisher_p", "BH_q_run", "tested", "testability_reason",
    ]
    with atomic_text_writer(path, gzipped=True) as out:
        out.write("\t".join(header) + "\n")
        for row, qvalue in zip(rows, qvalues):
            window = windows_by_index[row["window_index"]]
            counts = row["counts"]
            h1_mode, h2_mode = row["h1_mode"], row["h2_mode"]
            h1_mean, h2_mean = row["h1_mean"], row["h2_mean"]
            values = [
                sample, tissue, str(row["window_index"]), window.region_id, contig,
                str(window.cpg_positions0[0]), str(window.cpg_positions0[-1] + 2),
                ",".join(str(p + 1) for p in window.cpg_positions0), str(row["ps"]),
                *[str(value) for value in counts[:6]],
                str(row["h1_pure"]), str(row["h2_pure"]),
                str(REQUIRED_HAPLOTYPE_PURE_FRAGMENTS),
                format_float(h1_mode), format_float(h2_mode),
                format_float(h1_mode - h2_mode)
                if not math.isnan(h1_mode) and not math.isnan(h2_mode) else "NA",
                format_float(h1_mean), format_float(h2_mean),
                format_float(h1_mean - h2_mean)
                if not math.isnan(h1_mean) and not math.isnan(h2_mean) else "NA",
                format_float(row["pvalue"]), format_float(qvalue),
                "1" if not math.isnan(row["pvalue"]) else "0", row["reason"],
            ]
            out.write("\t".join(values) + "\n")
    return path



def file_fingerprint(path: str) -> Dict[str, object]:
    stat = os.stat(path)
    return {
        "path": os.path.abspath(path),
        "size_bytes": stat.st_size,
        "mtime": dt.datetime.fromtimestamp(stat.st_mtime, dt.timezone.utc).isoformat(),
    }


def tool_version(tool: str) -> str:
    proc = subprocess.run([tool, "--version"], text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    return (proc.stdout or "").splitlines()[0] if proc.returncode == 0 else "unknown"


def bam_contig_length(samtools: str, bam: str, contig: str) -> int:
    output = run_checked([samtools, "idxstats", bam])
    for line in output.splitlines():
        fields = line.split("\t")
        if fields and fields[0] == contig:
            return int(fields[1])
    die(f"contig {contig!r} not found in {bam}")


def bam_sort_order(samtools: str, bam: str) -> str:
    header = run_checked([samtools, "view", "-H", bam])
    for line in header.splitlines():
        if not line.startswith("@HD"):
            continue
        for field in line.split("\t")[1:]:
            if field.startswith("SO:"):
                return field[3:]
    return "unknown"


def sample_bam_tags(samtools: str, bam: str, region: Region, n: int = 1000) -> Dict[str, int]:
    counts = collections.Counter()
    command = [samtools, "view", bam, region.samtools]
    proc = subprocess.Popen(command, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    assert proc.stdout is not None
    for index, line in enumerate(proc.stdout):
        record = parse_sam_line(line)
        counts["records"] += 1
        for tag in ("XM", "XR", "XG", "XX"):
            if tag in record.tags:
                counts[tag] += 1
        if index + 1 >= n:
            proc.terminate()
            break
    proc.stdout.close()
    proc.wait()
    return dict(counts)


def preflight(args: argparse.Namespace, region: Region, contig_lengths: Mapping[str, int]) -> Dict[str, object]:
    required_paths = [
        args.reference, f"{args.reference}.fai", args.vcf, args.vcf_index,
        args.snpsplit_snps, args.unsplit_bam, args.g1_bam, args.g2_bam,
        args.unsplit_bai, args.g1_bai, args.g2_bai,
    ]
    for path in required_paths:
        if not os.path.isfile(path):
            die(f"required file does not exist: {path}")
    for tool in (args.samtools, args.bcftools):
        if shutil.which(tool) is None and not os.path.isfile(tool):
            die(f"required executable not found: {tool}")
    for bam in (args.unsplit_bam, args.g1_bam, args.g2_bam):
        run_checked([args.samtools, "quickcheck", "-v", bam])
        sort_order = bam_sort_order(args.samtools, bam)
        if sort_order != "coordinate":
            die(f"BAM must be coordinate-sorted (SO:coordinate), observed {sort_order!r}: {bam}")
        length = bam_contig_length(args.samtools, bam, region.contig)
        if length != contig_lengths[region.contig]:
            die(
                f"contig length mismatch: {bam} has {length}, reference has "
                f"{contig_lengths[region.contig]}"
            )
    tag_checks = {
        "unsplit": sample_bam_tags(args.samtools, args.unsplit_bam, region),
        "g1": sample_bam_tags(args.samtools, args.g1_bam, region),
        "g2": sample_bam_tags(args.samtools, args.g2_bam, region),
    }
    for label, counts in tag_checks.items():
        records = counts.get("records", 0)
        if records == 0:
            die(f"no records found in {label} BAM for {region.samtools}")
        for tag in ("XM", "XR", "XG"):
            if counts.get(tag, 0) != records:
                die(f"{label} BAM: {tag} present in {counts.get(tag, 0)}/{records} sampled records")
    samples = read_vcf_samples(args.vcf)
    if args.sample not in samples:
        die(f"VCF sample {args.sample!r} not found")
    return {
        "status": "PASS",
        "tag_checks": tag_checks,
        "tool_versions": {
            "samtools": tool_version(args.samtools),
            "bcftools": tool_version(args.bcftools),
            "python": sys.version.split()[0],
        },
    }


def write_json(path: str, payload: Mapping) -> None:
    with atomic_text_writer(path) as out:
        json.dump(payload, out, indent=2, sort_keys=True)
        out.write("\n")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Discover strict 000/111 three-CpG fragment bimodality and test "
            "WGS phase-set-linked haplotype imbalance."
        )
    )
    parser.add_argument("--sample", required=True)
    parser.add_argument("--tissue", required=True)
    parser.add_argument("--region", required=True, help="CONTIG or 1-based CONTIG:START-END")
    parser.add_argument("--reference", required=True)
    parser.add_argument("--vcf", required=True)
    parser.add_argument("--vcf-index", required=True)
    parser.add_argument("--snpsplit-snps", required=True)
    parser.add_argument("--unsplit-bam", required=True)
    parser.add_argument("--unsplit-bai", required=True)
    parser.add_argument("--g1-bam", required=True)
    parser.add_argument("--g1-bai", required=True)
    parser.add_argument("--g2-bam", required=True)
    parser.add_argument("--g2-bai", required=True)
    parser.add_argument("--output-prefix", required=True)
    parser.add_argument("--results-root", required=True)
    parser.add_argument("--tmp-root", required=True)
    parser.add_argument("--stage", choices=("all", "discover", "haplotype"), default="all")
    parser.add_argument("--candidate-windows", help="Required for --stage haplotype; defaults to OUTPUT_PREFIX.candidate_windows.tsv.gz")
    parser.add_argument("--preflight-only", action="store_true")
    parser.add_argument("--samtools", default="samtools")
    parser.add_argument("--bcftools", default="bcftools")
    parser.add_argument(
        "--fragment-cpgs",
        type=int,
        choices=ALLOWED_FRAGMENT_CPGS,
        default=DEFAULT_FRAGMENT_CPGS,
        help="Must be 3; no CpG4/CpG5 mode exists in this implementation",
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=4,
        help="Threads passed to each samtools view process (default: 4)",
    )
    parser.add_argument(
        "--workers",
        type=int,
        default=1,
        help=(
            "Python worker processes for the discovery pass; QNAME-sharded and "
            "deterministic (default: 1)"
        ),
    )
    parser.add_argument("--min-mapq", type=int, default=0)
    parser.add_argument("--pair-wait-bp", type=int, default=2000)
    parser.add_argument("--fetch-padding-bp", type=int, default=2000)
    parser.add_argument("--fetch-batch-gap-bp", type=int, default=50000)
    parser.add_argument("--min-fragments", type=int, default=10)
    parser.add_argument("--min-mode-count", type=int, default=3)
    parser.add_argument("--min-minor-mode-fraction", type=float, default=0.20)
    parser.add_argument("--max-mixed-fraction", type=float, default=0.20)
    parser.add_argument("--min-link-snps", type=int, default=1)
    parser.add_argument(
        "--min-haplotype-fragments",
        type=int,
        default=REQUIRED_HAPLOTYPE_PURE_FRAGMENTS,
        help="Fixed at 5 pure U+M fragments per haplotype",
    )
    return parser


def validate_args(args: argparse.Namespace) -> None:
    if args.threads < 1:
        die("--threads must be >= 1")
    if args.workers < 1:
        die("--workers must be >= 1")
    for name in ("min_minor_mode_fraction", "max_mixed_fraction"):
        value = getattr(args, name)
        if not 0 <= value <= 1:
            die(f"--{name.replace('_', '-')} must be between 0 and 1")
    for name in (
        "min_fragments", "min_mode_count", "min_link_snps", "min_haplotype_fragments",
        "pair_wait_bp", "fetch_padding_bp", "fetch_batch_gap_bp",
    ):
        if getattr(args, name) < 0:
            die(f"--{name.replace('_', '-')} must be non-negative")
    if args.fragment_cpgs != 3:
        die("strict CpG3 analysis supports exactly --fragment-cpgs 3")
    if args.min_haplotype_fragments != REQUIRED_HAPLOTYPE_PURE_FRAGMENTS:
        die(
            "--min-haplotype-fragments is fixed at "
            f"{REQUIRED_HAPLOTYPE_PURE_FRAGMENTS} for the primary analysis"
        )
    if not any(
        math.isclose(args.max_mixed_fraction, value, rel_tol=0.0, abs_tol=1e-12)
        for value in ALLOWED_MAX_MIXED_FRACTIONS
    ):
        die("--max-mixed-fraction must be one of 0.10, 0.20, or 0.30")


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    validate_args(args)
    tmp_root = configure_tmp_root(args.tmp_root)
    results_root = os.path.realpath(os.path.abspath(args.results_root))
    prefix = os.path.realpath(os.path.abspath(args.output_prefix))
    try:
        common = os.path.commonpath((results_root, prefix))
    except ValueError:
        common = ""
    if common != results_root or prefix == results_root:
        die(f"--output-prefix must be below --results-root: {results_root}")
    os.makedirs(os.path.dirname(prefix), exist_ok=True)
    fai_path = f"{args.reference}.fai"
    if not os.path.isfile(fai_path):
        die(f"reference index not found: {fai_path}")
    contig_lengths = read_fai(fai_path)
    region = parse_region(args.region, contig_lengths)
    log(
        f"cpg3_strict_bimodal_asm v{VERSION}; sample={args.sample}; tissue={args.tissue}; "
        f"region={region.samtools}; stage={args.stage}; "
        f"fragment_cpgs={args.fragment_cpgs}; workers={args.workers}"
    )

    preflight_result = preflight(args, region, contig_lengths)
    metadata: Dict[str, object] = {
        "program": "cpg3_strict_bimodal_asm",
        "version": VERSION,
        "started_at": dt.datetime.now(dt.timezone.utc).isoformat(),
        "sample": args.sample,
        "tissue": args.tissue,
        "region": dataclasses.asdict(region),
        "path_contract": {
            "results_root": results_root,
            "tmp_root": tmp_root,
            "system_tmp_forbidden": True,
        },
        "haplotype_contract": {
            "vcf_left_allele": "genome1/H1",
            "vcf_right_allele": "genome2/H2",
        },
        "coordinate_contract": {
            "internal_and_BED": "0-based half-open",
            "VCF_SNPsplit_and_cpg_positions_column": "1-based",
        },
        "analysis_contract": {
            "fragment_cpgs": 3,
            "fragment_pattern": {
                "U": "000 only",
                "M": "111 only",
                "X": "001,010,011,100,101,110",
                "continuous_fraction_cutoffs_used": False,
            },
            "candidate_window": {
                "minimum_complete_fragments": args.min_fragments,
                "minimum_U_fragments": args.min_mode_count,
                "minimum_M_fragments": args.min_mode_count,
                "minimum_minor_mode_fraction_among_U_plus_M": args.min_minor_mode_fraction,
                "maximum_X_fraction_among_all_complete_fragments": args.max_mixed_fraction,
            },
            "requires_complete_window": True,
            "m_bias_trim_bases": 0,
            "conversion_ambiguous_link_snps": "excluded",
            "association_test": "two-sided Fisher exact test on [[H1_M,H1_U],[H2_M,H2_U]]",
            "minimum_pure_U_plus_M_fragments_per_haplotype": REQUIRED_HAPLOTYPE_PURE_FRAGMENTS,
            "X_in_primary_fisher_test": False,
            "fdr_scope": "all tested window_by_PS rows in this requested run",
        },
        "parameters": vars(args),
        "inputs": {
            name: file_fingerprint(path)
            for name, path in {
                "reference": args.reference,
                "reference_fai": fai_path,
                "vcf": args.vcf,
                "vcf_index": args.vcf_index,
                "snpsplit_snps": args.snpsplit_snps,
                "unsplit_bam": args.unsplit_bam,
                "unsplit_bai": args.unsplit_bai,
                "g1_bam": args.g1_bam,
                "g1_bai": args.g1_bai,
                "g2_bam": args.g2_bam,
                "g2_bai": args.g2_bai,
            }.items()
        },
        "preflight": preflight_result,
    }
    metadata_path = f"{prefix}.run_metadata.json"
    if args.preflight_only:
        metadata["status"] = "PREFLIGHT_PASS"
        metadata["finished_at"] = dt.datetime.now(dt.timezone.utc).isoformat()
        write_json(metadata_path, metadata)
        log(f"preflight PASS; wrote {metadata_path}")
        return 0

    extended_region = Region(
        region.contig,
        max(0, region.start0 - args.fetch_padding_bp),
        min(contig_lengths[region.contig], region.end0 + args.fetch_padding_bp),
    )
    log("loading phased VCF/SNPsplit context")
    excluded_positions0, snp_by_pos0, variant_qc = load_variant_context(
        bcftools=args.bcftools,
        vcf_path=args.vcf,
        sample=args.sample,
        snpsplit_path=args.snpsplit_snps,
        extended_region=extended_region,
        include_conversion_snps=False,
    )
    log(
        f"variant context: het SNVs={variant_qc.get('heterozygous_snv', 0):,}; "
        f"primary-safe PS-link SNPs={len(snp_by_pos0):,}"
    )
    if args.stage in {"all", "haplotype"} and not snp_by_pos0:
        die(
            "no primary-safe phased SNPs linked consistently between the VCF and "
            "SNPsplit annotation in the requested interval"
        )
    log("enumerating shared reference CpGs")
    cpg_positions0, reference_qc = enumerate_shared_reference_cpgs(
        samtools=args.samtools,
        fasta=args.reference,
        region=region,
        excluded_het_positions0=excluded_positions0,
    )
    if len(cpg_positions0) < args.fragment_cpgs:
        die(
            f"fewer than {args.fragment_cpgs} eligible shared reference CpGs "
            "in requested region"
        )
    log(f"eligible shared reference CpGs={len(cpg_positions0):,}")

    windows: List[CandidateWindow]
    candidate_regions: List[CandidateRegion]
    qc_sections: Dict[str, object] = {
        "variant": variant_qc,
        "reference": reference_qc,
    }
    candidate_path = args.candidate_windows or f"{prefix}.candidate_windows.tsv.gz"

    if args.stage in {"all", "discover"}:
        log(
            f"starting unsplit-BAM {args.fragment_cpgs}-CpG discovery pass "
            f"with {args.workers} Python worker(s)"
        )
        pattern_counts, discovery_qc = discovery_pass(
            samtools=args.samtools,
            unsplit_bam=args.unsplit_bam,
            region=region,
            cpg_positions0=cpg_positions0,
            fragment_cpgs=args.fragment_cpgs,
            threads=args.threads,
            workers=args.workers,
            min_mapq=args.min_mapq,
            pair_wait_bp=args.pair_wait_bp,
        )
        windows, candidate_regions = build_candidates(
            cpg_positions0, pattern_counts, args
        )
        write_candidates(
            prefix, args.sample, args.tissue, region.contig, windows, candidate_regions
        )
        qc_sections["discovery"] = discovery_qc
        log(
            f"discovery complete: candidate windows={len(windows):,}; "
            f"merged regions={len(candidate_regions):,}"
        )
    else:
        if not os.path.isfile(candidate_path):
            die(f"candidate window file not found for haplotype stage: {candidate_path}")
        windows, candidate_regions = load_candidates(
            candidate_path, cpg_positions0, args.fragment_cpgs
        )
        log(f"loaded candidate windows={len(windows):,}; regions={len(candidate_regions):,}")

    if args.stage in {"all", "haplotype"}:
        stats: Dict[Tuple[int, str], List[int]] = {}
        hap_qc: MutableMapping[str, int] = collections.Counter()
        if windows:
            log("collecting H1/genome1 candidate-window fragments")
            collect_haplotype_counts_for_bam(
                source_haplotype="H1", samtools=args.samtools, bam=args.g1_bam,
                analysis_region=region, fetch_bounds=extended_region,
                candidate_regions=candidate_regions,
                cpg_positions0=cpg_positions0, fragment_cpgs=args.fragment_cpgs,
                snp_by_pos0=snp_by_pos0, stats=stats,
                threads=args.threads, min_mapq=args.min_mapq,
                pair_wait_bp=args.pair_wait_bp, fetch_padding_bp=args.fetch_padding_bp,
                fetch_batch_gap_bp=args.fetch_batch_gap_bp,
                min_link_snps=args.min_link_snps, qc=hap_qc,
            )
            log("collecting H2/genome2 candidate-window fragments")
            collect_haplotype_counts_for_bam(
                source_haplotype="H2", samtools=args.samtools, bam=args.g2_bam,
                analysis_region=region, fetch_bounds=extended_region,
                candidate_regions=candidate_regions,
                cpg_positions0=cpg_positions0, fragment_cpgs=args.fragment_cpgs,
                snp_by_pos0=snp_by_pos0, stats=stats,
                threads=args.threads, min_mapq=args.min_mapq,
                pair_wait_bp=args.pair_wait_bp, fetch_padding_bp=args.fetch_padding_bp,
                fetch_batch_gap_bp=args.fetch_batch_gap_bp,
                min_link_snps=args.min_link_snps, qc=hap_qc,
            )
        windows_by_index = {window.index: window for window in windows}
        write_haplotype_stats(
            prefix=prefix, sample=args.sample, tissue=args.tissue,
            contig=region.contig, windows_by_index=windows_by_index,
            stats=stats, min_haplotype_fragments=args.min_haplotype_fragments,
        )
        hap_qc["window_ps_rows"] = len(stats)
        qc_sections["haplotype"] = dict(hap_qc)
        log(f"haplotype stage complete: window-by-PS rows={len(stats):,}")

    qc_path = f"{prefix}.qc.json"
    write_json(qc_path, qc_sections)
    metadata["qc"] = qc_sections
    metadata["candidate_windows"] = len(windows)
    metadata["candidate_regions"] = len(candidate_regions)
    metadata["status"] = "PASS"
    metadata["finished_at"] = dt.datetime.now(dt.timezone.utc).isoformat()
    write_json(metadata_path, metadata)
    log(f"PASS; metadata={metadata_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
