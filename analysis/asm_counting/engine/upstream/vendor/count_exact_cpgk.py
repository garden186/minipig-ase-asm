#!/usr/bin/env python3
"""Count exact reference-consecutive CpG-K patterns by PS and haplotype."""

from __future__ import annotations

import argparse
import bisect
import collections
import multiprocessing as mp
import pickle
import queue
import sys
import traceback
import zlib
from pathlib import Path
from typing import Callable, Mapping, MutableMapping, Sequence

from cpgk_common import (
    COUNT_FIELDS,
    MAX_K,
    CpGKError,
    VERSION,
    exact_window_id,
    load_reference_cpgs,
    parse_positive_ints,
    write_json,
    write_tsv,
)
from cpgk_scan_helpers import (
    disrupted_reference_cpgs,
    load_core,
    merged_calls,
    sam_record_lines,
    stream_record_groups,
)


WORKER_BATCH_RECORDS = 4096
DEFAULT_WORKER_FLUSH_FRAGMENTS = 250_000


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--tissue", required=True)
    parser.add_argument("--chrom", required=True)
    parser.add_argument("--reference", required=True)
    parser.add_argument("--reference-cpgs", required=True)
    parser.add_argument("--genotype-vcf", required=True)
    parser.add_argument("--phased-vcf", required=True)
    parser.add_argument("--snpsplit-snps", required=True)
    parser.add_argument("--unsplit-bam", required=True)
    parser.add_argument("--output-prefix", required=True)
    parser.add_argument("--core", required=True)
    parser.add_argument("--k-values", default="3,4,5")
    parser.add_argument("--samtools", default="samtools")
    parser.add_argument("--bcftools", default="bcftools")
    parser.add_argument("--threads", type=int, default=2)
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument(
        "--worker-flush-fragments",
        type=int,
        default=DEFAULT_WORKER_FLUSH_FRAGMENTS,
        help="fragment groups per incremental worker count payload",
    )
    parser.add_argument("--min-mapq", type=int, default=0)
    parser.add_argument("--pair-wait-bp", type=int, default=2000)
    parser.add_argument("--min-link-snps", type=int, default=1)
    parser.add_argument("--force", action="store_true")
    return parser


def exact_windows_from_calls(
    calls_by_index: dict[int, int],
    k_values: Sequence[int],
    reference_positions0: Sequence[int],
):
    for first in sorted(calls_by_index):
        for k in k_values:
            last_exclusive = first + k
            if last_exclusive > len(reference_positions0):
                continue
            target = range(first, last_exclusive)
            if all(index in calls_by_index for index in target):
                yield k, first, sum(calls_by_index[index] for index in target)


def merge_counts(
    destination: MutableMapping[tuple[int, int, str], list[int]],
    source: Mapping[tuple[int, int, str], Sequence[int]],
) -> None:
    width = 2 * (MAX_K + 1)
    for key, source_counts in source.items():
        if len(source_counts) != width:
            raise CpGKError(f"invalid parallel count width for {key}: {len(source_counts)}")
        target = destination.setdefault(key, [0] * width)
        for index, value in enumerate(source_counts):
            target[index] += int(value)


def count_one_fragment(
    *,
    core,
    fragment_records: Sequence[object],
    cpg_positions0: Sequence[int],
    disrupted: Sequence[int],
    snp_by_pos0: Mapping[int, object],
    sorted_snp_positions0: Sequence[int],
    k_values: Sequence[int],
    min_link_snps: int,
    counts: MutableMapping[tuple[int, int, str], list[int]],
    qc: MutableMapping[str, int],
) -> None:
    qc["all_usable_fragments"] += 1
    if len(fragment_records) > 1:
        qc["paired_usable_fragments"] += 1

    snp_calls, snp_conflicts = merged_calls(
        core, fragment_records, core.extract_snp_calls, sorted_snp_positions0
    )
    qc["overlap_snp_conflicts"] += snp_conflicts
    phase_fragment = core.Fragment(
        fragment_records[0].qname,
        {},
        snp_calls,
        records=len(fragment_records),
        snp_conflicts=snp_conflicts,
    )
    states, votes = core.infer_ps_states(phase_fragment, snp_by_pos0, min_link_snps)
    if len(votes) > 1:
        qc["phase_excluded_multiple_ps_with_votes"] += 1
        return
    if len(states) != 1:
        qc["phase_not_uniquely_assigned"] += 1
        return
    phase_set, haplotype = next(iter(states.items()))
    if haplotype not in {"H1", "H2"}:
        qc["phase_invalid_haplotype"] += 1
        return
    qc["phase_assigned_fragments"] += 1
    qc[f"phase_{haplotype}_fragments"] += 1

    raw_calls, cpg_conflicts = merged_calls(
        core, fragment_records, core.extract_cpg_calls
    )
    qc["overlap_cpg_conflicts"] += cpg_conflicts
    calls_by_index: dict[int, int] = {}
    for pos0, call in raw_calls.items():
        cpg_index = bisect.bisect_left(cpg_positions0, pos0)
        if cpg_index >= len(cpg_positions0) or cpg_positions0[cpg_index] != pos0:
            qc["methylation_calls_not_reference_cpg"] += 1
            continue
        if disrupted[cpg_index]:
            qc["methylation_calls_genotype_disrupted"] += 1
            continue
        calls_by_index[cpg_index] = int(call)
    if not calls_by_index:
        return
    qc["phase_fragments_with_reference_cpg"] += 1

    hap_offset = 0 if haplotype == "H1" else MAX_K + 1
    fragment_windows = 0
    for k, first, n_methylated in exact_windows_from_calls(
        calls_by_index, k_values, cpg_positions0
    ):
        key = (k, first, str(phase_set))
        value = counts.setdefault(key, [0] * (2 * (MAX_K + 1)))
        value[hap_offset + n_methylated] += 1
        qc[f"K{k}_fragment_window_observations"] += 1
        fragment_windows += 1
    if fragment_windows:
        qc["fragments_with_exact_cpgk_window"] += 1
        qc["exact_cpgk_window_observations"] += fragment_windows


def count_fragment_records(
    *,
    core,
    records,
    cpg_positions0: Sequence[int],
    disrupted: Sequence[int],
    snp_by_pos0: Mapping[int, object],
    k_values: Sequence[int],
    min_mapq: int,
    pair_wait_bp: int,
    min_link_snps: int,
    counts: MutableMapping[tuple[int, int, str], list[int]],
    qc: MutableMapping[str, int],
    flush_fragments: int = 0,
    flush_callback: Callable[[Mapping[tuple[int, int, str], Sequence[int]]], None] | None = None,
) -> None:
    sorted_snp_positions0 = sorted(snp_by_pos0)
    groups = stream_record_groups(
        core=core,
        records=records,
        min_mapq=min_mapq,
        pair_wait_bp=pair_wait_bp,
        qc=qc,
    )
    since_flush = 0
    for fragment_records in groups:
        count_one_fragment(
            core=core,
            fragment_records=fragment_records,
            cpg_positions0=cpg_positions0,
            disrupted=disrupted,
            snp_by_pos0=snp_by_pos0,
            sorted_snp_positions0=sorted_snp_positions0,
            k_values=k_values,
            min_link_snps=min_link_snps,
            counts=counts,
            qc=qc,
        )
        since_flush += 1
        if flush_callback is not None and flush_fragments and since_flush >= flush_fragments:
            if counts:
                flush_callback(counts)
                counts.clear()
            since_flush = 0
    if flush_callback is not None and counts:
        flush_callback(counts)
        counts.clear()


def encode_counts(counts: Mapping[tuple[int, int, str], Sequence[int]]) -> bytes:
    return zlib.compress(pickle.dumps(counts, protocol=5), level=1)


def decode_counts(payload: bytes) -> dict[tuple[int, int, str], list[int]]:
    value = pickle.loads(zlib.decompress(payload))
    if not isinstance(value, dict):
        raise CpGKError("parallel worker returned a non-dictionary count payload")
    return value


def exact_count_worker(
    worker_id: int,
    input_queue,
    result_queue,
    core_path: str,
    cpg_positions0: Sequence[int],
    disrupted: Sequence[int],
    snp_by_pos0: Mapping[int, object],
    k_values: Sequence[int],
    min_mapq: int,
    pair_wait_bp: int,
    min_link_snps: int,
    flush_fragments: int,
) -> None:
    try:
        core = load_core(Path(core_path))

        def records():
            while True:
                item = input_queue.get()
                if item == "END":
                    return
                for line in item:
                    yield core.parse_sam_line(line)

        counts: dict[tuple[int, int, str], list[int]] = {}
        qc: MutableMapping[str, int] = collections.Counter()
        chunks_sent = 0

        def flush(local_counts: Mapping[tuple[int, int, str], Sequence[int]]) -> None:
            nonlocal chunks_sent
            result_queue.put(("COUNTS", worker_id, encode_counts(local_counts)))
            chunks_sent += 1

        count_fragment_records(
            core=core,
            records=records(),
            cpg_positions0=cpg_positions0,
            disrupted=disrupted,
            snp_by_pos0=snp_by_pos0,
            k_values=k_values,
            min_mapq=min_mapq,
            pair_wait_bp=pair_wait_bp,
            min_link_snps=min_link_snps,
            counts=counts,
            qc=qc,
            flush_fragments=flush_fragments,
            flush_callback=flush,
        )
        result_queue.put(("DONE", worker_id, dict(qc), chunks_sent))
    except BaseException:
        result_queue.put(("ERROR", worker_id, traceback.format_exc()))


def put_worker_item(input_queue, item, process: mp.Process, drain_results: Callable[[], None]) -> None:
    while True:
        try:
            input_queue.put(item, timeout=0.5)
            return
        except queue.Full:
            drain_results()
            if not process.is_alive():
                raise CpGKError(
                    f"exact-count worker {process.name} exited before consuming its input"
                )


def parallel_count_bam(
    *,
    core,
    core_path: str,
    samtools: str,
    bam: str,
    region,
    threads: int,
    workers: int,
    worker_flush_fragments: int,
    cpg_positions0: Sequence[int],
    disrupted: Sequence[int],
    snp_by_pos0: Mapping[int, object],
    k_values: Sequence[int],
    min_mapq: int,
    pair_wait_bp: int,
    min_link_snps: int,
    counts: MutableMapping[tuple[int, int, str], list[int]],
    qc: MutableMapping[str, int],
) -> None:
    if "fork" not in mp.get_all_start_methods():
        raise CpGKError("--workers > 1 requires multiprocessing fork support")
    context = mp.get_context("fork")
    input_queues = [context.Queue(maxsize=4) for _ in range(workers)]
    result_queue = context.Queue(maxsize=max(4, workers * 2))
    processes = [
        context.Process(
            target=exact_count_worker,
            name=f"exact-cpgk-{worker_id}",
            args=(
                worker_id,
                input_queues[worker_id],
                result_queue,
                core_path,
                cpg_positions0,
                disrupted,
                snp_by_pos0,
                tuple(k_values),
                min_mapq,
                pair_wait_bp,
                min_link_snps,
                worker_flush_fragments,
            ),
        )
        for worker_id in range(workers)
    ]
    for process in processes:
        process.start()

    done_workers: set[int] = set()
    chunks_received = 0

    def handle_result(result) -> None:
        nonlocal chunks_received
        status, worker_id, *payload = result
        if status == "COUNTS":
            merge_counts(counts, decode_counts(payload[0]))
            chunks_received += 1
        elif status == "DONE":
            if worker_id in done_workers:
                raise CpGKError(f"duplicate completion from worker {worker_id}")
            local_qc, chunks_sent = payload
            qc.update(local_qc)
            qc["parallel_count_chunks_sent"] += int(chunks_sent)
            done_workers.add(worker_id)
        elif status == "ERROR":
            raise CpGKError(f"parallel exact-count worker {worker_id} failed:\n{payload[0]}")
        else:
            raise CpGKError(f"unknown parallel worker message: {status}")

    def drain_results() -> None:
        while True:
            try:
                handle_result(result_queue.get_nowait())
            except queue.Empty:
                return

    batches: list[list[str]] = [[] for _ in range(workers)]
    records_dispatched = 0
    success = False
    try:
        for line in sam_record_lines(
            core=core, samtools=samtools, bam=bam, region=region, threads=threads
        ):
            qname = line.split("\t", 1)[0]
            if not qname or qname == line:
                raise CpGKError("SAM record lacks a QNAME delimiter")
            worker_id = zlib.crc32(qname.encode("utf-8")) % workers
            batches[worker_id].append(line)
            records_dispatched += 1
            if len(batches[worker_id]) >= WORKER_BATCH_RECORDS:
                put_worker_item(
                    input_queues[worker_id], batches[worker_id], processes[worker_id], drain_results
                )
                batches[worker_id] = []
            if records_dispatched % 10_000 == 0:
                drain_results()
            if records_dispatched % 5_000_000 == 0:
                print(
                    f"[SCAN] {bam} chr{region.contig}: dispatched {records_dispatched:,} SAM records",
                    file=sys.stderr,
                    flush=True,
                )
        for worker_id in range(workers):
            if batches[worker_id]:
                put_worker_item(
                    input_queues[worker_id], batches[worker_id], processes[worker_id], drain_results
                )
            put_worker_item(input_queues[worker_id], "END", processes[worker_id], drain_results)

        while len(done_workers) < workers:
            try:
                handle_result(result_queue.get(timeout=1.0))
            except queue.Empty:
                failed = [
                    process for process in processes
                    if process.exitcode is not None and process.exitcode != 0
                ]
                if failed:
                    raise CpGKError(
                        "parallel exact-count worker failed: "
                        + ", ".join(f"{process.name} exit={process.exitcode}" for process in failed)
                    )
        drain_results()
        for process in processes:
            process.join()
            if process.exitcode != 0:
                raise CpGKError(
                    f"parallel exact-count worker {process.name} exited with {process.exitcode}"
                )
        if chunks_received != int(qc["parallel_count_chunks_sent"]):
            raise CpGKError(
                f"parallel count chunk mismatch: received={chunks_received}; "
                f"sent={qc['parallel_count_chunks_sent']}"
            )
        if int(qc["sam_records_seen"]) != records_dispatched:
            raise CpGKError(
                f"parallel SAM record mismatch: parsed={qc['sam_records_seen']}; "
                f"dispatched={records_dispatched}"
            )
        qc["python_worker_processes"] = workers
        qc["sam_records_dispatched"] = records_dispatched
        qc["parallel_count_chunks_received"] = chunks_received
        success = True
    finally:
        if not success:
            for process in processes:
                if process.is_alive():
                    process.terminate()
            for process in processes:
                process.join()
        for input_queue in input_queues:
            if not success:
                input_queue.cancel_join_thread()
            input_queue.close()
            if success:
                input_queue.join_thread()
        if not success:
            result_queue.cancel_join_thread()
        result_queue.close()
        if success:
            result_queue.join_thread()


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    k_values = parse_positive_ints(args.k_values, label="K values")
    if min(k_values) < 3 or max(k_values) > MAX_K:
        raise CpGKError(f"K values must be between 3 and {MAX_K}")
    if (
        args.threads < 1
        or args.workers < 1
        or args.worker_flush_fragments < 1
        or args.min_link_snps < 1
        or args.pair_wait_bp < 0
    ):
        raise CpGKError("invalid scan resource or phase-link parameter")

    output_prefix = Path(args.output_prefix)
    count_path = Path(f"{output_prefix}.exact_cpgk_counts.tsv.gz")
    metadata_path = Path(f"{output_prefix}.exact_cpgk_metadata.json")
    if count_path.exists() and metadata_path.exists() and not args.force:
        print(f"[SKIP] completed exact CpG-K counts exist: {count_path}", file=sys.stderr)
        return 0

    core = load_core(Path(args.core))
    requested_chrom = args.chrom.removeprefix("chr")
    chrom, cpg_positions0 = load_reference_cpgs(Path(args.reference_cpgs))
    if chrom != requested_chrom:
        raise CpGKError("reference-CpG chromosome does not match request")

    contig_lengths = core.read_fai(f"{args.reference}.fai")
    analysis_region = core.parse_region(requested_chrom, contig_lengths)
    disrupted, genotype_qc = disrupted_reference_cpgs(
        core=core,
        bcftools=args.bcftools,
        genotype_vcf=args.genotype_vcf,
        sample=args.sample,
        region=analysis_region,
        cpg_positions0=cpg_positions0,
    )
    _excluded, snp_by_pos0, phased_qc = core.load_variant_context(
        bcftools=args.bcftools,
        vcf_path=args.phased_vcf,
        sample=args.sample,
        snpsplit_path=args.snpsplit_snps,
        extended_region=analysis_region,
        include_conversion_snps=False,
    )
    qc: MutableMapping[str, int] = collections.Counter()
    # key: K, first reference-CpG index, WGS phase set
    # value: H1_nM0..5,H2_nM0..5
    counts: dict[tuple[int, int, str], list[int]] = {}

    if args.workers > 1:
        parallel_count_bam(
            core=core,
            core_path=args.core,
            samtools=args.samtools,
            bam=args.unsplit_bam,
            region=analysis_region,
            threads=args.threads,
            workers=args.workers,
            worker_flush_fragments=args.worker_flush_fragments,
            cpg_positions0=cpg_positions0,
            disrupted=disrupted,
            snp_by_pos0=snp_by_pos0,
            k_values=k_values,
            min_mapq=args.min_mapq,
            pair_wait_bp=args.pair_wait_bp,
            min_link_snps=args.min_link_snps,
            counts=counts,
            qc=qc,
        )
    else:
        lines = sam_record_lines(
            core=core,
            samtools=args.samtools,
            bam=args.unsplit_bam,
            region=analysis_region,
            threads=args.threads,
        )
        records = (core.parse_sam_line(line) for line in lines)
        count_fragment_records(
            core=core,
            records=records,
            cpg_positions0=cpg_positions0,
            disrupted=disrupted,
            snp_by_pos0=snp_by_pos0,
            k_values=k_values,
            min_mapq=args.min_mapq,
            pair_wait_bp=args.pair_wait_bp,
            min_link_snps=args.min_link_snps,
            counts=counts,
            qc=qc,
        )
        qc["python_worker_processes"] = 1

    def rows():
        for (k, first, phase_set), values in sorted(counts.items()):
            last = first + k - 1
            row = {
                "sample": args.sample,
                "tissue": args.tissue,
                "chrom": requested_chrom,
                "k": k,
                "window_id": exact_window_id(requested_chrom, cpg_positions0, first, k),
                "first_cpg_index": first,
                "last_cpg_index": last,
                "cpg_positions_1based": ",".join(
                    str(cpg_positions0[index] + 1) for index in range(first, last + 1)
                ),
                "start0": cpg_positions0[first],
                "end0": cpg_positions0[last] + 2,
                "span_bp": cpg_positions0[last] - cpg_positions0[first] + 2,
                "ps": phase_set,
            }
            for hap_index, hap in enumerate(("H1", "H2")):
                offset = hap_index * (MAX_K + 1)
                for n_methylated in range(MAX_K + 1):
                    row[f"{hap}_nM{n_methylated}"] = values[offset + n_methylated]
                row[f"{hap}_all"] = sum(values[offset : offset + k + 1])
            yield row

    write_tsv(count_path, COUNT_FIELDS, rows())
    write_json(
        metadata_path,
        {
            "program": "count_exact_cpgk",
            "version": VERSION,
            "status": "PASS",
            "sample": args.sample,
            "tissue": args.tissue,
            "chrom": requested_chrom,
            "analysis_contract": {
                "coordinate_unit": "reference-consecutive exact CpG-K window",
                "fragment_requirement": "all K reference CpGs have non-conflicting calls",
                "phase_assignment": "one safe WGS PS and one H1/H2 state required",
                "raw_state": "number of methylated target CpGs; U/X/M is downstream",
                "parallelization": "deterministic CRC32 QNAME sharding; additive raw-count merge",
                "split_half_used": False,
                "balance_filter_used": False,
            },
            "parameters": vars(args),
            "genotype_qc": genotype_qc,
            "phased_qc": phased_qc,
            "qc": dict(qc),
            "observed_window_ps_rows": len(counts),
            "output": str(count_path),
        },
    )
    print(
        f"[PASS] {args.sample} {args.tissue} chr{requested_chrom}: "
        f"{len(counts):,} exact window-PS rows",
        file=sys.stderr,
    )
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except CpGKError as error:
        print(f"[ERROR] {error}", file=sys.stderr)
        raise SystemExit(2)
