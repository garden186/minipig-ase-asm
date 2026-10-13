"""Conversion-aware adapter around the preserved exact-CpG4 counting engine."""
from __future__ import annotations
from collections import Counter
from pathlib import Path
import sys
import traceback
import types

from support import PACKAGE, rss_kib
sys.path.insert(0, str(PACKAGE / "vendor"))
import cpgk_scan_helpers as helpers
import count_exact_cpgk as exact

CORE_PATH = PACKAGE / "vendor/bimodal_fragment_asm.py"
CORE = helpers.load_core(CORE_PATH)


def allele_support(base, h1, h2, xg):
    matches = [i for i, allele in enumerate((h1, h2))
               if base == allele or (xg == "CT" and allele == "C" and base == "T")
               or (xg == "GA" and allele == "G" and base == "A")]
    return matches[0] if len(matches) == 1 else None


def conversion_core(snps, qc):
    proxy = types.SimpleNamespace(**vars(CORE))

    def extract(record, positions):
        calls = CORE.extract_snp_calls(record, positions)
        xg = record.tags.get("XG")
        if xg not in {"CT", "GA"}:
            raise ValueError("Missing or invalid XG on a usable alignment")
        for pos, base in calls.items():
            info = snps[pos]
            support = allele_support(base, info.h1, info.h2, xg)
            if support is not None:
                allele = (info.h1, info.h2)[support]
                if base != allele:
                    qc["conversion_normalized_read_site_observations"] += 1
                calls[pos] = allele
            # Nonallelic bases must survive until mate merging to preserve true conflicts.
        return calls

    proxy.extract_snp_calls = extract
    return proxy


def checked_records(lines):
    for line in lines:
        record = CORE.parse_sam_line(line)
        if CORE.record_is_usable(record, 0):
            if record.tags.get("XG") not in {"CT", "GA"}:
                raise ValueError("Missing/invalid XG: " + record.qname)
            if len(record.tags.get("XM", "")) != len(record.seq) or record.seq == "*":
                raise ValueError("Missing/invalid XM or SEQ: " + record.qname)
        yield record


def converted_worker(worker_id, input_queue, result_queue, core_path, cpg_positions0,
                     disrupted, snp_by_pos0, k_values, min_mapq, pair_wait_bp,
                     min_link_snps, flush_fragments):
    """Apply the adapter in every child, including spawn-based parity tests."""
    try:
        def lines():
            while True:
                item = input_queue.get()
                if item == "END":
                    return
                yield from item
        qc = Counter()
        core = conversion_core(snp_by_pos0, qc)
        counts = {}
        chunks = 0

        def flush(values):
            nonlocal chunks
            result_queue.put(("COUNTS", worker_id, exact.encode_counts(values)))
            chunks += 1

        exact.count_fragment_records(core=core, records=checked_records(lines()),
            cpg_positions0=cpg_positions0, disrupted=disrupted, snp_by_pos0=snp_by_pos0,
            k_values=k_values, min_mapq=min_mapq, pair_wait_bp=pair_wait_bp,
            min_link_snps=min_link_snps, counts=counts, qc=qc,
            flush_fragments=flush_fragments, flush_callback=flush)
        qc[f"worker_{worker_id}_peak_rss_kib"] = rss_kib()
        result_queue.put(("DONE", worker_id, dict(qc), chunks))
    except BaseException:
        result_queue.put(("ERROR", worker_id, traceback.format_exc()))


def count(core, samtools, bam, region, threads, workers, positions, disrupted, snps):
    counts, qc = {}, Counter()
    if workers > 1:
        original = exact.exact_count_worker
        exact.exact_count_worker = converted_worker
        try:
            exact.parallel_count_bam(core=core, core_path=str(CORE_PATH), samtools=samtools,
                bam=bam, region=region, threads=threads, workers=workers,
                worker_flush_fragments=250000, cpg_positions0=positions, disrupted=disrupted,
                snp_by_pos0=snps, k_values=[4], min_mapq=0, pair_wait_bp=2000,
                min_link_snps=1, counts=counts, qc=qc)
        finally:
            exact.exact_count_worker = original
    else:
        proxy = conversion_core(snps, qc)
        lines = helpers.sam_record_lines(core=core, samtools=samtools, bam=bam,
                                         region=region, threads=threads)
        exact.count_fragment_records(core=proxy, records=checked_records(lines),
            cpg_positions0=positions, disrupted=disrupted, snp_by_pos0=snps,
            k_values=[4], min_mapq=0, pair_wait_bp=2000, min_link_snps=1,
            counts=counts, qc=qc)
    return counts, dict(qc)
