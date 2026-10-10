#!/usr/bin/env python3
"""
Collect Stage 2 RNA-seq QC for one minipig sample.

The script recursively discovers and summarizes:
  * Trim Galore/Cutadapt trimming reports (R1 and R2, all tissues)
  * STAR Log.final.out files
  * Picard/GATK MarkDuplicates metrics
  * the pipeline alignment_summary.tsv, when present
  * optionally, WASP vW tags in deduplicated BAMs (slow; requires samtools)

No input file is modified.

Example
-------
python stage2_qc_collect.py \
  --project-dir /path/to/project \
  --sample 0326

Optional exact WASP-tag scan (can take a long time):
python stage2_qc_collect.py \
  --project-dir /path/to/project \
  --sample 0326 \
  --scan-vw --threads 4
"""

from __future__ import annotations

import argparse
import csv
import re
import shutil
import subprocess
import sys
from collections import defaultdict
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple


TRIM_FIELDS = [
    "sample",
    "tissue",
    "mate",
    "report",
    "trimming_mode",
    "quality_cutoff",
    "min_length_bp",
    "adapter_overlap_bp",
    "total_reads",
    "adapter_reads",
    "adapter_reads_pct",
    "cutadapt_reads_written",
    "cutadapt_reads_written_pct",
    "short_reads_removed",
    "short_reads_removed_pct",
    "estimated_final_reads",
    "estimated_final_reads_pct",
    "total_bases",
    "quality_trimmed_bases",
    "quality_trimmed_bases_pct",
    "written_bases",
    "written_bases_pct",
    "polyg_reads",
    "polyg_reads_pct",
    "adapter_calls_1_2bp",
    "adapter_calls_1_2bp_pct_of_adapter_calls",
    "adapter_calls_ge3bp",
    "adapter_calls_ge3bp_pct_of_total_reads",
    "adapter_calls_ge5bp",
    "adapter_calls_ge5bp_pct_of_total_reads",
]


ALIGN_FIELDS = [
    "sample",
    "tissue",
    "star_log",
    "input_reads",
    "average_input_read_length",
    "uniquely_mapped_reads",
    "uniquely_mapped_pct",
    "multimapped_reads",
    "multimapped_pct",
    "too_many_loci_reads",
    "too_many_loci_pct",
    "unmapped_too_short_reads",
    "unmapped_too_short_pct",
    "unmapped_other_reads",
    "unmapped_other_pct",
    "mismatch_rate_pct",
    "total_splices",
    "annotated_splices",
    "noncanonical_splices",
    "duplicate_metrics",
    "read_pairs_examined",
    "read_pair_duplicates",
    "optical_duplicate_pairs",
    "percent_duplication",
    "estimated_library_size",
    "dedup_bam_reads",
    "phaseready_reads",
    "phaseready_pct",
    "vw_bam",
    "vw_total_alignments",
    "vw_primary_alignments",
    "vw_pass_1",
    "vw_fail_2_7",
    "vw_no_tag",
    "vw_other",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Collect Trim Galore, STAR, duplicate and optional WASP-tag QC."
    )
    parser.add_argument("--project-dir", required=True, type=Path)
    parser.add_argument("--sample", required=True)
    parser.add_argument(
        "--outdir",
        type=Path,
        default=None,
        help=(
            "Output directory (default: "
            "<project>/results/qc/rnaseq/stage2_qc/<sample>)"
        ),
    )
    parser.add_argument(
        "--scan-vw",
        action="store_true",
        help="Stream all sample dedup BAMs with samtools and count vW tags (slow).",
    )
    parser.add_argument("--samtools", default="samtools")
    parser.add_argument("--threads", type=int, default=4)
    return parser.parse_args()


def read_text(path: Path) -> str:
    return path.read_text(encoding="utf-8", errors="replace")


def numeric(value: object) -> str:
    """Return a TSV-friendly representation without converting missing to zero."""
    if value is None:
        return ""
    if isinstance(value, float):
        return f"{value:.6f}".rstrip("0").rstrip(".")
    return str(value)


def int_clean(value: str) -> Optional[int]:
    match = re.search(r"[-+]?\d[\d,]*", value)
    return int(match.group(0).replace(",", "")) if match else None


def float_clean(value: str) -> Optional[float]:
    match = re.search(r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)", value.replace(",", ""))
    return float(match.group(0)) if match else None


def pct_clean(value: str) -> Optional[float]:
    match = re.search(r"\(([-+]?\d+(?:\.\d+)?)%\)", value)
    if not match:
        match = re.search(r"([-+]?\d+(?:\.\d+)?)%", value)
    return float(match.group(1)) if match else None


def line_value(text: str, label: str) -> Optional[str]:
    for line in text.splitlines():
        if line.strip().startswith(label):
            return line.split(":", 1)[1].strip() if ":" in line else line.strip()
    return None


def safe_pct(numerator: Optional[int], denominator: Optional[int]) -> Optional[float]:
    if numerator is None or denominator in (None, 0):
        return None
    return 100.0 * numerator / denominator


def infer_trim_identity(path: Path, sample: str) -> Tuple[Optional[str], Optional[str]]:
    name = path.name
    patterns = [
        rf"^{re.escape(sample)}[-_](.+?)_R_([12])\.fastq\.gz_trimming_report\.txt$",
        rf"^{re.escape(sample)}[-_](.+?)[-_]?R([12]).*trimming_report\.txt$",
        rf"^{re.escape(sample)}[-_](.+?)_val_([12]).*trimming_report\.txt$",
    ]
    for pattern in patterns:
        match = re.match(pattern, name, flags=re.IGNORECASE)
        if match:
            return match.group(1), match.group(2)
    return None, None


def parse_adapter_lengths(text: str) -> Dict[int, int]:
    """Aggregate Cutadapt adapter-removal length tables across adapter sections."""
    counts: Dict[int, int] = defaultdict(int)
    lines = text.splitlines()
    in_table = False
    for line in lines:
        stripped = line.strip()
        if stripped.lower().startswith("length") and "count" in stripped.lower():
            in_table = True
            continue
        if not in_table:
            continue
        match = re.match(r"^(\d+)\s+([\d,]+)(?:\s+|$)", stripped)
        if match:
            counts[int(match.group(1))] += int(match.group(2).replace(",", ""))
            continue
        if stripped == "" or stripped.startswith("===") or stripped.startswith("RUN "):
            in_table = False
    return dict(counts)


def parse_trim_report(path: Path, sample: str) -> Optional[Dict[str, object]]:
    tissue, mate = infer_trim_identity(path, sample)
    if tissue is None:
        return None
    text = read_text(path)

    total_reads_s = line_value(text, "Total reads processed")
    adapter_reads_s = line_value(text, "Reads with adapters")
    reads_written_s = line_value(text, "Reads written (passing filters)")
    total_bases_s = line_value(text, "Total basepairs processed")
    quality_bases_s = line_value(text, "Quality-trimmed")
    written_bases_s = line_value(text, "Total written (filtered)")
    short_s = line_value(text, "Sequences removed because they became shorter")
    polyg_s = line_value(text, "Reads with poly-G/C trimmed")

    total_reads = int_clean(total_reads_s or "")
    adapter_reads = int_clean(adapter_reads_s or "")
    reads_written = int_clean(reads_written_s or "")
    total_bases = int_clean(total_bases_s or "")
    quality_bases = int_clean(quality_bases_s or "")
    written_bases = int_clean(written_bases_s or "")
    short_reads = int_clean(short_s or "")
    polyg_reads = int_clean(polyg_s or "")
    final_reads = None
    if total_reads is not None and short_reads is not None:
        final_reads = max(0, total_reads - short_reads)

    adapter_lengths = parse_adapter_lengths(text)
    calls_1_2 = sum(v for k, v in adapter_lengths.items() if k <= 2)
    calls_ge3 = sum(v for k, v in adapter_lengths.items() if k >= 3)
    calls_ge5 = sum(v for k, v in adapter_lengths.items() if k >= 5)
    all_calls = sum(adapter_lengths.values())

    quality_cutoff = None
    quality_line = line_value(text, "Quality Phred score cutoff")
    if quality_line:
        quality_cutoff = int_clean(quality_line)
    min_length = None
    for line in text.splitlines():
        if line.startswith("Minimum required sequence length"):
            min_length = int_clean(line.split(":", 1)[-1])
            break
    overlap = None
    overlap_line = line_value(text, "Minimum required adapter overlap (stringency)")
    if overlap_line:
        overlap = int_clean(overlap_line)
    mode = line_value(text, "Trimming mode")

    return {
        "sample": sample,
        "tissue": tissue,
        "mate": mate,
        "report": str(path),
        "trimming_mode": mode,
        "quality_cutoff": quality_cutoff,
        "min_length_bp": min_length,
        "adapter_overlap_bp": overlap,
        "total_reads": total_reads,
        "adapter_reads": adapter_reads,
        "adapter_reads_pct": pct_clean(adapter_reads_s or ""),
        "cutadapt_reads_written": reads_written,
        "cutadapt_reads_written_pct": pct_clean(reads_written_s or ""),
        "short_reads_removed": short_reads,
        "short_reads_removed_pct": pct_clean(short_s or ""),
        "estimated_final_reads": final_reads,
        "estimated_final_reads_pct": safe_pct(final_reads, total_reads),
        "total_bases": total_bases,
        "quality_trimmed_bases": quality_bases,
        "quality_trimmed_bases_pct": safe_pct(quality_bases, total_bases),
        "written_bases": written_bases,
        "written_bases_pct": (
            pct_clean(written_bases_s or "")
            if written_bases_s
            else safe_pct(written_bases, total_bases)
        ),
        "polyg_reads": polyg_reads,
        "polyg_reads_pct": pct_clean(polyg_s or ""),
        "adapter_calls_1_2bp": calls_1_2 if all_calls else None,
        "adapter_calls_1_2bp_pct_of_adapter_calls": safe_pct(calls_1_2, all_calls),
        "adapter_calls_ge3bp": calls_ge3 if all_calls else None,
        "adapter_calls_ge3bp_pct_of_total_reads": safe_pct(calls_ge3, total_reads),
        "adapter_calls_ge5bp": calls_ge5 if all_calls else None,
        "adapter_calls_ge5bp_pct_of_total_reads": safe_pct(calls_ge5, total_reads),
    }


def discover_trim_reports(project: Path, sample: str) -> Tuple[List[Dict[str, object]], List[str]]:
    roots = [
        project / "results" / "qc" / "rnaseq" / "trim",
        project / "results" / "rnaseq" / "trim",
    ]
    candidates: List[Path] = []
    for root in roots:
        if root.is_dir():
            candidates.extend(root.rglob(f"*{sample}*trimming_report.txt"))

    # Reports are often copied from the trim tree into the QC tree. Keep one
    # report per tissue/mate, preferring the QC copy.
    chosen: Dict[Tuple[str, str], Path] = {}
    warnings: List[str] = []
    for path in sorted(set(candidates)):
        tissue, mate = infer_trim_identity(path, sample)
        if tissue is None or mate is None:
            continue
        key = (tissue, mate)
        old = chosen.get(key)
        if old is None:
            chosen[key] = path
        else:
            old_qc = "/results/qc/" in str(old)
            new_qc = "/results/qc/" in str(path)
            if new_qc and not old_qc:
                chosen[key] = path
            if read_text(old) != read_text(path):
                warnings.append(
                    f"Non-identical duplicate trimming reports for {tissue} R{mate}: "
                    f"{old} ; {path}"
                )

    rows = []
    for key, path in sorted(chosen.items()):
        parsed = parse_trim_report(path, sample)
        if parsed:
            rows.append(parsed)
    return rows, warnings


def infer_tissue_from_path(path: Path, sample: str, suffix_pattern: str) -> Optional[str]:
    match = re.match(
        rf"^{re.escape(sample)}[-_](.+?){suffix_pattern}$",
        path.name,
        flags=re.IGNORECASE,
    )
    if match:
        return match.group(1)
    parts = list(path.parts)
    for index, part in enumerate(parts):
        if part == sample and index > 0:
            previous = parts[index - 1]
            if previous.lower() not in {"align", "rnaseq", "results"}:
                return previous
    return None


def parse_star_log(path: Path, sample: str) -> Optional[Dict[str, object]]:
    tissue = infer_tissue_from_path(path, sample, r"\.STAR_Log\.final\.out")
    if tissue is None:
        tissue = infer_tissue_from_path(path, sample, r"\.Log\.final\.out")
    if tissue is None:
        return None
    values: Dict[str, str] = {}
    for line in read_text(path).splitlines():
        if "|" not in line:
            continue
        key, value = line.split("|", 1)
        values[key.strip()] = value.strip()

    def i(key: str) -> Optional[int]:
        return int_clean(values.get(key, ""))

    def f(key: str) -> Optional[float]:
        return float_clean(values.get(key, ""))

    return {
        "sample": sample,
        "tissue": tissue,
        "star_log": str(path),
        "input_reads": i("Number of input reads"),
        "average_input_read_length": f("Average input read length"),
        "uniquely_mapped_reads": i("Uniquely mapped reads number"),
        "uniquely_mapped_pct": f("Uniquely mapped reads %"),
        "multimapped_reads": i("Number of reads mapped to multiple loci"),
        "multimapped_pct": f("% of reads mapped to multiple loci"),
        "too_many_loci_reads": i("Number of reads mapped to too many loci"),
        "too_many_loci_pct": f("% of reads mapped to too many loci"),
        "unmapped_too_short_reads": i("Number of reads unmapped: too short"),
        "unmapped_too_short_pct": f("% of reads unmapped: too short"),
        "unmapped_other_reads": i("Number of reads unmapped: other"),
        "unmapped_other_pct": f("% of reads unmapped: other"),
        "mismatch_rate_pct": f("Mismatch rate per base, %"),
        "total_splices": i("Number of splices: Total"),
        "annotated_splices": i("Number of splices: Annotated (sjdb)"),
        "noncanonical_splices": i("Number of splices: Non-canonical"),
    }


def discover_star_logs(project: Path, sample: str) -> Tuple[Dict[str, Dict[str, object]], List[str]]:
    roots = [
        project / "results" / "qc" / "rnaseq" / "align",
        project / "results" / "rnaseq" / "align",
    ]
    candidates: List[Path] = []
    for root in roots:
        if root.is_dir():
            candidates.extend(root.rglob("*Log.final.out"))

    by_tissue: Dict[str, List[Tuple[int, Path, Dict[str, object]]]] = defaultdict(list)
    warnings: List[str] = []
    for path in sorted(set(candidates)):
        if "pass1" in {part.lower() for part in path.parts}:
            continue
        if sample not in str(path):
            continue
        parsed = parse_star_log(path, sample)
        if not parsed:
            continue
        score = 0
        if f"{sample}_" in path.name and ".STAR_Log.final.out" in path.name:
            score += 10
        if "/results/qc/" in str(path):
            score += 5
        by_tissue[str(parsed["tissue"])].append((score, path, parsed))

    selected: Dict[str, Dict[str, object]] = {}
    for tissue, items in by_tissue.items():
        items.sort(key=lambda item: (item[0], item[1].stat().st_mtime), reverse=True)
        selected[tissue] = items[0][2]
        if len(items) > 1:
            warnings.append(
                f"Multiple final STAR logs found for {tissue}; selected {items[0][1]}"
            )
    return selected, warnings


def parse_picard_metrics(path: Path, sample: str) -> Optional[Dict[str, object]]:
    lines = read_text(path).splitlines()
    header: Optional[List[str]] = None
    data: Optional[List[str]] = None
    in_metrics = False
    for line in lines:
        if line.startswith("## METRICS CLASS") and "DuplicationMetrics" in line:
            in_metrics = True
            continue
        if not in_metrics or not line.strip() or line.startswith("#"):
            continue
        if header is None:
            header = line.split("\t")
            continue
        data = line.split("\t")
        break
    if not header or not data:
        return None
    values = dict(zip(header, data))
    tissue = infer_tissue_from_path(path, sample, r"\.dedup\.metrics\.txt")
    if tissue is None:
        library = values.get("LIBRARY", "")
        match = re.match(rf"^{re.escape(sample)}[-_](.+)$", library)
        tissue = match.group(1) if match else None
    if tissue is None:
        return None

    return {
        "sample": sample,
        "tissue": tissue,
        "duplicate_metrics": str(path),
        "read_pairs_examined": int_clean(values.get("READ_PAIRS_EXAMINED", "")),
        "read_pair_duplicates": int_clean(values.get("READ_PAIR_DUPLICATES", "")),
        "optical_duplicate_pairs": int_clean(
            values.get("READ_PAIR_OPTICAL_DUPLICATES", "")
        ),
        "percent_duplication": float_clean(values.get("PERCENT_DUPLICATION", "")),
        "estimated_library_size": int_clean(values.get("ESTIMATED_LIBRARY_SIZE", "")),
    }


def discover_picard_metrics(
    project: Path, sample: str
) -> Tuple[Dict[str, Dict[str, object]], List[str]]:
    roots = [
        project / "results" / "qc" / "rnaseq" / "align",
        project / "results" / "rnaseq" / "align",
    ]
    candidates: List[Path] = []
    for root in roots:
        if root.is_dir():
            candidates.extend(root.rglob(f"*{sample}*.dedup.metrics.txt"))

    by_tissue: Dict[str, List[Tuple[int, Path, Dict[str, object]]]] = defaultdict(list)
    warnings: List[str] = []
    for path in sorted(set(candidates)):
        parsed = parse_picard_metrics(path, sample)
        if not parsed:
            continue
        score = 5 if "/results/qc/" in str(path) else 0
        by_tissue[str(parsed["tissue"])].append((score, path, parsed))

    selected: Dict[str, Dict[str, object]] = {}
    for tissue, items in by_tissue.items():
        items.sort(key=lambda item: (item[0], item[1].stat().st_mtime), reverse=True)
        selected[tissue] = items[0][2]
        if len(items) > 1:
            warnings.append(
                f"Multiple duplicate metrics found for {tissue}; selected {items[0][1]}"
            )
    return selected, warnings


def load_alignment_summary(project: Path, sample: str) -> Dict[str, Dict[str, object]]:
    candidates = list(
        (project / "results").rglob(f"{sample}.alignment_summary.tsv")
    ) if (project / "results").is_dir() else []
    if not candidates:
        return {}
    # Prefer the canonical QC path, otherwise newest.
    candidates.sort(
        key=lambda path: ("/results/qc/" in str(path), path.stat().st_mtime),
        reverse=True,
    )
    output: Dict[str, Dict[str, object]] = {}
    with candidates[0].open(encoding="utf-8", errors="replace", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            tissue = row.get("tissue", "").strip()
            if not tissue:
                continue
            output[tissue] = {
                "dedup_bam_reads": int_clean(row.get("dedup_bam_reads", "")),
                "phaseready_reads": int_clean(row.get("phaseready_reads", "")),
                "phaseready_pct": float_clean(row.get("phaseready_pct", "")),
            }
    return output


def discover_dedup_bams(project: Path, sample: str) -> Dict[str, Path]:
    root = project / "results" / "rnaseq" / "align"
    output: Dict[str, Path] = {}
    if not root.is_dir():
        return output
    for path in sorted(root.rglob(f"{sample}_*.dedup.bam")):
        tissue = infer_tissue_from_path(path, sample, r"\.dedup\.bam")
        if tissue:
            output[tissue] = path
    return output


def scan_vw_tags(
    bam: Path, samtools: str, threads: int
) -> Dict[str, object]:
    counts = defaultdict(int)
    command = [samtools, "view", "-@", str(max(1, threads)), str(bam)]
    process = subprocess.Popen(
        command,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        encoding="utf-8",
        errors="replace",
    )
    assert process.stdout is not None
    for line in process.stdout:
        fields = line.rstrip("\n").split("\t")
        if len(fields) < 11:
            continue
        counts["total"] += 1
        flag = int(fields[1])
        if not flag & 2304:  # neither secondary (256) nor supplementary (2048)
            counts["primary"] += 1
        tag_value = None
        for field in fields[11:]:
            if field.startswith("vW:i:"):
                try:
                    tag_value = int(field.rsplit(":", 1)[1])
                except ValueError:
                    tag_value = -1
                break
        if tag_value is None:
            counts["no_tag"] += 1
        elif tag_value == 1:
            counts["pass_1"] += 1
        elif 2 <= tag_value <= 7:
            counts["fail_2_7"] += 1
        else:
            counts["other"] += 1
    stderr = process.stderr.read() if process.stderr else ""
    return_code = process.wait()
    if return_code != 0:
        raise RuntimeError(
            f"samtools failed for {bam} (exit {return_code}): {stderr.strip()}"
        )
    return {
        "vw_bam": str(bam),
        "vw_total_alignments": counts["total"],
        "vw_primary_alignments": counts["primary"],
        "vw_pass_1": counts["pass_1"],
        "vw_fail_2_7": counts["fail_2_7"],
        "vw_no_tag": counts["no_tag"],
        "vw_other": counts["other"],
    }


def write_tsv(path: Path, fields: Sequence[str], rows: Iterable[Dict[str, object]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: numeric(row.get(field)) for field in fields})


def markdown_value(value: object, digits: int = 2) -> str:
    if value is None or value == "":
        return "NA"
    if isinstance(value, float):
        return f"{value:.{digits}f}"
    if isinstance(value, int):
        return f"{value:,}"
    return str(value)


def build_report(
    sample: str,
    trim_rows: List[Dict[str, object]],
    align_rows: List[Dict[str, object]],
    warnings: List[str],
    scan_vw: bool,
) -> str:
    tissues = sorted(
        {str(row["tissue"]) for row in trim_rows + align_rows if row.get("tissue")}
    )
    mates: Dict[str, set] = defaultdict(set)
    for row in trim_rows:
        mates[str(row["tissue"])].add(str(row["mate"]))

    lines = [
        f"# Stage 2 RNA-seq QC: {sample}",
        "",
        "## Discovery",
        "",
        f"- Tissues discovered: {len(tissues)} ({', '.join(tissues) if tissues else 'none'})",
        f"- Trimming reports: {len(trim_rows)}",
        f"- Final STAR logs: {sum(1 for row in align_rows if row.get('star_log'))}",
        f"- Duplicate metric files: {sum(1 for row in align_rows if row.get('duplicate_metrics'))}",
        f"- Exact WASP vW scan: {'performed' if scan_vw else 'not requested'}",
        "",
        "## Trimming summary",
        "",
        "| Tissue | R1 bases retained | R2 bases retained | R1 short removed | R2 short removed | R1 poly-G/C | R2 poly-G/C | overlap |",
        "|---|---:|---:|---:|---:|---:|---:|---:|",
    ]
    by_trim = {(str(row["tissue"]), str(row["mate"])): row for row in trim_rows}
    for tissue in sorted(mates):
        r1 = by_trim.get((tissue, "1"), {})
        r2 = by_trim.get((tissue, "2"), {})
        overlap = r1.get("adapter_overlap_bp") or r2.get("adapter_overlap_bp")
        lines.append(
            "| {t} | {r1b}% | {r2b}% | {r1s}% | {r2s}% | {r1p}% | {r2p}% | {o} |".format(
                t=tissue,
                r1b=markdown_value(r1.get("written_bases_pct")),
                r2b=markdown_value(r2.get("written_bases_pct")),
                r1s=markdown_value(r1.get("short_reads_removed_pct")),
                r2s=markdown_value(r2.get("short_reads_removed_pct")),
                r1p=markdown_value(r1.get("polyg_reads_pct")),
                r2p=markdown_value(r2.get("polyg_reads_pct")),
                o=markdown_value(overlap, 0),
            )
        )

    lines.extend(
        [
            "",
            "Adapter percentages reported by Trim Galore include very short matches. "
            "The detailed TSV separately reports 1–2 bp calls and calls of at least 3 or 5 bp.",
            "",
            "## Alignment and duplicate summary",
            "",
            "| Tissue | input reads | unique % | multimapped % | mismatch % | duplicate % | phase-ready % | WASP-fail vW2–7 |",
            "|---|---:|---:|---:|---:|---:|---:|---:|",
        ]
    )
    for row in sorted(align_rows, key=lambda item: str(item.get("tissue", ""))):
        duplicate = row.get("percent_duplication")
        duplicate_pct = 100.0 * duplicate if isinstance(duplicate, float) else None
        lines.append(
            "| {t} | {inp} | {uniq} | {multi} | {mm} | {dup} | {phase} | {vw} |".format(
                t=row.get("tissue", "NA"),
                inp=markdown_value(row.get("input_reads"), 0),
                uniq=markdown_value(row.get("uniquely_mapped_pct")),
                multi=markdown_value(row.get("multimapped_pct")),
                mm=markdown_value(row.get("mismatch_rate_pct")),
                dup=markdown_value(duplicate_pct),
                phase=markdown_value(row.get("phaseready_pct")),
                vw=markdown_value(row.get("vw_fail_2_7"), 0),
            )
        )

    observations: List[str] = []
    missing_mates = [t for t in sorted(mates) if mates[t] != {"1", "2"}]
    if missing_mates:
        observations.append("Missing R1/R2 mate report: " + ", ".join(missing_mates))
    low_base_retention = [
        f"{row['tissue']}-R{row['mate']}"
        for row in trim_rows
        if isinstance(row.get("written_bases_pct"), float)
        and float(row["written_bases_pct"]) < 95.0
    ]
    if low_base_retention:
        observations.append("Base retention below 95%: " + ", ".join(low_base_retention))
    overlap_one = [
        f"{row['tissue']}-R{row['mate']}"
        for row in trim_rows
        if row.get("adapter_overlap_bp") == 1
    ]
    if overlap_one:
        observations.append(
            "Adapter overlap/stringency is 1 bp; interpret the raw adapter-hit percentage "
            "together with the ≥3 bp and ≥5 bp columns."
        )
    low_unique = [
        str(row["tissue"])
        for row in align_rows
        if isinstance(row.get("uniquely_mapped_pct"), float)
        and float(row["uniquely_mapped_pct"]) < 70.0
    ]
    if low_unique:
        observations.append("Uniquely mapped fraction below 70%: " + ", ".join(low_unique))
    high_dup = [
        str(row["tissue"])
        for row in align_rows
        if isinstance(row.get("percent_duplication"), float)
        and float(row["percent_duplication"]) >= 0.50
    ]
    if high_dup:
        observations.append(
            "Duplicate fraction is at least 50% and requires ASE-specific review: "
            + ", ".join(high_dup)
            + ". phASER removes reads flagged as duplicate by default."
        )
    missing_star = [str(row["tissue"]) for row in align_rows if not row.get("star_log")]
    missing_dup = [
        str(row["tissue"]) for row in align_rows if not row.get("duplicate_metrics")
    ]
    if missing_star:
        observations.append("Missing final STAR log: " + ", ".join(missing_star))
    if missing_dup:
        observations.append("Missing duplicate metrics: " + ", ".join(missing_dup))

    lines.extend(["", "## Review indicators", ""])
    if observations:
        lines.extend([f"- {item}" for item in observations])
    else:
        lines.append("- No automatic review indicator was triggered.")
    if warnings:
        lines.extend(["", "## File-discovery warnings", ""])
        lines.extend([f"- {item}" for item in warnings])
    lines.extend(
        [
            "",
            "> Thresholds in this report are review indicators, not automatic biological "
            "pass/fail criteria. Final acceptance should consider all tissues and the ASE design.",
            "",
        ]
    )
    return "\n".join(lines)


def main() -> int:
    args = parse_args()
    project = args.project_dir.resolve()
    if not project.is_dir():
        print(f"[ERROR] Project directory not found: {project}", file=sys.stderr)
        return 2
    outdir = (
        args.outdir.resolve()
        if args.outdir
        else project / "results" / "qc" / "rnaseq" / "stage2_qc" / args.sample
    )
    outdir.mkdir(parents=True, exist_ok=True)

    trim_rows, trim_warnings = discover_trim_reports(project, args.sample)
    star_rows, star_warnings = discover_star_logs(project, args.sample)
    duplicate_rows, duplicate_warnings = discover_picard_metrics(project, args.sample)
    pipeline_summary = load_alignment_summary(project, args.sample)

    all_tissues = sorted(
        set(str(row["tissue"]) for row in trim_rows)
        | set(star_rows)
        | set(duplicate_rows)
        | set(pipeline_summary)
    )
    align_rows: List[Dict[str, object]] = []
    for tissue in all_tissues:
        row: Dict[str, object] = {"sample": args.sample, "tissue": tissue}
        row.update(star_rows.get(tissue, {}))
        row.update(duplicate_rows.get(tissue, {}))
        row.update(pipeline_summary.get(tissue, {}))
        align_rows.append(row)

    warnings = trim_warnings + star_warnings + duplicate_warnings
    if args.scan_vw:
        if shutil.which(args.samtools) is None:
            print(f"[ERROR] samtools not found: {args.samtools}", file=sys.stderr)
            return 2
        dedup_bams = discover_dedup_bams(project, args.sample)
        align_by_tissue = {str(row["tissue"]): row for row in align_rows}
        for tissue, bam in sorted(dedup_bams.items()):
            print(f"[vW] Scanning {tissue}: {bam}", file=sys.stderr, flush=True)
            result = scan_vw_tags(bam, args.samtools, args.threads)
            if tissue not in align_by_tissue:
                new_row: Dict[str, object] = {"sample": args.sample, "tissue": tissue}
                align_rows.append(new_row)
                align_by_tissue[tissue] = new_row
            align_by_tissue[tissue].update(result)
        missing_vw = sorted(set(all_tissues) - set(dedup_bams))
        if missing_vw:
            warnings.append("No dedup BAM found for vW scan: " + ", ".join(missing_vw))

    trim_path = outdir / f"{args.sample}.trim_qc_summary.tsv"
    align_path = outdir / f"{args.sample}.alignment_qc_summary.tsv"
    report_path = outdir / f"{args.sample}.stage2_qc_report.md"
    manifest_path = outdir / f"{args.sample}.stage2_qc_manifest.tsv"

    write_tsv(trim_path, TRIM_FIELDS, sorted(trim_rows, key=lambda r: (str(r["tissue"]), str(r["mate"]))))
    write_tsv(align_path, ALIGN_FIELDS, sorted(align_rows, key=lambda r: str(r["tissue"])))
    report_path.write_text(
        build_report(args.sample, trim_rows, align_rows, warnings, args.scan_vw),
        encoding="utf-8",
    )

    manifest_rows = []
    for row in trim_rows:
        manifest_rows.append({"type": "trim_report", "tissue": row["tissue"], "path": row["report"]})
    for row in align_rows:
        for field, file_type in [
            ("star_log", "star_log"),
            ("duplicate_metrics", "duplicate_metrics"),
            ("vw_bam", "vw_bam"),
        ]:
            if row.get(field):
                manifest_rows.append({"type": file_type, "tissue": row["tissue"], "path": row[field]})
    write_tsv(manifest_path, ["type", "tissue", "path"], manifest_rows)

    print("Stage 2 QC collection complete")
    print(f"  Trim reports:      {len(trim_rows)}")
    print(f"  Tissues:           {len(all_tissues)}")
    print(f"  Trim summary:      {trim_path}")
    print(f"  Alignment summary: {align_path}")
    print(f"  Markdown report:   {report_path}")
    print(f"  Manifest:          {manifest_path}")
    if not args.scan_vw:
        print("  vW scan:           skipped (use --scan-vw only if needed)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
