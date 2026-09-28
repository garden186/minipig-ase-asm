#!/usr/bin/env python3
"""Reconstruct SNP-level allelic imbalance from ASE fragment assignments.

The existing v0.2.6 candidate archive contains block-level Haplotype-A/B
counts but not per-SNP allele counts.  This audit reads the full v0.2.6 block
table and the server-only fragment-assignment audit, reconstructs REF/ALT
fragment counts for each sample-tissue-SNP, applies an exact two-sided
binomial test, and performs BH correction within each sample-tissue family.

This is a site-level audit within the existing ASE block-analysis universe;
it is not a new genome-wide ASE caller.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import math
import re
import sys
from collections import defaultdict
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, Iterator, Mapping, MutableMapping, Sequence, Tuple


VERSION = "ase_snp_level_v0.1.0"
DEFAULT_SAMPLES = "0235,0242,0276,0309,0326,0349,0355,0377,0384,0385"
TISSUE_NORMALIZATION = {
    "Tenderlo": "Tenderloin",
    "blood": "Blood",
    "L-blood": "Blood",
}
VARIANT_RE = re.compile(r"^([^_]+)_(\d+)_([ACGTN]+)_([ACGTN]+)$", re.I)

BLOCK_REQUIRED = {
    "sample",
    "tissue",
    "gene_id",
    "gene_name",
    "haplotype_block_id",
    "measurement_id",
    "analysis_set",
    "block_snps",
    "haplotypeA",
    "haplotypeB",
    "gene_variants_with_reads",
    "v026_technical_qc_status",
    "v026_is_core_ase",
    "v026_ase_direction",
}
FRAGMENT_REQUIRED = {
    "sample",
    "tissue",
    "haplotype_block_id",
    "qname",
    "gene_id",
    "haplotype",
    "final_status",
    "gene_supported_variants",
}

SITE_HEADER = [
    "sample",
    "tissue",
    "variant_id",
    "chrom",
    "position",
    "ref",
    "alt",
    "ref_fragment_count",
    "alt_fragment_count",
    "total_informative_fragments",
    "minor_fragment_count",
    "major_allele_fraction",
    "major_allele",
    "site_p_exact",
    "site_q_bh",
    "site_bh_scope",
    "site_testable",
    "site_fdr_significant",
    "site_strong_effect",
    "site_level_ase_signal",
    "has_core_parent_block",
    "core_parent_direction_status",
    "core_concordant_site_level_ase_signal",
    "n_all_parent_blocks",
    "all_parent_gene_ids",
    "all_parent_gene_names",
    "all_parent_block_ids",
    "n_core_parent_blocks",
    "core_parent_gene_ids",
    "core_parent_gene_names",
    "core_parent_block_ids",
    "core_parent_measurement_ids",
    "core_parent_major_alleles",
    "fragment_source",
]

GENE_SITE_HEADER = [
    "gene_id",
    "gene_name",
    "tissue",
    "variant_id",
    "chrom",
    "position",
    "ref",
    "alt",
    "n_samples_rna_covered_in_core_block",
    "samples_rna_covered_in_core_block",
    "n_samples_site_testable",
    "samples_site_testable",
    "n_samples_site_level_ase_signal",
    "samples_site_level_ase_signal",
    "n_samples_core_direction_concordant_site_ase",
    "samples_core_direction_concordant_site_ase",
]


class AnalysisError(RuntimeError):
    pass


def normalize_sample(value: str) -> str:
    text = str(value).strip()
    return text.zfill(4) if text.isdigit() else text


def normalize_tissue(value: str) -> str:
    text = str(value).strip()
    return TISSUE_NORMALIZATION.get(text, text)


def split_values(value: str) -> list[str]:
    return [item.strip() for item in str(value or "").split(",") if item.strip()]


def parse_variant(value: str) -> tuple[str, int, str, str] | None:
    match = VARIANT_RE.match(str(value).strip())
    if not match:
        return None
    chrom, position, ref, alt = match.groups()
    ref, alt = ref.upper(), alt.upper()
    if len(ref) != 1 or len(alt) != 1 or ref not in "ACGT" or alt not in "ACGT":
        return None
    return chrom.removeprefix("chr"), int(position), ref, alt


def open_text(path: Path, mode: str = "rt"):
    if path.suffix == ".gz":
        return gzip.open(path, mode, newline="")
    return path.open(mode.replace("t", ""), encoding="utf-8", newline="") if "b" not in mode else path.open(mode)


def iter_tsv(path: Path, required: set[str]) -> Iterator[dict[str, str]]:
    with open_text(path, "rt") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise AnalysisError(f"{path}: missing required columns: {sorted(missing)}")
        yield from reader


def find_existing(root: Path, names: Sequence[str]) -> Path:
    for name in names:
        path = root / name
        if path.is_file() and path.stat().st_size > 0:
            return path
    raise AnalysisError(f"none of the required files exists below {root}: {list(names)}")


def exact_binom_two_sided_half(ref_count: int, alt_count: int) -> float | None:
    total = ref_count + alt_count
    if total <= 0:
        return None
    smaller = min(ref_count, alt_count)
    base = math.lgamma(total + 1) - total * math.log(2.0)
    log_probabilities = [
        base - math.lgamma(index + 1) - math.lgamma(total - index + 1)
        for index in range(smaller + 1)
    ]
    maximum = max(log_probabilities)
    log_cdf = maximum + math.log(sum(math.exp(value - maximum) for value in log_probabilities))
    return max(sys.float_info.min, min(1.0, 2.0 * math.exp(log_cdf)))


def bh_adjust(indexed_pvalues: Sequence[Tuple[int, float]]) -> Dict[int, float]:
    ranked = sorted(indexed_pvalues, key=lambda item: (item[1], item[0]))
    result: Dict[int, float] = {}
    running = 1.0
    total = len(ranked)
    for rank_index in range(total - 1, -1, -1):
        row_index, pvalue = ranked[rank_index]
        rank = rank_index + 1
        running = min(running, pvalue * total / rank)
        result[row_index] = min(1.0, running)
    return result


def fmt_float(value: float | None) -> str:
    return "" if value is None or not math.isfinite(value) else f"{value:.12g}"


@dataclass(frozen=True)
class Parent:
    gene_id: str
    gene_name: str
    block_id: str
    measurement_id: str
    is_core: bool
    major_allele: str


@dataclass
class BlockInfo:
    sample: str
    tissue: str
    gene_id: str
    gene_name: str
    block_id: str
    measurement_id: str
    is_core: bool
    allele_a: dict[str, str]
    allele_b: dict[str, str]
    parent_major_allele: dict[str, str]


@dataclass
class SiteAccumulator:
    ref_count: int = 0
    alt_count: int = 0

    def add(self, allele: str, ref: str, alt: str) -> None:
        if allele == ref:
            self.ref_count += 1
        elif allele == alt:
            self.alt_count += 1


def block_path(block_root: Path, sample: str) -> Path:
    return find_existing(
        block_root / sample,
        ("ase_haplotype_block_reclassified.tsv.gz", "ase_haplotype_block_reclassified.tsv"),
    )


def fragment_path(fragment_root: Path, v025_root: Path, sample: str) -> Path:
    names = ("audit/ase_fragment_assignment.tsv.gz", "audit/ase_fragment_assignment.tsv")
    for root in (fragment_root / sample, v025_root / sample):
        try:
            return find_existing(root, names)
        except AnalysisError:
            continue
    raise AnalysisError(
        f"sample {sample}: fragment assignment missing in both {fragment_root / sample} "
        f"and {v025_root / sample}"
    )


def load_blocks(path: Path, sample: str) -> tuple[dict[tuple[str, str, str], BlockInfo], dict[tuple[str, str], set[Parent]], dict[str, int]]:
    blocks: dict[tuple[str, str, str], BlockInfo] = {}
    parents: MutableMapping[tuple[str, str], set[Parent]] = defaultdict(set)
    qc = defaultdict(int)
    for row in iter_tsv(path, BLOCK_REQUIRED):
        qc["n_block_rows"] += 1
        if normalize_sample(row["sample"]) != sample:
            raise AnalysisError(f"{path}: row sample {row['sample']} does not match directory sample {sample}")
        if row["analysis_set"] != "MAIN":
            qc["n_non_main_block_rows"] += 1
            continue
        if row["v026_technical_qc_status"] != "PASS":
            qc["n_technical_fail_block_rows"] += 1
            continue

        block_snps = split_values(row["block_snps"])
        hap_a = split_values(row["haplotypeA"])
        hap_b = split_values(row["haplotypeB"])
        if not (len(block_snps) == len(hap_a) == len(hap_b)):
            raise AnalysisError(
                f"{path}: unaligned block SNP/haplotype columns in {row['measurement_id']}"
            )
        aligned = dict(zip(block_snps, zip(hap_a, hap_b)))
        observed = []
        for variant in split_values(row["gene_variants_with_reads"]):
            parsed = parse_variant(variant)
            if parsed is None:
                qc["n_non_snp_or_malformed_variants"] += 1
                continue
            if variant not in aligned:
                qc["n_read_variants_missing_from_block"] += 1
                continue
            allele_a, allele_b = (str(x).upper() for x in aligned[variant])
            _, _, ref, alt = parsed
            if allele_a not in {ref, alt} or allele_b not in {ref, alt} or allele_a == allele_b:
                qc["n_invalid_block_allele_variants"] += 1
                continue
            observed.append((variant, allele_a, allele_b))
        if not observed:
            qc["n_blocks_without_usable_snp"] += 1
            continue

        tissue = normalize_tissue(row["tissue"])
        gene_id = row["gene_id"].split(".")[0]
        key = (tissue, gene_id, row["haplotype_block_id"])
        if key in blocks:
            raise AnalysisError(f"{path}: duplicate block key: {key}")
        is_core = str(row["v026_is_core_ase"]).strip().upper() in {"1", "TRUE", "YES"}
        direction = row["v026_ase_direction"]
        allele_a_map = {variant: a for variant, a, _ in observed}
        allele_b_map = {variant: b for variant, _, b in observed}
        if direction == "HAP_A":
            major = dict(allele_a_map)
        elif direction == "HAP_B":
            major = dict(allele_b_map)
        else:
            major = {}
        info = BlockInfo(
            sample=sample,
            tissue=tissue,
            gene_id=gene_id,
            gene_name=row["gene_name"],
            block_id=row["haplotype_block_id"],
            measurement_id=row["measurement_id"],
            is_core=is_core,
            allele_a=allele_a_map,
            allele_b=allele_b_map,
            parent_major_allele=major,
        )
        blocks[key] = info
        for variant, _, _ in observed:
            parents[(tissue, variant)].add(
                Parent(
                    gene_id=gene_id,
                    gene_name=row["gene_name"],
                    block_id=info.block_id,
                    measurement_id=info.measurement_id,
                    is_core=is_core,
                    major_allele=major.get(variant, ""),
                )
            )
            qc["n_block_variant_links"] += 1
            if is_core:
                qc["n_core_block_variant_links"] += 1
        qc["n_loaded_blocks"] += 1
        if is_core:
            qc["n_loaded_core_blocks"] += 1
    return blocks, dict(parents), dict(qc)


def count_sites(
    path: Path,
    sample: str,
    blocks: Mapping[tuple[str, str, str], BlockInfo],
) -> tuple[dict[tuple[str, str], SiteAccumulator], dict[str, int]]:
    sites: MutableMapping[tuple[str, str], SiteAccumulator] = defaultdict(SiteAccumulator)
    qc = defaultdict(int)
    for row in iter_tsv(path, FRAGMENT_REQUIRED):
        qc["n_fragment_rows"] += 1
        if row["final_status"] != "ASSIGNED" or row["haplotype"] not in {"A", "B"}:
            qc["n_unassigned_or_invalid_haplotype_rows"] += 1
            continue
        tissue = normalize_tissue(row["tissue"])
        gene_id = row["gene_id"].split(".")[0]
        block_key = (tissue, gene_id, row["haplotype_block_id"])
        block = blocks.get(block_key)
        if block is None:
            qc["n_fragment_rows_outside_selected_block_universe"] += 1
            continue
        qc["n_target_assigned_fragment_rows"] += 1
        variants = set(split_values(row["gene_supported_variants"]))
        if not variants:
            qc["n_target_fragments_without_supported_variant"] += 1
            continue
        allele_map = block.allele_a if row["haplotype"] == "A" else block.allele_b
        for variant in variants:
            allele = allele_map.get(variant)
            parsed = parse_variant(variant)
            if allele is None or parsed is None:
                qc["n_fragment_variant_links_not_usable"] += 1
                continue
            _, _, ref, alt = parsed
            sites[(tissue, variant)].add(allele, ref, alt)
            qc["n_usable_fragment_variant_links"] += 1
    return dict(sites), dict(qc)


def parent_direction_status(core_parents: Sequence[Parent], major_allele: str) -> str:
    if not core_parents:
        return "NOT_APPLICABLE"
    alleles = {parent.major_allele for parent in core_parents if parent.major_allele}
    if not major_allele:
        return "SITE_TIE"
    if len(alleles) != 1:
        return "PARENT_DIRECTION_CONFLICT"
    return "CONSISTENT" if major_allele in alleles else "OPPOSITE"


def make_site_rows(
    sample: str,
    fragment_source: Path,
    sites: Mapping[tuple[str, str], SiteAccumulator],
    parents_by_site: Mapping[tuple[str, str], set[Parent]],
    min_total: int,
    min_major_fraction: float,
    fdr: float,
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for tissue, variant in sorted(sites, key=lambda item: (item[0], parse_variant(item[1]) or ("", 0, "", ""))):
        parsed = parse_variant(variant)
        if parsed is None:
            continue
        chrom, position, ref, alt = parsed
        accumulator = sites[(tissue, variant)]
        ref_count = accumulator.ref_count
        alt_count = accumulator.alt_count
        total = ref_count + alt_count
        major_fraction = max(ref_count, alt_count) / total if total else None
        major_allele = ref if ref_count > alt_count else alt if alt_count > ref_count else ""
        pvalue = exact_binom_two_sided_half(ref_count, alt_count) if total >= min_total else None
        parents = sorted(
            parents_by_site.get((tissue, variant), set()),
            key=lambda parent: (parent.gene_id, parent.block_id, parent.measurement_id),
        )
        core_parents = [parent for parent in parents if parent.is_core]
        status = parent_direction_status(core_parents, major_allele)
        rows.append(
            {
                "sample": sample,
                "tissue": tissue,
                "variant_id": variant,
                "chrom": chrom,
                "position": position,
                "ref": ref,
                "alt": alt,
                "ref_fragment_count": ref_count,
                "alt_fragment_count": alt_count,
                "total_informative_fragments": total,
                "minor_fragment_count": min(ref_count, alt_count),
                "major_allele_fraction": major_fraction,
                "major_allele": major_allele,
                "site_p_exact": pvalue,
                "site_q_bh": None,
                "site_bh_scope": f"{sample}|{tissue}",
                "site_testable": pvalue is not None,
                "site_fdr_significant": False,
                "site_strong_effect": bool(major_fraction is not None and major_fraction >= min_major_fraction),
                "site_level_ase_signal": False,
                "has_core_parent_block": bool(core_parents),
                "core_parent_direction_status": status,
                "core_concordant_site_level_ase_signal": False,
                "n_all_parent_blocks": len({parent.block_id for parent in parents}),
                "all_parent_gene_ids": ",".join(sorted({parent.gene_id for parent in parents})),
                "all_parent_gene_names": ",".join(sorted({parent.gene_name for parent in parents})),
                "all_parent_block_ids": ",".join(sorted({parent.block_id for parent in parents})),
                "n_core_parent_blocks": len({parent.block_id for parent in core_parents}),
                "core_parent_gene_ids": ",".join(sorted({parent.gene_id for parent in core_parents})),
                "core_parent_gene_names": ",".join(sorted({parent.gene_name for parent in core_parents})),
                "core_parent_block_ids": ",".join(sorted({parent.block_id for parent in core_parents})),
                "core_parent_measurement_ids": ",".join(sorted({parent.measurement_id for parent in core_parents})),
                "core_parent_major_alleles": ",".join(sorted({parent.major_allele for parent in core_parents if parent.major_allele})),
                "fragment_source": str(fragment_source),
                "_core_parents": core_parents,
            }
        )

    families: MutableMapping[str, list[tuple[int, float]]] = defaultdict(list)
    for index, row in enumerate(rows):
        if row["site_p_exact"] is not None:
            families[str(row["tissue"])].append((index, float(row["site_p_exact"])))
    for indexed in families.values():
        for index, qvalue in bh_adjust(indexed).items():
            row = rows[index]
            row["site_q_bh"] = qvalue
            row["site_fdr_significant"] = qvalue <= fdr
            row["site_level_ase_signal"] = bool(row["site_fdr_significant"] and row["site_strong_effect"])
            row["core_concordant_site_level_ase_signal"] = bool(
                row["site_level_ase_signal"] and row["core_parent_direction_status"] == "CONSISTENT"
            )
    return rows


def serializable_site_row(row: Mapping[str, object]) -> dict[str, object]:
    result = {key: row.get(key, "") for key in SITE_HEADER}
    for field_name in ("major_allele_fraction", "site_p_exact", "site_q_bh"):
        value = result[field_name]
        result[field_name] = fmt_float(float(value)) if value is not None and value != "" else ""
    for field_name in (
        "site_testable",
        "site_fdr_significant",
        "site_strong_effect",
        "site_level_ase_signal",
        "has_core_parent_block",
        "core_concordant_site_level_ase_signal",
    ):
        result[field_name] = int(bool(result[field_name]))
    return result


def add_gene_site_aggregate(
    aggregate: MutableMapping[tuple[str, str, str, str], dict[str, object]],
    row: Mapping[str, object],
) -> None:
    core_parents: Sequence[Parent] = row.get("_core_parents", [])  # type: ignore[assignment]
    by_gene: MutableMapping[tuple[str, str], set[str]] = defaultdict(set)
    for parent in core_parents:
        if parent.major_allele:
            by_gene[(parent.gene_id, parent.gene_name)].add(parent.major_allele)
        else:
            by_gene[(parent.gene_id, parent.gene_name)]
    sample = str(row["sample"])
    for (gene_id, gene_name), parent_major_alleles in by_gene.items():
        key = (gene_id, gene_name, str(row["tissue"]), str(row["variant_id"]))
        target = aggregate.setdefault(
            key,
            {
                "rna": set(),
                "testable": set(),
                "signal": set(),
                "concordant": set(),
                "chrom": row["chrom"],
                "position": row["position"],
                "ref": row["ref"],
                "alt": row["alt"],
            },
        )
        target["rna"].add(sample)  # type: ignore[union-attr]
        if row["site_testable"]:
            target["testable"].add(sample)  # type: ignore[union-attr]
        if row["site_level_ase_signal"]:
            target["signal"].add(sample)  # type: ignore[union-attr]
        if (
            row["site_level_ase_signal"]
            and len(parent_major_alleles) == 1
            and row["major_allele"] in parent_major_alleles
        ):
            target["concordant"].add(sample)  # type: ignore[union-attr]


def write_gene_site_summary(
    path: Path,
    aggregate: Mapping[tuple[str, str, str, str], Mapping[str, object]],
) -> None:
    with gzip.open(path, "wt", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=GENE_SITE_HEADER, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for (gene_id, gene_name, tissue, variant), value in sorted(
            aggregate.items(), key=lambda item: (item[0][2], item[0][0], item[1]["chrom"], item[1]["position"])
        ):
            rna = sorted(value["rna"])  # type: ignore[arg-type]
            testable = sorted(value["testable"])  # type: ignore[arg-type]
            signal = sorted(value["signal"])  # type: ignore[arg-type]
            concordant = sorted(value["concordant"])  # type: ignore[arg-type]
            writer.writerow(
                {
                    "gene_id": gene_id,
                    "gene_name": gene_name,
                    "tissue": tissue,
                    "variant_id": variant,
                    "chrom": value["chrom"],
                    "position": value["position"],
                    "ref": value["ref"],
                    "alt": value["alt"],
                    "n_samples_rna_covered_in_core_block": len(rna),
                    "samples_rna_covered_in_core_block": ",".join(rna),
                    "n_samples_site_testable": len(testable),
                    "samples_site_testable": ",".join(testable),
                    "n_samples_site_level_ase_signal": len(signal),
                    "samples_site_level_ase_signal": ",".join(signal),
                    "n_samples_core_direction_concordant_site_ase": len(concordant),
                    "samples_core_direction_concordant_site_ase": ",".join(concordant),
                }
            )


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--block-root", required=True, type=Path)
    parser.add_argument("--fragment-root", required=True, type=Path)
    parser.add_argument("--v025-root", required=True, type=Path)
    parser.add_argument("--out-dir", required=True, type=Path)
    parser.add_argument("--samples", default=DEFAULT_SAMPLES)
    parser.add_argument("--min-total", type=int, default=15)
    parser.add_argument("--min-major-fraction", type=float, default=0.65)
    parser.add_argument("--fdr", type=float, default=0.05)
    return parser.parse_args(argv)


def validate_args(args: argparse.Namespace) -> list[str]:
    samples = [normalize_sample(item) for item in args.samples.split(",") if item.strip()]
    if not samples or len(samples) != len(set(samples)):
        raise AnalysisError("--samples must contain unique sample IDs")
    if args.min_total < 1:
        raise AnalysisError("--min-total must be >=1")
    if not 0.5 <= args.min_major_fraction <= 1.0:
        raise AnalysisError("--min-major-fraction must be between 0.5 and 1")
    if not 0.0 < args.fdr < 1.0:
        raise AnalysisError("--fdr must be between 0 and 1")
    for root in (args.block_root, args.fragment_root, args.v025_root):
        if not root.is_dir():
            raise AnalysisError(f"missing input directory: {root}")
    args.out_dir.mkdir(parents=True, exist_ok=True)
    expected = (
        "site_level_ase_measurements.tsv.gz",
        "core_block_site_level_ase_measurements.tsv.gz",
        "core_gene_variant_tissue_summary.tsv.gz",
        "snp_level_ase_summary.tsv",
        "run_metadata.json",
    )
    existing = [str(args.out_dir / name) for name in expected if (args.out_dir / name).exists()]
    if existing:
        raise AnalysisError(f"refusing to overwrite existing outputs: {existing}")
    return samples


def main(argv: Sequence[str] | None = None) -> None:
    args = parse_args(argv)
    samples = validate_args(args)
    all_path = args.out_dir / "site_level_ase_measurements.tsv.gz"
    core_path = args.out_dir / "core_block_site_level_ase_measurements.tsv.gz"
    summary_path = args.out_dir / "snp_level_ase_summary.tsv"
    gene_summary_path = args.out_dir / "core_gene_variant_tissue_summary.tsv.gz"

    summary_rows: list[dict[str, object]] = []
    sample_qc: dict[str, object] = {}
    gene_aggregate: MutableMapping[tuple[str, str, str, str], dict[str, object]] = {}

    with gzip.open(all_path, "wt", newline="") as all_handle, gzip.open(core_path, "wt", newline="") as core_handle:
        all_writer = csv.DictWriter(all_handle, fieldnames=SITE_HEADER, delimiter="\t", lineterminator="\n")
        core_writer = csv.DictWriter(core_handle, fieldnames=SITE_HEADER, delimiter="\t", lineterminator="\n")
        all_writer.writeheader()
        core_writer.writeheader()

        for sample in samples:
            blocks_file = block_path(args.block_root, sample)
            fragments_file = fragment_path(args.fragment_root, args.v025_root, sample)
            print(f"[{sample}] loading blocks: {blocks_file}", flush=True)
            blocks, parents, block_qc = load_blocks(blocks_file, sample)
            print(f"[{sample}] streaming fragments: {fragments_file}", flush=True)
            sites, fragment_qc = count_sites(fragments_file, sample, blocks)
            rows = make_site_rows(
                sample,
                fragments_file,
                sites,
                parents,
                args.min_total,
                args.min_major_fraction,
                args.fdr,
            )
            per_tissue: MutableMapping[str, dict[str, int]] = defaultdict(lambda: defaultdict(int))
            for row in rows:
                serialized = serializable_site_row(row)
                all_writer.writerow(serialized)
                tissue_counts = per_tissue[str(row["tissue"])]
                tissue_counts["n_rna_covered_sites"] += 1
                if row["site_testable"]:
                    tissue_counts["n_site_testable"] += 1
                if row["site_level_ase_signal"]:
                    tissue_counts["n_site_level_ase_signal"] += 1
                if row["has_core_parent_block"]:
                    core_writer.writerow(serialized)
                    tissue_counts["n_core_block_rna_covered_sites"] += 1
                    if row["site_testable"]:
                        tissue_counts["n_core_block_site_testable"] += 1
                    if row["site_level_ase_signal"]:
                        tissue_counts["n_core_block_site_level_ase_signal"] += 1
                    if row["core_concordant_site_level_ase_signal"]:
                        tissue_counts["n_core_concordant_site_level_ase_signal"] += 1
                    add_gene_site_aggregate(gene_aggregate, row)
            for tissue, counts in sorted(per_tissue.items()):
                summary_rows.append({"sample": sample, "tissue": tissue, **counts})
            sample_qc[sample] = {
                "block_file": str(blocks_file),
                "fragment_file": str(fragments_file),
                "block_qc": block_qc,
                "fragment_qc": fragment_qc,
                "n_output_site_rows": len(rows),
            }

    summary_header = [
        "sample",
        "tissue",
        "n_rna_covered_sites",
        "n_site_testable",
        "n_site_level_ase_signal",
        "n_core_block_rna_covered_sites",
        "n_core_block_site_testable",
        "n_core_block_site_level_ase_signal",
        "n_core_concordant_site_level_ase_signal",
    ]
    with summary_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=summary_header, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in summary_rows:
            writer.writerow({key: row.get(key, 0) for key in summary_header})

    write_gene_site_summary(gene_summary_path, gene_aggregate)

    metadata = {
        "version": VERSION,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "analysis_scope": (
            "SNP-level exact-binomial audit within MAIN, technical-PASS rows of the existing "
            "v0.2.6 ASE block-analysis universe"
        ),
        "unit": "sample-tissue-SNP",
        "counts": "distinct assigned RNA fragment qnames supporting REF or ALT at each SNP",
        "multiple_testing": "BH across site-testable SNPs within each sample-tissue",
        "site_signal_definition": {
            "minimum_total_informative_fragments": args.min_total,
            "maximum_q_bh": args.fdr,
            "minimum_major_allele_fraction": args.min_major_fraction,
        },
        "core_concordant_definition": (
            "site-level ASE signal whose major allele matches the major allele implied by its core parent block"
        ),
        "samples": samples,
        "inputs": {
            "block_root": str(args.block_root),
            "fragment_root": str(args.fragment_root),
            "v025_fallback_root": str(args.v025_root),
        },
        "outputs": {
            "all_site_measurements": str(all_path),
            "core_block_site_measurements": str(core_path),
            "core_gene_variant_tissue_summary": str(gene_summary_path),
            "sample_tissue_summary": str(summary_path),
        },
        "sample_qc": sample_qc,
        "status": "PASS",
    }
    (args.out_dir / "run_metadata.json").write_text(
        json.dumps(metadata, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(
        json.dumps(
            {
                "version": VERSION,
                "status": "PASS",
                "n_samples": len(samples),
                "n_core_gene_variant_tissue_rows": len(gene_aggregate),
                "out_dir": str(args.out_dir),
            },
            sort_keys=True,
        ),
        flush=True,
    )


if __name__ == "__main__":
    try:
        main()
    except AnalysisError as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        raise SystemExit(2)
