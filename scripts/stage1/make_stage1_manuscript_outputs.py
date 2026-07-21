#!/usr/bin/env python3

"""Create concise Stage 1 supplementary figures and reporting tables."""

from __future__ import annotations

import argparse
import csv
import shutil
import statistics
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch


REQUIRED_FIELDS = [
    "sample",
    "mean_autosome_depth",
    "autosome_fraction_ge_10x",
    "primary_mapped_fraction",
    "properly_paired_fraction",
    "duplicate_fraction",
    "heterozygous_snv_before_dp10",
    "final_heterozygous_snv_count",
    "dp10_retention_fraction",
    "phased_snv_count",
    "phased_fraction",
    "multi_variant_blocks",
    "block_span_n50_bp",
    "variants_per_block_n50",
]

DICTIONARY = [
    ("sample", "Sample identifier", "text"),
    ("mean_autosome_depth", "Length-weighted mean depth across autosomes 1-18", "x"),
    ("autosome_fraction_ge_10x", "Fraction of autosomal bases covered at >=10x", "fraction"),
    ("primary_mapped_fraction", "Primary mapped reads reported by samtools flagstat", "fraction"),
    ("properly_paired_fraction", "Properly paired reads reported by samtools flagstat", "fraction"),
    ("duplicate_fraction", "Duplicate fraction reported by GATK MarkDuplicates", "fraction"),
    ("heterozygous_snv_before_dp10", "PASS biallelic heterozygous SNVs before DP filtering", "count"),
    ("final_heterozygous_snv_count", "Primary heterozygous SNVs retained at DP>=10", "count"),
    ("dp10_retention_fraction", "Fraction of pre-DP heterozygous SNVs retained at DP>=10", "fraction"),
    ("phased_snv_count", "Heterozygous SNVs represented with phased genotypes", "count"),
    ("phased_fraction", "Phased SNVs divided by input heterozygous SNVs", "fraction"),
    ("multi_variant_blocks", "WhatsHap phase blocks containing multiple variants", "count"),
    ("block_span_n50_bp", "N50 of observed phase-block spans", "bp"),
    ("variants_per_block_n50", "N50 of variant counts per phase block", "variants"),
]

COLORS = {
    "navy": "#173F5F",
    "blue": "#20639B",
    "teal": "#2A9D8F",
    "gold": "#E9C46A",
    "orange": "#F4A261",
    "red": "#E76F51",
    "gray": "#59636E",
    "light": "#F4F7F9",
}


def read_cohort(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = reader.fieldnames or []
        missing = [field for field in REQUIRED_FIELDS if field not in fields]
        if missing:
            raise ValueError(f"Missing cohort-QC fields: {', '.join(missing)}")
        rows = list(reader)
    if not rows:
        raise ValueError("Cohort-QC table contains no samples")
    return rows


def box(ax: plt.Axes, x: float, y: float, width: float, height: float, title: str,
        body: str, color: str) -> None:
    patch = FancyBboxPatch(
        (x, y), width, height,
        boxstyle="round,pad=0.015,rounding_size=0.025",
        linewidth=1.5, edgecolor=color, facecolor="white",
    )
    ax.add_patch(patch)
    ax.text(x + width / 2, y + height * 0.72, title, ha="center", va="center",
            fontsize=11, fontweight="bold", color=color)
    ax.text(x + width / 2, y + height * 0.37, body, ha="center", va="center",
            fontsize=8.6, color="#27323A", linespacing=1.35)


def arrow(ax: plt.Axes, start: tuple[float, float], end: tuple[float, float],
          color: str = "#68737D") -> None:
    ax.add_patch(FancyArrowPatch(start, end, arrowstyle="-|>", mutation_scale=14,
                                linewidth=1.4, color=color))


def workflow_figure(output_dir: Path) -> None:
    fig, ax = plt.subplots(figsize=(12.5, 7.0))
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    fig.patch.set_facecolor("white")

    ax.text(0.04, 0.95, "Stage 1: WGS preprocessing and local read-backed phasing",
            fontsize=17, fontweight="bold", color=COLORS["navy"], va="top")
    ax.text(0.04, 0.905, "Inputs, primary analysis outputs, and manuscript-level QC",
            fontsize=9.5, color=COLORS["gray"], va="top")

    w, h = 0.245, 0.225
    coords = {
        "input": (0.045, 0.585),
        "align": (0.375, 0.585),
        "call": (0.705, 0.585),
        "filter": (0.705, 0.205),
        "phase": (0.375, 0.205),
        "deliver": (0.045, 0.205),
    }
    box(ax, *coords["input"], w, h, "INPUT",
        "Paired-end WGS FASTQ\nSscrofa11.1 reference", COLORS["blue"])
    box(ax, *coords["align"], w, h, "READ PROCESSING + ALIGNMENT",
        "fastp -> BWA-MEM -> MarkDuplicates\nOutput: deduplicated BAM", COLORS["teal"])
    box(ax, *coords["call"], w, h, "VARIANT CALLING",
        "DeepVariant 1.10.0 (WGS)\nAutosomes 1-18\nOutput: raw VCF", COLORS["orange"])
    box(ax, *coords["filter"], w, h, "PRIMARY HET-SNV SET",
        "PASS; biallelic SNV; heterozygous\nFORMAT/DP >=10\nOutput: het-SNV VCF", COLORS["red"])
    box(ax, *coords["phase"], w, h, "LOCAL PHASING",
        "WhatsHap 2.8; PS tags\nProper-pair, MAPQ >=30 BAM\nOutput: phased VCF + phase QC", COLORS["blue"])
    box(ax, *coords["deliver"], w, h, "STAGE 1 DELIVERABLES",
        "Phased genotype backbone for ASE/ASM\nCohort QC table and figures", COLORS["navy"])

    arrow(ax, (0.29, 0.697), (0.375, 0.697))
    arrow(ax, (0.62, 0.697), (0.705, 0.697))
    arrow(ax, (0.827, 0.585), (0.827, 0.43))
    arrow(ax, (0.705, 0.317), (0.62, 0.317))
    arrow(ax, (0.375, 0.317), (0.29, 0.317))

    qc_x, qc_y, qc_w, qc_h = 0.375, 0.495, 0.245, 0.055
    qc = FancyBboxPatch((qc_x, qc_y), qc_w, qc_h,
                        boxstyle="round,pad=0.01,rounding_size=0.015",
                        linewidth=1, edgecolor=COLORS["gray"],
                        facecolor=COLORS["light"])
    ax.add_patch(qc)
    ax.text(qc_x + qc_w / 2, qc_y + qc_h / 2,
            "mosdepth + flagstat + duplicate metrics -> alignment QC TSV",
            ha="center", va="center", fontsize=7.8, color=COLORS["gray"])
    arrow(ax, (0.497, 0.585), (0.497, 0.55), COLORS["gray"])

    ax.text(0.045, 0.095,
            "Interpretation: WhatsHap phase orientation is valid within each PS block; "
            "Stage 1 does not provide chromosome-scale or parental haplotypes.",
            fontsize=8.5, color=COLORS["gray"], va="center")
    fig.tight_layout(pad=0.8)
    save_figure(fig, output_dir / "Supplementary_Figure_S1_Stage1_workflow")
    plt.close(fig)


def save_figure(fig: plt.Figure, stem: Path) -> None:
    fig.savefig(stem.with_suffix(".pdf"), bbox_inches="tight")
    fig.savefig(stem.with_suffix(".svg"), bbox_inches="tight")
    fig.savefig(stem.with_suffix(".png"), dpi=300, bbox_inches="tight")


def cohort_qc_figure(rows: list[dict[str, str]], output_dir: Path) -> None:
    samples = [row["sample"] for row in rows]
    x = list(range(len(rows)))
    depth = [float(row["mean_autosome_depth"]) for row in rows]
    coverage = [100 * float(row["autosome_fraction_ge_10x"]) for row in rows]
    het_million = [float(row["final_heterozygous_snv_count"]) / 1e6 for row in rows]
    phased = [100 * float(row["phased_fraction"]) for row in rows]
    block_n50_kb = [float(row["block_span_n50_bp"]) / 1000 for row in rows]

    width = max(10.5, 0.38 * len(rows) + 5.2)
    fig, axes = plt.subplots(2, 2, figsize=(width, 8.2), sharex=True)
    fig.patch.set_facecolor("white")
    fig.suptitle("Stage 1 cohort quality control", x=0.06, y=0.995,
                 ha="left", fontsize=16, fontweight="bold", color=COLORS["navy"])

    ax = axes[0, 0]
    ax.bar(x, depth, color=COLORS["blue"], alpha=0.86, label="Mean depth")
    ax.set_ylabel("Mean autosomal depth (x)")
    ax2 = ax.twinx()
    ax2.plot(x, coverage, "o", color=COLORS["orange"], markersize=4.5,
             label="Bases >=10x")
    ax2.set_ylabel("Autosomal bases >=10x (%)", color=COLORS["orange"])
    ax2.tick_params(axis="y", colors=COLORS["orange"])
    ax.set_title("A  Sequencing depth and breadth", loc="left", fontweight="bold")

    ax = axes[0, 1]
    ax.vlines(x, 0, het_million, color=COLORS["teal"], linewidth=2)
    ax.scatter(x, het_million, color=COLORS["teal"], s=28, zorder=3)
    ax.set_ylabel("Final heterozygous SNVs (million)")
    ax.set_ylim(bottom=0)
    ax.set_title("B  Primary heterozygous-SNV set", loc="left", fontweight="bold")

    ax = axes[1, 0]
    ax.scatter(x, phased, color=COLORS["red"], s=32, zorder=3)
    ax.plot(x, phased, color=COLORS["red"], alpha=0.35, linewidth=1)
    low, high = min(phased), max(phased)
    margin = max(0.5, (high - low) * 0.25)
    ax.set_ylim(max(0, low - margin), min(100.5, high + margin))
    ax.set_ylabel("Phased heterozygous SNVs (%)")
    ax.set_title("C  Local phasing rate", loc="left", fontweight="bold")

    ax = axes[1, 1]
    ax.bar(x, block_n50_kb, color=COLORS["gold"], edgecolor="#B58627")
    ax.set_ylabel("Phase-block span N50 (kb)")
    ax.set_ylim(bottom=0)
    ax.set_title("D  Local phase-block continuity", loc="left", fontweight="bold")

    for ax in axes.flat:
        ax.grid(axis="y", color="#D9E0E5", linewidth=0.7, alpha=0.8)
        ax.set_axisbelow(True)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.set_xticks(x)
        ax.set_xticklabels(samples, rotation=60, ha="right", fontsize=8)
    axes[0, 0].tick_params(labelbottom=False)
    axes[0, 1].tick_params(labelbottom=False)

    fig.tight_layout(rect=(0.03, 0.04, 0.99, 0.96), h_pad=2.0, w_pad=2.3)
    save_figure(fig, output_dir / "Supplementary_Figure_S2_Stage1_cohort_QC")
    plt.close(fig)


def write_dictionary(output: Path) -> None:
    with output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(("column", "definition", "unit"))
        writer.writerows(DICTIONARY)


def write_statistics(rows: list[dict[str, str]], output: Path) -> None:
    metrics = [field for field in REQUIRED_FIELDS if field != "sample"]
    with output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(("metric", "n", "median", "minimum", "maximum"))
        for metric in metrics:
            values = [float(row[metric]) for row in rows]
            writer.writerow((metric, len(values), statistics.median(values), min(values), max(values)))


def write_legends(output: Path) -> None:
    output.write_text(
        "# Stage 1 supplementary legends\n\n"
        "## Supplementary Figure S1\n\n"
        "Stage 1 workflow for autosomal WGS preprocessing, germline variant calling, "
        "heterozygous-SNV filtering, and local read-backed phasing. Paired-end WGS reads "
        "were processed with fastp, aligned to Sscrofa11.1 with BWA-MEM, and duplicate "
        "records were marked with GATK MarkDuplicates. Autosomal variants were called with "
        "DeepVariant 1.10.0. PASS biallelic heterozygous SNVs with FORMAT/DP >=10 were "
        "retained and phased with WhatsHap 2.8 using proper-pair reads with MAPQ >=30. "
        "Phase orientation is local to each PS block.\n\n"
        "## Supplementary Figure S2\n\n"
        "Stage 1 cohort quality control. (A) Length-weighted mean autosomal depth and the "
        "percentage of autosomal bases covered at >=10x. (B) Number of PASS biallelic "
        "heterozygous SNVs retained after the DP >=10 filter. (C) Percentage of input "
        "heterozygous SNVs phased by WhatsHap. (D) N50 of observed local phase-block spans.\n\n"
        "## Supplementary Table S1\n\n"
        "Per-sample WGS alignment, primary heterozygous-SNV, and local phasing QC metrics. "
        "Fractions are stored on a 0-1 scale.\n",
        encoding="utf-8",
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--cohort-qc", type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--workflow-only", action="store_true")
    args = parser.parse_args()
    if not args.workflow_only and args.cohort_qc is None:
        parser.error("--cohort-qc is required unless --workflow-only is used")

    args.output_dir.mkdir(parents=True, exist_ok=True)
    workflow_figure(args.output_dir)
    write_legends(args.output_dir / "stage1_supplementary_legends.md")
    if args.workflow_only:
        return

    assert args.cohort_qc is not None
    rows = read_cohort(args.cohort_qc)
    cohort_qc_figure(rows, args.output_dir)
    shutil.copyfile(
        args.cohort_qc,
        args.output_dir / "Supplementary_Table_S1_Stage1_cohort_QC.tsv",
    )
    write_dictionary(args.output_dir / "Supplementary_Table_S1_metric_dictionary.tsv")
    write_statistics(rows, args.output_dir / "stage1_cohort_descriptive_statistics.tsv")


if __name__ == "__main__":
    main()
