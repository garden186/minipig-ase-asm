#!/usr/bin/env bash

set -euo pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
# shellcheck source=common.sh
source "${SCRIPT_DIR}/common.sh"

usage() {
    echo "Usage: bash $0 --sample SAMPLE --config FILE [--threads N]"
}

parse_stage1_args "$@"
if [[ "${SHOW_HELP}" -eq 1 ]]; then usage; exit 0; fi
load_stage1_config
for tool in bcftools python3; do require_command "${tool}"; done

TRIM_DIR="${PROJECT_DIR}/results/wgs/trim/${SAMPLE}"
ALIGN_DIR="${PROJECT_DIR}/results/wgs/align"
VARIANT_DIR="${PROJECT_DIR}/results/wgs/variants"
QC_DIR="${PROJECT_DIR}/results/qc/wgs"
PHASING_DIR="${PROJECT_DIR}/results/phasing"

DEEPVARIANT_VCF="${VARIANT_DIR}/${SAMPLE}.deepvariant.vcf.gz"
HET_VCF="${VARIANT_DIR}/${SAMPLE}.het.snp.vcf.gz"
PHASED_PREFIX="${PHASING_DIR}/${SAMPLE}.phased.wgs.honest"
PHASED_VCF="${PHASED_PREFIX}.vcf.gz"

for path in \
    "${TRIM_DIR}/${SAMPLE}.fastp.json" \
    "${ALIGN_DIR}/${SAMPLE}.dedup.bam" \
    "${ALIGN_DIR}/${SAMPLE}.dedup.metrics.txt" \
    "${QC_DIR}/${SAMPLE}.mosdepth.mosdepth.summary.txt" \
    "${QC_DIR}/${SAMPLE}.flagstat.txt" \
    "${QC_DIR}/${SAMPLE}.stage1_alignment_qc.tsv" \
    "${QC_DIR}/${SAMPLE}.het_snp_qc.tsv" \
    "${QC_DIR}/${SAMPLE}.het_dp_retention_qc.tsv" \
    "${PHASED_PREFIX}.stats.tsv" \
    "${PHASED_PREFIX}.blocks.tsv" \
    "${PHASED_PREFIX}.summary.tsv"; do
    require_nonempty "${path}"
done

vcf_is_valid "${DEEPVARIANT_VCF}" || die "Invalid DeepVariant VCF"
vcf_is_valid "${HET_VCF}" || die "Invalid heterozygous-SNV VCF"
vcf_is_valid "${PHASED_VCF}" || die "Invalid phased VCF"

HET_COUNT=$(bcftools index -n "${HET_VCF}")
PHASED_COUNT=$(bcftools index -n "${PHASED_VCF}")
[[ "${HET_COUNT}" -eq "${PHASED_COUNT}" ]] \
    || die "Het and phased VCF counts differ: ${HET_COUNT} vs ${PHASED_COUNT}"

python3 - \
    "${PHASED_PREFIX}.stats.tsv" \
    "${PHASED_PREFIX}.summary.tsv" \
    "${HET_COUNT}" <<'PY'
import csv
import sys

stats_path, summary_path, expected_text = sys.argv[1:]
expected = int(expected_text)

with open(stats_path, encoding="utf-8", newline="") as handle:
    stats = list(csv.DictReader(handle, delimiter="\t"))
all_stats = next((row for row in stats if row["chromosome"] == "ALL"), None)
if all_stats is None or int(all_stats["variants"]) != expected:
    raise SystemExit("WhatsHap ALL stats do not match the phased VCF count")

with open(summary_path, encoding="utf-8", newline="") as handle:
    summary = list(csv.DictReader(handle, delimiter="\t"))
all_summary = next((row for row in summary if row["chrom"] == "ALL"), None)
if all_summary is None or int(all_summary["variants"]) != expected:
    raise SystemExit("WhatsHap ALL summary does not match the phased VCF count")
if len(summary) != 19:
    raise SystemExit(f"Expected 18 autosomes plus ALL, observed {len(summary)} rows")
PY

if [[ ! -s "${ALIGN_DIR}/${SAMPLE}.dedup.phaseready.bam" ]]; then
    echo "[WARN] Phase-ready BAM is not on this filesystem; restore its archived copy only to rerun phasing."
fi
if [[ ! -s "${PHASED_PREFIX}.used_reads.txt" ]]; then
    echo "[WARN] WhatsHap used-read list is absent; phased VCF and summary validation still passed."
fi

echo "[PASS] Stage 1 validated for ${SAMPLE}"
echo "  DeepVariant records: $(bcftools index -n "${DEEPVARIANT_VCF}")"
echo "  Het/phased records:  ${HET_COUNT}"
awk -F'\t' 'NR==1 || $1=="ALL"' "${PHASED_PREFIX}.summary.tsv"
