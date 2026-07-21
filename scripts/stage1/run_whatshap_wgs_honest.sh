#!/usr/bin/env bash

set -euo pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
# shellcheck source=common.sh
source "${SCRIPT_DIR}/common.sh"

usage() {
    cat <<'EOF'
Usage: bash run_whatshap_wgs_honest.sh --sample SAMPLE --config FILE [options]

Options:
  --threads N       Accepted for interface consistency; WhatsHap phase is not
                    parallelized by this script.
  --summary-only    Rebuild stats, block list, and summary from an existing
                    phased VCF. The phase-ready BAM is not required.
EOF
}

parse_stage1_args "$@"
if [[ "${SHOW_HELP}" -eq 1 ]]; then usage; exit 0; fi
load_stage1_config
make_autosome_files
for tool in whatshap bcftools python3; do require_command "${tool}"; done

VARIANT_DIR="${PROJECT_DIR}/results/wgs/variants"
ALIGN_DIR="${PROJECT_DIR}/results/wgs/align"
OUT_DIR="${PROJECT_DIR}/results/phasing"
mkdir -p "${OUT_DIR}"

HET_VCF="${VARIANT_DIR}/${SAMPLE}.het.snp.vcf.gz"
WGS_BAM="${ALIGN_DIR}/${SAMPLE}.dedup.phaseready.bam"
PREFIX="${OUT_DIR}/${SAMPLE}.phased.wgs.honest"
OUT_VCF="${PREFIX}.vcf.gz"
STATS_TSV="${PREFIX}.stats.tsv"
BLOCKS_TSV="${PREFIX}.blocks.tsv"
BLOCKS_GTF="${PREFIX}.blocks.gtf"
READS_TXT="${PREFIX}.used_reads.txt"
RUN_LOG="${PREFIX}.log"
SUMMARY_TSV="${PREFIX}.summary.tsv"

if ! vcf_is_valid "${OUT_VCF}"; then
    [[ "${SUMMARY_ONLY}" -eq 0 ]] \
        || die "--summary-only requires a valid phased VCF: ${OUT_VCF}"
    vcf_is_valid "${HET_VCF}" || die "Invalid het-SNV VCF: ${HET_VCF}"
    require_nonempty "${WGS_BAM}"
    bam_index "${WGS_BAM}" >/dev/null || die "BAM index missing for ${WGS_BAM}"

    TMP_VCF="${OUT_DIR}/.${SAMPLE}.phased.wgs.honest.$$.vcf.gz"
    TMP_READS="${OUT_DIR}/.${SAMPLE}.phased.wgs.honest.used_reads.$$.txt"
    trap 'rm -f "${TMP_VCF}" "${TMP_VCF}.tbi" "${TMP_READS}"' EXIT
    note "Running WhatsHap local read-backed phasing for ${SAMPLE}"
    whatshap phase \
        -o "${TMP_VCF}" \
        --reference "${REFERENCE_FASTA}" \
        --tag PS \
        --only-snvs \
        --sample "${SAMPLE}" \
        --mapping-quality "${PHASING_MIN_MAPQ}" \
        --output-read-list "${TMP_READS}" \
        "${HET_VCF}" "${WGS_BAM}" \
        2>&1 | tee "${RUN_LOG}"
    bcftools index -t "${TMP_VCF}"
    vcf_is_valid "${TMP_VCF}" || die "WhatsHap phased VCF failed validation"
    require_nonempty "${TMP_READS}"
    mv "${TMP_VCF}" "${OUT_VCF}"
    mv "${TMP_VCF}.tbi" "${OUT_VCF}.tbi"
    mv "${TMP_READS}" "${READS_TXT}"
    trap - EXIT
else
    note "[SKIP] Valid phased VCF exists; checking downstream QC artifacts"
fi

if [[ ! -s "${STATS_TSV}" || ! -s "${BLOCKS_TSV}" || ! -s "${BLOCKS_GTF}" ]]; then
    TMP_STATS="${STATS_TSV}.tmp.$$"
    TMP_BLOCKS="${BLOCKS_TSV}.tmp.$$"
    TMP_GTF="${BLOCKS_GTF}.tmp.$$"
    trap 'rm -f "${TMP_STATS}" "${TMP_BLOCKS}" "${TMP_GTF}"' EXIT
    note "Generating WhatsHap stats and block list"
    whatshap stats \
        --tsv "${TMP_STATS}" \
        --block-list "${TMP_BLOCKS}" \
        --gtf "${TMP_GTF}" \
        --chr-lengths "${AUTOSOME_LENGTHS}" \
        --only-snvs \
        "${OUT_VCF}"
    require_nonempty "${TMP_STATS}"
    require_nonempty "${TMP_BLOCKS}"
    require_nonempty "${TMP_GTF}"
    mv "${TMP_STATS}" "${STATS_TSV}"
    mv "${TMP_BLOCKS}" "${BLOCKS_TSV}"
    mv "${TMP_GTF}" "${BLOCKS_GTF}"
    trap - EXIT
fi

python3 "${SCRIPT_DIR}/summarize_whatshap.py" \
    --stats "${STATS_TSV}" \
    --blocks "${BLOCKS_TSV}" \
    --output "${SUMMARY_TSV}" \
    --autosomes ${AUTOSOMES}
require_nonempty "${SUMMARY_TSV}"

note "WhatsHap phasing/QC complete: ${OUT_VCF}"
