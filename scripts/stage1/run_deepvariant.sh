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
for tool in docker bcftools; do require_command "${tool}"; done

require_nonempty "${REFERENCE_FASTA}"
require_nonempty "${REFERENCE_FASTA}.fai"
case "${REFERENCE_FASTA}" in
    "${PROJECT_DIR}"/*) ;;
    *) die "REFERENCE_FASTA must be located below PROJECT_DIR for Docker mounting" ;;
esac

ALIGN_DIR="${PROJECT_DIR}/results/wgs/align"
OUT_DIR="${PROJECT_DIR}/results/wgs/variants"
QC_DIR="${PROJECT_DIR}/results/qc/wgs"
INTERMEDIATE_DIR="${OUT_DIR}/intermediate_${SAMPLE}"
DV_LOG_DIR="${OUT_DIR}/logs_${SAMPLE}"
BAM="${ALIGN_DIR}/${SAMPLE}.dedup.bam"
OUT_VCF="${OUT_DIR}/${SAMPLE}.deepvariant.vcf.gz"
OUT_GVCF="${OUT_DIR}/${SAMPLE}.deepvariant.g.vcf.gz"
STATS="${QC_DIR}/${SAMPLE}.deepvariant.bcftools.stats.txt"

require_nonempty "${BAM}"
bam_index "${BAM}" >/dev/null || die "BAM index missing for ${BAM}"
mkdir -p "${OUT_DIR}" "${QC_DIR}" "${INTERMEDIATE_DIR}" "${DV_LOG_DIR}"

gvcf_is_valid() {
    [[ -s "${OUT_GVCF}" && -s "${OUT_GVCF}.tbi" ]] || return 1
    bcftools view -h "${OUT_GVCF}" >/dev/null 2>&1
}

if vcf_is_valid "${OUT_VCF}" && gvcf_is_valid; then
    note "[SKIP] Complete DeepVariant VCF and gVCF outputs exist"
else
    note "Running ${DEEPVARIANT_IMAGE} for ${SAMPLE}"
    REF_IN_CONTAINER="/data/${REFERENCE_FASTA#${PROJECT_DIR}/}"
    docker run \
        --rm \
        -v "${PROJECT_DIR}:/data" \
        "${DEEPVARIANT_IMAGE}" \
        /opt/deepvariant/bin/run_deepvariant \
        --model_type=WGS \
        --ref="${REF_IN_CONTAINER}" \
        --reads="/data/results/wgs/align/${SAMPLE}.dedup.bam" \
        --output_vcf="/data/results/wgs/variants/${SAMPLE}.deepvariant.vcf.gz" \
        --output_gvcf="/data/results/wgs/variants/${SAMPLE}.deepvariant.g.vcf.gz" \
        --num_shards="${THREADS}" \
        --sample_name="${SAMPLE}" \
        --regions="${AUTOSOMES}" \
        --intermediate_results_dir="/data/results/wgs/variants/intermediate_${SAMPLE}" \
        --logging_dir="/data/results/wgs/variants/logs_${SAMPLE}" \
        --vcf_stats_report=true
fi

vcf_is_valid "${OUT_VCF}" || die "DeepVariant VCF failed validation: ${OUT_VCF}"
gvcf_is_valid || die "DeepVariant gVCF failed validation: ${OUT_GVCF}"

TMP_STATS="${STATS}.tmp.$$"
bcftools stats "${OUT_VCF}" > "${TMP_STATS}"
require_nonempty "${TMP_STATS}"
mv "${TMP_STATS}" "${STATS}"

note "DeepVariant complete: ${OUT_VCF}"
