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
require_nonempty "${REFERENCE_FASTA}"
require_nonempty "${REFERENCE_FASTA}.fai"

VCF_DIR="${PROJECT_DIR}/results/wgs/variants"
QC_DIR="${PROJECT_DIR}/results/qc/wgs"
VCF_IN="${VCF_DIR}/${SAMPLE}.deepvariant.vcf.gz"
VCF_OUT="${VCF_DIR}/${SAMPLE}.het.snp.vcf.gz"
QC_OUT="${QC_DIR}/${SAMPLE}.het_snp_qc.tsv"
DP_QC_OUT="${QC_DIR}/${SAMPLE}.het_dp_retention_qc.tsv"
mkdir -p "${VCF_DIR}" "${QC_DIR}"
vcf_is_valid "${VCF_IN}" || die "Input DeepVariant VCF is invalid: ${VCF_IN}"

if vcf_is_valid "${VCF_OUT}"; then
    note "[SKIP] Complete heterozygous-SNV VCF exists"
else
    AUTOSOME_CSV=${AUTOSOMES// /,}
    TMP_VCF="${VCF_DIR}/.${SAMPLE}.het.snp.$$.vcf.gz"
    trap 'rm -f "${TMP_VCF}" "${TMP_VCF}.tbi"' EXIT
    note "Preparing PASS biallelic heterozygous SNVs with DP>=${MIN_HET_DP}"

    # Multiallelic sites are excluded rather than split. This retains the
    # primary call-set definition used in the cohort analysis.
    bcftools view -r "${AUTOSOME_CSV}" -m2 -M2 "${VCF_IN}" -Ou \
    | bcftools norm -f "${REFERENCE_FASTA}" -c w -Ou \
    | bcftools view \
        -v snps \
        -f PASS \
        -g het \
        -i "FORMAT/DP>=${MIN_HET_DP}" \
        -Oz -o "${TMP_VCF}"

    bcftools index -t "${TMP_VCF}"
    vcf_is_valid "${TMP_VCF}" || die "Filtered het-SNV VCF failed validation"
    mv "${TMP_VCF}" "${VCF_OUT}"
    mv "${TMP_VCF}.tbi" "${VCF_OUT}.tbi"
    trap - EXIT
fi

python3 "${SCRIPT_DIR}/summarize_het_vcf.py" \
    --sample "${SAMPLE}" \
    --vcf "${VCF_OUT}" \
    --output "${QC_OUT}"

python3 "${SCRIPT_DIR}/summarize_dp_retention.py" \
    --sample "${SAMPLE}" \
    --deepvariant-vcf "${VCF_IN}" \
    --filtered-vcf "${VCF_OUT}" \
    --autosomes ${AUTOSOMES} \
    --minimum-dp "${MIN_HET_DP}" \
    --output "${DP_QC_OUT}"

note "Heterozygous-SNV preparation complete: ${VCF_OUT}"
