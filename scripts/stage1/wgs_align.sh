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
load_sample_fastqs
make_autosome_files

for tool in bwa samtools gatk mosdepth python3; do require_command "${tool}"; done
require_nonempty "${REFERENCE_FASTA}"
require_nonempty "${REFERENCE_FASTA}.fai"
for suffix in amb ann bwt pac sa; do require_nonempty "${REFERENCE_FASTA}.${suffix}"; done
REFERENCE_DICT="${REFERENCE_FASTA%.*}.dict"
require_nonempty "${REFERENCE_DICT}"

TRIM_DIR="${PROJECT_DIR}/results/wgs/trim/${SAMPLE}"
R1="${TRIM_DIR}/${SAMPLE}_R1.fastq.gz"
R2="${TRIM_DIR}/${SAMPLE}_R2.fastq.gz"
require_nonempty "${R1}"
require_nonempty "${R2}"

ALIGN_DIR="${PROJECT_DIR}/results/wgs/align"
QC_DIR="${PROJECT_DIR}/results/qc/wgs"
TMP_DIR="${PROJECT_DIR}/tmp/${SAMPLE}"
mkdir -p "${ALIGN_DIR}" "${QC_DIR}" "${TMP_DIR}"

SORTED_BAM="${ALIGN_DIR}/${SAMPLE}.sorted.bam"
DEDUP_BAM="${ALIGN_DIR}/${SAMPLE}.dedup.bam"
DEDUP_METRICS="${ALIGN_DIR}/${SAMPLE}.dedup.metrics.txt"
PHASEREADY_BAM="${ALIGN_DIR}/${SAMPLE}.dedup.phaseready.bam"
MOSDEPTH_PREFIX="${QC_DIR}/${SAMPLE}.mosdepth"
FLAGSTAT="${QC_DIR}/${SAMPLE}.flagstat.txt"
ALIGN_QC="${QC_DIR}/${SAMPLE}.stage1_alignment_qc.tsv"

if [[ -n "${OPTICAL_DUPLICATE_PIXEL_DISTANCE:-}" ]]; then
    PIXEL_DISTANCE=${OPTICAL_DUPLICATE_PIXEL_DISTANCE}
else
    FIRST_READ=$(gzip -cd "${R1}" 2>/dev/null | awk 'NR==1 {print; exit}' || true)
    case "${FIRST_READ}" in
        @D*|@M*) PIXEL_DISTANCE=100 ;;
        *)       PIXEL_DISTANCE=2500 ;;
    esac
fi
[[ "${PIXEL_DISTANCE}" =~ ^[0-9]+$ ]] || die "Invalid optical duplicate pixel distance"

if [[ -s "${SORTED_BAM}" ]] && bam_index "${SORTED_BAM}" >/dev/null && samtools quickcheck "${SORTED_BAM}"; then
    note "[SKIP] Sorted BAM exists"
else
    TMP_SORTED="${TMP_DIR}/${SAMPLE}.sorted.$$.bam"
    trap 'rm -f "${TMP_SORTED}" "${TMP_SORTED}.bai"' EXIT
    note "Aligning ${SAMPLE} with BWA-MEM"
    bwa mem \
        -t "${THREADS}" \
        -K 100000000 \
        -Y \
        -R "@RG\tID:${SAMPLE}\tSM:${SAMPLE}\tPL:ILLUMINA\tLB:${SAMPLE}" \
        "${REFERENCE_FASTA}" "${R1}" "${R2}" \
    | samtools sort \
        -@ "${SAMTOOLS_SORT_THREADS}" \
        -m "${SAMTOOLS_SORT_MEMORY}" \
        -T "${TMP_DIR}/sort_${SAMPLE}" \
        -o "${TMP_SORTED}" -
    samtools quickcheck "${TMP_SORTED}"
    samtools index -@ "${SAMTOOLS_SORT_THREADS}" "${TMP_SORTED}"
    mv "${TMP_SORTED}" "${SORTED_BAM}"
    mv "${TMP_SORTED}.bai" "${SORTED_BAM}.bai"
    trap - EXIT
fi

if [[ -s "${DEDUP_BAM}" && -s "${DEDUP_METRICS}" ]] \
    && bam_index "${DEDUP_BAM}" >/dev/null && samtools quickcheck "${DEDUP_BAM}"; then
    note "[SKIP] Duplicate-marked BAM exists"
else
    TMP_DEDUP="${TMP_DIR}/${SAMPLE}.dedup.$$.bam"
    TMP_METRICS="${TMP_DIR}/${SAMPLE}.dedup.metrics.$$.txt"
    trap 'rm -f "${TMP_DEDUP}" "${TMP_DEDUP}.bai" "${TMP_METRICS}"' EXIT
    note "Marking duplicates for ${SAMPLE}"
    gatk --java-options "-Djava.io.tmpdir=${TMP_DIR} ${GATK_JAVA_OPTIONS}" \
        MarkDuplicates \
        -I "${SORTED_BAM}" \
        -O "${TMP_DEDUP}" \
        -M "${TMP_METRICS}" \
        --TMP_DIR "${TMP_DIR}" \
        --REMOVE_DUPLICATES false \
        --VALIDATION_STRINGENCY LENIENT \
        --CREATE_INDEX false \
        --ASSUME_SORT_ORDER coordinate \
        --OPTICAL_DUPLICATE_PIXEL_DISTANCE "${PIXEL_DISTANCE}" \
        --TAGGING_POLICY OpticalOnly
    samtools quickcheck "${TMP_DEDUP}"
    samtools index -@ "${SAMTOOLS_SORT_THREADS}" "${TMP_DEDUP}"
    require_nonempty "${TMP_METRICS}"
    mv "${TMP_DEDUP}" "${DEDUP_BAM}"
    mv "${TMP_DEDUP}.bai" "${DEDUP_BAM}.bai"
    mv "${TMP_METRICS}" "${DEDUP_METRICS}"
    trap - EXIT
fi

if [[ -s "${MOSDEPTH_PREFIX}.mosdepth.summary.txt" \
      && -s "${MOSDEPTH_PREFIX}.mosdepth.global.dist.txt" \
      && -s "${MOSDEPTH_PREFIX}.regions.bed.gz" ]]; then
    note "[SKIP] mosdepth outputs exist"
else
    note "Calculating alignment depth with mosdepth"
    mosdepth \
        --threads "${MOSDEPTH_THREADS}" \
        --by 1000 \
        --no-per-base \
        --fast-mode \
        --mapq "${MOSDEPTH_MIN_MAPQ}" \
        --flag 3844 \
        "${MOSDEPTH_PREFIX}" "${DEDUP_BAM}"
fi

if [[ ! -s "${FLAGSTAT}" ]]; then
    TMP_FLAGSTAT="${FLAGSTAT}.tmp.$$"
    samtools flagstat -@ "${SAMTOOLS_SORT_THREADS}" "${DEDUP_BAM}" > "${TMP_FLAGSTAT}"
    require_nonempty "${TMP_FLAGSTAT}"
    mv "${TMP_FLAGSTAT}" "${FLAGSTAT}"
fi

if [[ -s "${PHASEREADY_BAM}" ]] && bam_index "${PHASEREADY_BAM}" >/dev/null \
    && samtools quickcheck "${PHASEREADY_BAM}"; then
    note "[SKIP] Phase-ready BAM exists"
else
    TMP_PHASE="${TMP_DIR}/${SAMPLE}.dedup.phaseready.$$.bam"
    trap 'rm -f "${TMP_PHASE}" "${TMP_PHASE}.bai"' EXIT
    note "Creating proper-pair MAPQ>=${PHASING_MIN_MAPQ} phase-ready BAM"
    samtools view \
        -@ "${THREADS}" \
        -b \
        -f 2 \
        -F 3844 \
        -q "${PHASING_MIN_MAPQ}" \
        -o "${TMP_PHASE}" "${DEDUP_BAM}"
    samtools quickcheck "${TMP_PHASE}"
    samtools index -@ "${SAMTOOLS_SORT_THREADS}" "${TMP_PHASE}"
    mv "${TMP_PHASE}" "${PHASEREADY_BAM}"
    mv "${TMP_PHASE}.bai" "${PHASEREADY_BAM}.bai"
    trap - EXIT
fi

python3 "${SCRIPT_DIR}/summarize_alignment_qc.py" \
    --sample "${SAMPLE}" \
    --mosdepth-summary "${MOSDEPTH_PREFIX}.mosdepth.summary.txt" \
    --mosdepth-global-dist "${MOSDEPTH_PREFIX}.mosdepth.global.dist.txt" \
    --autosome-lengths "${AUTOSOME_LENGTHS}" \
    --flagstat "${FLAGSTAT}" \
    --duplication-metrics "${DEDUP_METRICS}" \
    --output "${ALIGN_QC}"

note "Alignment and QC complete: ${ALIGN_QC}"
