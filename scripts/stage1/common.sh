#!/usr/bin/env bash

set -euo pipefail

STAGE1_SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
STAGE1_REPO_ROOT=$(cd "${STAGE1_SCRIPT_DIR}/../.." && pwd)

die() {
    echo "[ERROR] $*" >&2
    exit 1
}

note() {
    echo "[$(date '+%F %T')] $*"
}

require_command() {
    command -v "$1" >/dev/null 2>&1 || die "Required command not found: $1"
}

require_nonempty() {
    [[ -s "$1" ]] || die "Required non-empty file not found: $1"
}

resolve_from_project() {
    local path=$1
    if [[ "${path}" = /* ]]; then
        printf '%s\n' "${path}"
    else
        printf '%s/%s\n' "${PROJECT_DIR%/}" "${path}"
    fi
}

bam_index() {
    local bam=$1
    if [[ -s "${bam}.bai" ]]; then
        printf '%s\n' "${bam}.bai"
    elif [[ -s "${bam%.bam}.bai" ]]; then
        printf '%s\n' "${bam%.bam}.bai"
    else
        return 1
    fi
}

vcf_is_valid() {
    local vcf=$1
    [[ -s "${vcf}" && -s "${vcf}.tbi" ]] || return 1
    bcftools view -h "${vcf}" >/dev/null 2>&1 || return 1
    [[ $(bcftools index -n "${vcf}" 2>/dev/null || echo 0) -gt 0 ]]
}

parse_stage1_args() {
    SAMPLE=""
    CONFIG_FILE="${STAGE1_REPO_ROOT}/config/stage1.env"
    THREADS=16
    SUMMARY_ONLY=0

    while [[ $# -gt 0 ]]; do
        case "$1" in
            --sample)
                [[ $# -ge 2 ]] || die "--sample requires a value"
                SAMPLE=$2
                shift 2
                ;;
            --config)
                [[ $# -ge 2 ]] || die "--config requires a value"
                CONFIG_FILE=$2
                shift 2
                ;;
            --threads)
                [[ $# -ge 2 ]] || die "--threads requires a value"
                THREADS=$2
                shift 2
                ;;
            --summary-only)
                SUMMARY_ONLY=1
                shift
                ;;
            -h|--help)
                SHOW_HELP=1
                shift
                ;;
            *)
                die "Unknown argument: $1"
                ;;
        esac
    done

    SHOW_HELP=${SHOW_HELP:-0}
    [[ "${SHOW_HELP}" -eq 1 ]] && return 0
    [[ -n "${SAMPLE}" ]] || die "--sample is required"
    [[ "${SAMPLE}" =~ ^[A-Za-z0-9._-]+$ ]] || die "Unsafe sample identifier: ${SAMPLE}"
    [[ "${THREADS}" =~ ^[1-9][0-9]*$ ]] || die "--threads must be a positive integer"
    require_nonempty "${CONFIG_FILE}"
}

load_stage1_config() {
    # shellcheck disable=SC1090
    source "${CONFIG_FILE}"

    : "${PROJECT_DIR:?PROJECT_DIR must be set in the Stage 1 config}"
    : "${REFERENCE_FASTA:?REFERENCE_FASTA must be set in the Stage 1 config}"
    : "${SAMPLE_SHEET:?SAMPLE_SHEET must be set in the Stage 1 config}"
    : "${AUTOSOMES:?AUTOSOMES must be set in the Stage 1 config}"

    MIN_HET_DP=${MIN_HET_DP:-10}
    MOSDEPTH_MIN_MAPQ=${MOSDEPTH_MIN_MAPQ:-20}
    PHASING_MIN_MAPQ=${PHASING_MIN_MAPQ:-30}
    DEEPVARIANT_IMAGE=${DEEPVARIANT_IMAGE:-google/deepvariant:1.10.0}
    FASTP_THREADS=${FASTP_THREADS:-8}
    MOSDEPTH_THREADS=${MOSDEPTH_THREADS:-4}
    SAMTOOLS_SORT_THREADS=${SAMTOOLS_SORT_THREADS:-8}
    SAMTOOLS_SORT_MEMORY=${SAMTOOLS_SORT_MEMORY:-4G}
    GATK_JAVA_OPTIONS=${GATK_JAVA_OPTIONS:-"-Xmx16g -XX:ParallelGCThreads=4"}

    [[ "${MIN_HET_DP}" =~ ^[0-9]+$ ]] || die "MIN_HET_DP must be an integer"
    [[ "${MOSDEPTH_MIN_MAPQ}" =~ ^[0-9]+$ ]] || die "MOSDEPTH_MIN_MAPQ must be an integer"
    [[ "${PHASING_MIN_MAPQ}" =~ ^[0-9]+$ ]] || die "PHASING_MIN_MAPQ must be an integer"

    mkdir -p "${PROJECT_DIR}"
    PROJECT_DIR=$(cd "${PROJECT_DIR}" && pwd)
    REFERENCE_FASTA=$(resolve_from_project "${REFERENCE_FASTA}")
    SAMPLE_SHEET=$(resolve_from_project "${SAMPLE_SHEET}")

    require_nonempty "${SAMPLE_SHEET}"
}

load_sample_fastqs() {
    local row
    row=$(awk -F'\t' -v s="${SAMPLE}" '
        $0 !~ /^#/ && $1 == s { print $2 "\t" $3; found=1; exit }
        END { if (!found) exit 1 }
    ' "${SAMPLE_SHEET}") || die "Sample ${SAMPLE} not found in ${SAMPLE_SHEET}"

    IFS=$'\t' read -r WGS_R1 WGS_R2 <<< "${row}"
    [[ -n "${WGS_R1}" && -n "${WGS_R2}" ]] || die "Missing FASTQ path for ${SAMPLE}"
    WGS_R1=$(resolve_from_project "${WGS_R1}")
    WGS_R2=$(resolve_from_project "${WGS_R2}")
    require_nonempty "${WGS_R1}"
    require_nonempty "${WGS_R2}"
}

make_autosome_files() {
    require_nonempty "${REFERENCE_FASTA}.fai"
    REF_META_DIR="${PROJECT_DIR}/results/reference"
    AUTOSOME_LENGTHS="${REF_META_DIR}/autosome.lengths.tsv"
    AUTOSOME_BED="${REF_META_DIR}/autosomes.bed"
    mkdir -p "${REF_META_DIR}"

    local tmp_lengths="${AUTOSOME_LENGTHS}.tmp.$$"
    local tmp_bed="${AUTOSOME_BED}.tmp.$$"
    : > "${tmp_lengths}"
    : > "${tmp_bed}"
    local chrom length
    for chrom in ${AUTOSOMES}; do
        length=$(awk -F'\t' -v c="${chrom}" '$1==c {print $2; exit}' "${REFERENCE_FASTA}.fai")
        [[ -n "${length}" ]] || die "Autosome ${chrom} is absent from ${REFERENCE_FASTA}.fai"
        printf '%s\t%s\n' "${chrom}" "${length}" >> "${tmp_lengths}"
        printf '%s\t0\t%s\n' "${chrom}" "${length}" >> "${tmp_bed}"
    done
    mv "${tmp_lengths}" "${AUTOSOME_LENGTHS}"
    mv "${tmp_bed}" "${AUTOSOME_BED}"
}

write_software_versions() {
    local output=$1
    {
        printf 'captured_at\t%s\n' "$(date --iso-8601=seconds 2>/dev/null || date)"
        printf 'deepvariant_image\t%s\n' "${DEEPVARIANT_IMAGE}"
        for tool in fastp bwa samtools bcftools mosdepth gatk whatshap docker python3; do
            if command -v "${tool}" >/dev/null 2>&1; then
                version=$("${tool}" --version 2>&1 | awk 'NR==1 {print; exit}' || true)
                printf '%s\t%s\n' "${tool}" "${version:-version unavailable}"
            else
                printf '%s\tNOT_FOUND\n' "${tool}"
            fi
        done
    } > "${output}"
}
