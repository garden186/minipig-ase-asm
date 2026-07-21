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
require_command fastp
require_command python3

OUT_DIR="${PROJECT_DIR}/results/wgs/trim/${SAMPLE}"
R1_OUT="${OUT_DIR}/${SAMPLE}_R1.fastq.gz"
R2_OUT="${OUT_DIR}/${SAMPLE}_R2.fastq.gz"
JSON_OUT="${OUT_DIR}/${SAMPLE}.fastp.json"
HTML_OUT="${OUT_DIR}/${SAMPLE}.fastp.html"
LOG_OUT="${OUT_DIR}/${SAMPLE}.fastp.log"
mkdir -p "${OUT_DIR}"

if [[ -s "${R1_OUT}" && -s "${R2_OUT}" && -s "${JSON_OUT}" && -s "${HTML_OUT}" ]]; then
    note "[SKIP] Complete fastp output exists for ${SAMPLE}"
    exit 0
fi

TMP_PREFIX="${OUT_DIR}/.${SAMPLE}.fastp.$$"
TMP_R1="${TMP_PREFIX}.R1.fastq.gz"
TMP_R2="${TMP_PREFIX}.R2.fastq.gz"
TMP_JSON="${TMP_PREFIX}.json"
TMP_HTML="${TMP_PREFIX}.html"
trap 'rm -f "${TMP_R1}" "${TMP_R2}" "${TMP_JSON}" "${TMP_HTML}"' EXIT

note "Running fastp for ${SAMPLE}"
fastp \
    --in1 "${WGS_R1}" \
    --in2 "${WGS_R2}" \
    --out1 "${TMP_R1}" \
    --out2 "${TMP_R2}" \
    --adapter_sequence AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
    --adapter_sequence_r2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
    --cut_right \
    --cut_right_window_size 4 \
    --cut_right_mean_quality 15 \
    --qualified_quality_phred 15 \
    --unqualified_percent_limit 40 \
    --n_base_limit 5 \
    --length_required 50 \
    --trim_poly_g \
    --poly_g_min_len 10 \
    --json "${TMP_JSON}" \
    --html "${TMP_HTML}" \
    --thread "${FASTP_THREADS}" \
    2>&1 | tee "${LOG_OUT}"

for output in "${TMP_R1}" "${TMP_R2}" "${TMP_JSON}" "${TMP_HTML}"; do
    require_nonempty "${output}"
done

python3 - "${TMP_JSON}" <<'PY'
import json
import sys

with open(sys.argv[1], encoding="utf-8") as handle:
    report = json.load(handle)

before = report["summary"]["before_filtering"]
after = report["summary"]["after_filtering"]
if before["total_reads"] <= 0 or after["total_reads"] <= 0:
    raise SystemExit("fastp reported zero reads")
if after["total_reads"] > before["total_reads"]:
    raise SystemExit("fastp after-filter read count exceeds input count")
PY

mv "${TMP_R1}" "${R1_OUT}"
mv "${TMP_R2}" "${R2_OUT}"
mv "${TMP_JSON}" "${JSON_OUT}"
mv "${TMP_HTML}" "${HTML_OUT}"
trap - EXIT

note "fastp complete: ${JSON_OUT}"
