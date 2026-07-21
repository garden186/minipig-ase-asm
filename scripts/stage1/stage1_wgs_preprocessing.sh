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
[[ "${SUMMARY_ONLY}" -eq 0 ]] || die "--summary-only is only valid for the WhatsHap script"
load_stage1_config

LOG_DIR="${PROJECT_DIR}/logs/${SAMPLE}"
MASTER_LOG="${LOG_DIR}/stage1.log"
VERSION_LOG="${LOG_DIR}/stage1_software_versions.txt"
mkdir -p "${LOG_DIR}"
write_software_versions "${VERSION_LOG}"

master_log() {
    echo "[$(date '+%F %T')] $*" | tee -a "${MASTER_LOG}"
}

run_step() {
    local name=$1
    local script=$2
    local step_log="${LOG_DIR}/stage1_${name}.log"
    local start end elapsed
    start=$(date +%s)
    master_log ">>> ${name} starting"
    if bash "${script}" \
        --sample "${SAMPLE}" \
        --config "${CONFIG_FILE}" \
        --threads "${THREADS}" > "${step_log}" 2>&1; then
        end=$(date +%s)
        elapsed=$((end - start))
        master_log "<<< ${name} complete ($((elapsed / 3600))h$(((elapsed % 3600) / 60))m)"
    else
        local rc=$?
        end=$(date +%s)
        elapsed=$((end - start))
        master_log "<<< ${name} FAILED (exit ${rc}; $((elapsed / 3600))h$(((elapsed % 3600) / 60))m)"
        master_log "See ${step_log}"
        tail -n 20 "${step_log}" | sed 's/^/    /' | tee -a "${MASTER_LOG}"
        return "${rc}"
    fi
}

master_log "Stage 1 started for ${SAMPLE}"
master_log "Project: ${PROJECT_DIR}"
master_log "Config: ${CONFIG_FILE}"

run_step fastp "${SCRIPT_DIR}/run_fastp_wgs.sh"
run_step align "${SCRIPT_DIR}/wgs_align.sh"
run_step deepvariant "${SCRIPT_DIR}/run_deepvariant.sh"
run_step het_vcf "${SCRIPT_DIR}/prepare_het_vcf.sh"
run_step whatshap "${SCRIPT_DIR}/run_whatshap_wgs_honest.sh"
run_step validate "${SCRIPT_DIR}/validate_stage1_outputs.sh"

master_log "Stage 1 complete for ${SAMPLE}"
master_log "Primary output: ${PROJECT_DIR}/results/phasing/${SAMPLE}.phased.wgs.honest.vcf.gz"
