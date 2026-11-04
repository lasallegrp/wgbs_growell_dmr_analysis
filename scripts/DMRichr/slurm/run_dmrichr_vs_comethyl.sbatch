#!/bin/bash

# ============================================================================
# run_dmrichr_vs_comethyl.sbatch
#
# Run 04_dmrichr_vs_comethyl.R: DMRichR DMRs (outcome x timepoint) vs comethyl
# consensus modules, each DMR set paired with its own timepoint's kME and
# ME-trait table.
#
# Submit FROM THE PROJECT ROOT (log paths are relative to it):
#   cd /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_DBS_dmrich_George
#   sbatch scripts/slurm/run_dmrichr_vs_comethyl.sbatch
#
# Optional overrides (lists use ":" instead of ","):
#   sbatch --export=ALL,OUTCOMES=egwg:ppwr scripts/slurm/run_dmrichr_vs_comethyl.sbatch
#   sbatch --export=ALL,REGION_FILTER=covMin10_methSD0p09 scripts/slurm/run_dmrichr_vs_comethyl.sbatch
#   sbatch --export=ALL,COMETHYL_ROOT=/quobyte/lasallegrp/projects/GROWELL/WGBS/2025_DBS_dmrich_George/results/comethyl_output/consensus scripts/slurm/run_dmrichr_vs_comethyl.sbatch
#   sbatch --export=ALL,EXCLUDE_SEX_CHR=1,OVERLAP_MODE=coord scripts/slurm/run_dmrichr_vs_comethyl.sbatch
# ============================================================================

#SBATCH --job-name=DMRichR_vs_comethyl
#SBATCH --mail-user=gekuodza@ucdavis.edu
#SBATCH --account=lasallegrp
#SBATCH --partition=high
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=48G
#SBATCH --time=04:00:00
#SBATCH --output=results/logs/DMRichR_vs_comethyl_%j.out
#SBATCH --error=results/logs/DMRichR_vs_comethyl_%j.err
#SBATCH --mail-type=ALL

set -euo pipefail

START_TIME=$(date +%s)
echo "Job ID: ${SLURM_JOB_ID:-interactive}"
echo "Host: $(hostname)"
echo "Working directory: $(pwd)"
echo "Start time: $(date)"

# ============================================================================
# Project paths
# ============================================================================
PROJECT_DIR="${PROJECT_DIR:-${SLURM_SUBMIT_DIR:-$(pwd)}}"
PROJECT_DIR="$(cd "${PROJECT_DIR}" && pwd)"
SCRIPT_DIR="${SCRIPT_DIR:-${PROJECT_DIR}/scripts/DMRichr}"      # 04_dmrichr_vs_comethyl.R + helper_overlap.R
SCRIPT="${SCRIPT_DIR}/04_dmrichr_vs_comethyl.R"
HELPER="${SCRIPT_DIR}/helper_overlap.R"
RESULTS_DIR="${RESULTS_DIR:-${PROJECT_DIR}/results/dmrichr_output}"
COMETHYL_ROOT="${COMETHYL_ROOT:-}"                               # empty = auto-detect under results/

for f in "${SCRIPT}" "${HELPER}"; do
    [[ -f "${f}" ]] || { echo "ERROR: file not found: ${f}" >&2; exit 1; }
done
[[ -d "${RESULTS_DIR}" ]] || { echo "ERROR: results folder not found: ${RESULTS_DIR}" >&2; exit 1; }

# ============================================================================
# Analysis settings (override with sbatch --export=ALL,VAR=value)
# ============================================================================
OUTCOMES="${OUTCOMES:-}"                                         # empty = all; e.g. egwg:ppwr
TIMEPOINTS="${TIMEPOINTS:-Baseline:36-38wk:Postpartum}"
TRAIT_MAP="${TRAIT_MAP:-apo=apo:egwg=egwg:ppwr=ppwr_e:experimental_group=group}"
CPG_FILTER="${CPG_FILTER:-cov3_75pct}"
REGION_FILTER="${REGION_FILTER:-covMin14_methSD0p06}"
VARIANT="${VARIANT:-v2_exclude_protected_pcs}"
OVERLAP_MODE="${OVERLAP_MODE:-both}"                             # coord | gene | both
MAX_GAP="${MAX_GAP:-0}"
MIN_OVERLAP_BP="${MIN_OVERLAP_BP:-1}"
BACKGROUND_FILE="${BACKGROUND_FILE:-DMRs/background_annotated.xlsx}"
EXCLUDE_SEX_CHR="${EXCLUDE_SEX_CHR:-0}"
FDR="${FDR:-0.05}"
MODULE_TRAIT_P="${MODULE_TRAIT_P:-0.05}"
MIN_DMRS_PLOT="${MIN_DMRS_PLOT:-3}"
MODULE_GENES="${MODULE_GENES:-10}"

OUTCOMES="${OUTCOMES//:/,}"; TIMEPOINTS="${TIMEPOINTS//:/,}"; TRAIT_MAP="${TRAIT_MAP//:/,}"

SUFFIX=""
[[ -n "${OUTCOMES}" ]]           && SUFFIX="${SUFFIX}_${OUTCOMES//,/_}"
[[ "${EXCLUDE_SEX_CHR}" == "1" ]] && SUFFIX="${SUFFIX}_noXY"
[[ "${MAX_GAP}" != "0" ]]         && SUFFIX="${SUFFIX}_gap${MAX_GAP}"
[[ "${MIN_OVERLAP_BP}" != "1" ]]  && SUFFIX="${SUFFIX}_minov${MIN_OVERLAP_BP}"

# ============================================================================
# Pixi environment
# ============================================================================
export PIXI_HOME="/quobyte/lasallegrp/home/George/.pixi"
export PATH="${PIXI_HOME}/bin:${PATH}"
export PIXI_CACHE_DIR="/tmp/${USER}/pixi-cache"; mkdir -p "${PIXI_CACHE_DIR}"
PIXI_MANIFEST="${PIXI_MANIFEST:-${PROJECT_DIR}/reproducibility/pixi/pixi.toml}"
[[ -f "${PIXI_MANIFEST}" ]] || { echo "ERROR: Pixi manifest not found: ${PIXI_MANIFEST}" >&2; exit 1; }

pixi run --manifest-path "${PIXI_MANIFEST}" Rscript --vanilla -e '
cat("R home:", R.home(), "\n")
req <- c("optparse", "GenomicRanges", "openxlsx", "ggplot2")
miss <- req[!vapply(req, requireNamespace, logical(1), quietly = TRUE)]
if (length(miss)) stop("Missing required R packages: ", paste(miss, collapse = ", "))
opt <- c("ggrepel", "data.table")
mo <- opt[!vapply(opt, requireNamespace, logical(1), quietly = TRUE)]
if (length(mo)) message("Optional packages missing (labels / fast reading of the large membership tables): ",
                        paste(mo, collapse = ", "))
cat("Package check OK\n")'

# ============================================================================
# Run
# ============================================================================
echo
echo "Running DMRichR vs comethyl consensus"
echo "  Results dir:    ${RESULTS_DIR}"
echo "  Comethyl root:  ${COMETHYL_ROOT:-auto}"
echo "  Network:        ${CPG_FILTER} / ${REGION_FILTER} / ${VARIANT}"
echo "  Outcomes:       ${OUTCOMES:-all}   Timepoints: ${TIMEPOINTS}"
echo "  Trait map:      ${TRAIT_MAP}"
echo "  Overlap:        ${OVERLAP_MODE}, gap ${MAX_GAP}, min overlap ${MIN_OVERLAP_BP}, noXY=${EXCLUDE_SEX_CHR}"
echo

cmd=(pixi run --manifest-path "${PIXI_MANIFEST}" Rscript --vanilla "${SCRIPT}" \
    --helper "${HELPER}" \
    --project_dir "${PROJECT_DIR}" \
    --results_dir "${RESULTS_DIR}" \
    --timepoints "${TIMEPOINTS}" \
    --trait_map "${TRAIT_MAP}" \
    --cpg_filter "${CPG_FILTER}" \
    --region_filter "${REGION_FILTER}" \
    --variant "${VARIANT}" \
    --overlap_mode "${OVERLAP_MODE}" \
    --max_gap "${MAX_GAP}" \
    --min_overlap_bp "${MIN_OVERLAP_BP}" \
    --background_file "${BACKGROUND_FILE}" \
    --fdr "${FDR}" \
    --module_trait_p "${MODULE_TRAIT_P}" \
    --min_dmrs_plot "${MIN_DMRS_PLOT}" \
    --module_genes "${MODULE_GENES}" \
    --out_suffix "${SUFFIX}" )
[[ -n "${OUTCOMES}" ]]            && cmd+=(--outcomes "${OUTCOMES}")
[[ -n "${COMETHYL_ROOT}" ]]       && cmd+=(--comethyl_root "${COMETHYL_ROOT}")
[[ "${EXCLUDE_SEX_CHR}" == "1" ]] && cmd+=(--exclude_sex_chr)

printf '%q ' "${cmd[@]}"; echo
"${cmd[@]}"

END_TIME=$(date +%s)
echo
echo "DMRichR vs comethyl completed. Elapsed: $((END_TIME - START_TIME)) s"
echo "Output: ${RESULTS_DIR}/DMRichR_vs_comethyl/${REGION_FILTER}/${VARIANT}${SUFFIX}"
