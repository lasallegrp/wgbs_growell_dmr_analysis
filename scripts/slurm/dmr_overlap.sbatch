#!/bin/bash

# ============================================================================
# dmr_overlap.sbatch
#
# Run 03_dmrichr_overlap.R on the GROWELL outcome x timepoint layout:
#   results/dmrichr_output/<outcome>/<timepoint>/DMRs/DMRs_annotated.xlsx
# Sets are discovered automatically (no config table needed).
#
# Comparisons:
#   within_outcome_across_timepoints   e.g. egwg: Baseline / 36-38wk / Postpartum
#   within_timepoint_across_outcomes   e.g. Baseline: egwg / ppwr / apo / ...
# Overlap by coordinates and/or gene symbol, optional regioneR permutation test.
# Output: results/dmrichr_output/DMRichR_overlap/<RUN_TAG>/
#
# Submit FROM THE PROJECT ROOT (the log paths below are relative to it):
#   cd /quobyte/lasallegrp/projects/GROWELL/WGBS/<project>
#   sbatch scripts/slurm/dmr_overlap.sbatch
#
# Optional overrides at submit time, e.g.:
#   sbatch --export=ALL,OUTCOMES=egwg,ppwr,COMPARISONS=within_timepoint_across_outcomes scripts/slurm/dmr_overlap.sbatch
#   sbatch --export=ALL,PERM_TEST=1 scripts/slurm/dmr_overlap.sbatch                       # 100 x resampleRegions
#   sbatch --export=ALL,PERM_TEST=1,PERM_N=1000 scripts/slurm/dmr_overlap.sbatch           # final numbers
#   sbatch --export=ALL,PERM_TEST=1,PERM_RANDOMIZE=randomize,PERM_GENOME=hg38 scripts/slurm/dmr_overlap.sbatch
#   sbatch --export=ALL,REQUIRE_DIRECTION=1,EXCLUDE_SEX_CHR=1 scripts/slurm/dmr_overlap.sbatch
#   sbatch --export=ALL,OVERLAP_MODE=coord,MAX_GAP=500 scripts/slurm/dmr_overlap.sbatch
#   sbatch --export=ALL,UPSET_GENES=5 scripts/slurm/dmr_overlap.sbatch
#   sbatch --export=ALL,OUTCOMES=egwg:ppwr,COMPARISONS=within_timepoint_across_outcomes,PERM_TEST=1        scripts/slurm/dmr_overlap.sbatch
#   sbatch --export=ALL,OUTCOMES=egwg:ppwr,COMPARISONS=within_outcome_across_timepoints,PERM_TEST=1 scripts/slurm/dmr_overlap.sbatch
#   sbatch --export=ALL,COMPARISONS=within_timepoint_across_outcomes,PERM_TEST=1        scripts/slurm/dmr_overlap.sbatch
#   sbatch --export=ALL,COMPARISONS=within_outcome_across_timepoints,PERM_TEST=1 scripts/slurm/dmr_overlap.sbatch
# Note: --export values are comma-separated, so for OUTCOMES/TIMEPOINTS with
# several entries use ":" instead of "," (e.g. OUTCOMES=egwg:ppwr); it is
# converted below.
# ============================================================================

#SBATCH --job-name=DMRichR_overlap
#SBATCH --mail-user=gekuodza@ucdavis.edu
#SBATCH --account=lasallegrp
#SBATCH --partition=high
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=03:00:00
#SBATCH --output=results/logs/DMRichR_overlap_%j.out
#SBATCH --error=results/logs/DMRichR_overlap_%j.err
#SBATCH --mail-type=ALL

set -euo pipefail

# ============================================================================
# Run information
# ============================================================================

START_TIME=$(date +%s)

echo "Job ID: ${SLURM_JOB_ID:-interactive}"
echo "Host: $(hostname)"
echo "Working directory: $(pwd)"
echo "Start time: $(date)"

# ============================================================================
# Project paths
# ============================================================================

# Project root = submit directory unless PROJECT_DIR is given
PROJECT_DIR="${PROJECT_DIR:-${SLURM_SUBMIT_DIR:-$(pwd)}}"
PROJECT_DIR="$(cd "${PROJECT_DIR}" && pwd)"

SCRIPT_DIR="${SCRIPT_DIR:-${PROJECT_DIR}/scripts/DMRichr}"      # 03_dmrichr_overlap.R + helper_overlap.R
OVERLAP_SCRIPT="${SCRIPT_DIR}/03_dmrichr_overlap.R"
HELPER="${SCRIPT_DIR}/helper_overlap.R"
RESULTS_DIR="${RESULTS_DIR:-${PROJECT_DIR}/results/dmrichr_output}"   # holds <outcome>/<timepoint>/

if [[ ! -d "${RESULTS_DIR}" ]]; then
    echo "ERROR: results folder not found: ${RESULTS_DIR}" >&2
    echo "       Submit from the project root, or set PROJECT_DIR / RESULTS_DIR." >&2
    exit 1
fi
for f in "${OVERLAP_SCRIPT}" "${HELPER}"; do
    if [[ ! -f "${f}" ]]; then
        echo "ERROR: file not found: ${f}" >&2
        exit 1
    fi
done

# ============================================================================
# Analysis settings (override with sbatch --export=ALL,VAR=value)
# ============================================================================

OUTCOMES="${OUTCOMES:-}"                            # empty = all outcome folders; e.g. egwg:ppwr
TIMEPOINTS="${TIMEPOINTS:-Baseline:36-38wk:Postpartum}"   # order used in tables/plots
SUBDIR="${SUBDIR:-}"                                # e.g. cutoff_005 if there is a level above DMRs/
EXCLUDE_DIRS="${EXCLUDE_DIRS:-logs:metadata:DMRichR_overlap}"
COMPARISONS="${COMPARISONS:-within_outcome_across_timepoints:within_timepoint_across_outcomes}"
OVERLAP_MODE="${OVERLAP_MODE:-both}"                # coord | gene | both
MAX_GAP="${MAX_GAP:-0}"                             # bp; 0 = DMRs must overlap
MIN_OVERLAP_BP="${MIN_OVERLAP_BP:-1}"
BACKGROUND_FILE="${BACKGROUND_FILE:-DMRs/background_annotated.xlsx}"   # or 'none'
REQUIRE_DIRECTION="${REQUIRE_DIRECTION:-0}"         # 1 = count only same-direction overlaps
EXCLUDE_SEX_CHR="${EXCLUDE_SEX_CHR:-0}"             # 1 = drop chrX/chrY
PER_OVERLAP_FOLDERS="${PER_OVERLAP_FOLDERS:-TRUE}"  # each overlap's annotated DMRs in its own folder
UPSET_GENES="${UPSET_GENES:-5}"                     # gene names per bar on the coloured UpSet (0 = none)
UPSET_GENES_SCOPE="${UPSET_GENES_SCOPE:-multi}"     # multi = overlap bars only | all
UPSET_HIGHLIGHT_GENES="${UPSET_HIGHLIGHT_GENES:-}"  # optional gene list file; only these get labeled

# regioneR permutation test (coord mode)
PERM_TEST="${PERM_TEST:-0}"                         # 1 = run
PERM_N="${PERM_N:-100}"                             # min p = 1/(PERM_N+1); use >= 1000 for final numbers
PERM_RANDOMIZE="${PERM_RANDOMIZE:-resample}"        # resample (resampleRegions) | randomize (randomizeRegions)
PERM_GENOME="${PERM_GENOME:-hg38}"                  # randomize only: hg38 (needs BSgenome) or chrom.sizes file
PERM_MASK="${PERM_MASK:-}"                          # randomize only: BED of regions to exclude
PERM_ALTERNATIVE="${PERM_ALTERNATIVE:-greater}"     # greater | less | auto
PERM_LOCAL_Z="${PERM_LOCAL_Z:-TRUE}"
PERM_SEED="${PERM_SEED:-1}"
PERM_CORES="${PERM_CORES:-${SLURM_CPUS_PER_TASK:-1}}"

# ":" -> "," for list-valued settings (sbatch --export splits on commas)
OUTCOMES="${OUTCOMES//:/,}"
TIMEPOINTS="${TIMEPOINTS//:/,}"
EXCLUDE_DIRS="${EXCLUDE_DIRS//:/,}"
COMPARISONS="${COMPARISONS//:/,}"

# Separate output folders for any non-default setting, so runs never overwrite each other
RUN_TAG="${OUTCOMES:+${OUTCOMES//,/_}}"
RUN_TAG="${RUN_TAG:-all_outcomes}"
[[ -n "${SUBDIR}" ]]                 && RUN_TAG="${RUN_TAG}_${SUBDIR}"
[[ "${REQUIRE_DIRECTION}" == "1" ]]  && RUN_TAG="${RUN_TAG}_samedir"
[[ "${EXCLUDE_SEX_CHR}" == "1" ]]    && RUN_TAG="${RUN_TAG}_noXY"
[[ "${MAX_GAP}" != "0" ]]            && RUN_TAG="${RUN_TAG}_gap${MAX_GAP}"
[[ "${MIN_OVERLAP_BP}" != "1" ]]     && RUN_TAG="${RUN_TAG}_minov${MIN_OVERLAP_BP}"
[[ "${PERM_TEST}" == "1" ]]          && RUN_TAG="${RUN_TAG}_perm${PERM_N}_${PERM_RANDOMIZE}"
OUTDIR="${OUTDIR:-${RESULTS_DIR}/DMRichR_overlap/${RUN_TAG}}"
mkdir -p "${OUTDIR}"

# ============================================================================
# Pixi environment
# ============================================================================

export PIXI_HOME="/quobyte/lasallegrp/home/George/.pixi"
export PATH="${PIXI_HOME}/bin:${PATH}"

export PIXI_CACHE_DIR="/tmp/${USER}/pixi-cache"
mkdir -p "${PIXI_CACHE_DIR}"

PIXI_MANIFEST="${PIXI_MANIFEST:-${PROJECT_DIR}/reproducibility/pixi/pixi.toml}"

if [[ ! -f "${PIXI_MANIFEST}" ]]; then
    echo "ERROR: Pixi manifest not found: ${PIXI_MANIFEST}" >&2
    exit 1
fi

# Environment check: fail now if a required package is missing,
# rather than after reading all DMR sets.
NEED_REGIONER=0; [[ "${PERM_TEST}" == "1" ]] && NEED_REGIONER=1
NEED_BSGENOME=""
if [[ "${PERM_TEST}" == "1" && "${PERM_RANDOMIZE}" == "randomize" && ! -f "${PERM_GENOME}" ]]; then
    NEED_BSGENOME="BSgenome.Hsapiens.UCSC.${PERM_GENOME}"
fi
export NEED_REGIONER NEED_BSGENOME

pixi run --manifest-path "${PIXI_MANIFEST}" Rscript --vanilla -e '
cat("R home:", R.home(), "\n")
req <- c("optparse", "GenomicRanges", "openxlsx", "ggplot2")
if (Sys.getenv("NEED_REGIONER") == "1") req <- c(req, "regioneR")
if (nzchar(Sys.getenv("NEED_BSGENOME"))) req <- c(req, Sys.getenv("NEED_BSGENOME"))
opt <- c("UpSetR", "ggrepel")
miss <- req[!vapply(req, requireNamespace, logical(1), quietly = TRUE)]
if (length(miss)) stop("Missing required R packages: ", paste(miss, collapse = ", "),
                       "\n  (pixi: bioconductor-regioner / bioconductor-bsgenome.hsapiens.ucsc.hg38)")
mo <- opt[!vapply(opt, requireNamespace, logical(1), quietly = TRUE)]
if (length(mo)) message("Optional packages missing (UpSet plots / labels skipped): ", paste(mo, collapse = ", "))
cat("Package check OK\n")'

# ============================================================================
# Run overlap
# ============================================================================

echo
echo "Running DMR overlap"
echo "  Project:        ${PROJECT_DIR}"
echo "  Results dir:    ${RESULTS_DIR}"
echo "  Outcomes:       ${OUTCOMES:-all}"
echo "  Timepoints:     ${TIMEPOINTS}"
echo "  Subdir:         ${SUBDIR:-none}"
echo "  Comparisons:    ${COMPARISONS}"
echo "  Overlap mode:   ${OVERLAP_MODE}"
echo "  Max gap / min overlap (bp): ${MAX_GAP} / ${MIN_OVERLAP_BP}"
echo "  Same direction: ${REQUIRE_DIRECTION}"
echo "  Exclude chrX/Y: ${EXCLUDE_SEX_CHR}"
echo "  Background:     ${BACKGROUND_FILE}"
echo "  UpSet genes:    ${UPSET_GENES} (${UPSET_GENES_SCOPE})${UPSET_HIGHLIGHT_GENES:+, highlight: ${UPSET_HIGHLIGHT_GENES}}"
if [[ "${PERM_TEST}" == "1" ]]; then
echo "  regioneR:       ${PERM_N} x ${PERM_RANDOMIZE}Regions, alternative=${PERM_ALTERNATIVE}, cores=${PERM_CORES}, seed=${PERM_SEED}"
[[ "${PERM_RANDOMIZE}" == "randomize" ]] && \
echo "                  genome=${PERM_GENOME}${PERM_MASK:+, mask=${PERM_MASK}}"
else
echo "  regioneR:       off (PERM_TEST=1 to run)"
fi
echo "  Output:         ${OUTDIR}"
echo

cmd=(pixi run --manifest-path "${PIXI_MANIFEST}" Rscript --vanilla \
    "${OVERLAP_SCRIPT}" \
    --helper "${HELPER}" \
    --layout outcome_timepoint \
    --project_dir "${PROJECT_DIR}" \
    --results_dir "${RESULTS_DIR}" \
    --timepoints "${TIMEPOINTS}" \
    --exclude_dirs "${EXCLUDE_DIRS}" \
    --comparisons "${COMPARISONS}" \
    --overlap_mode "${OVERLAP_MODE}" \
    --max_gap "${MAX_GAP}" \
    --min_overlap_bp "${MIN_OVERLAP_BP}" \
    --background_file "${BACKGROUND_FILE}" \
    --per_overlap_folders "${PER_OVERLAP_FOLDERS}" \
    --upset_genes "${UPSET_GENES}" \
    --upset_genes_scope "${UPSET_GENES_SCOPE}" \
    --out_dir "${OUTDIR}" )

[[ -n "${OUTCOMES}" ]]              && cmd+=(--outcomes "${OUTCOMES}")
[[ -n "${SUBDIR}" ]]                && cmd+=(--subdir "${SUBDIR}")
[[ "${REQUIRE_DIRECTION}" == "1" ]] && cmd+=(--require_direction)
[[ "${EXCLUDE_SEX_CHR}" == "1" ]]   && cmd+=(--exclude_sex_chr)
[[ -n "${UPSET_HIGHLIGHT_GENES}" ]] && cmd+=(--upset_highlight_genes "${UPSET_HIGHLIGHT_GENES}")
if [[ "${PERM_TEST}" == "1" ]]; then
    cmd+=(--perm_test --perm_n "${PERM_N}" --perm_randomize "${PERM_RANDOMIZE}" \
          --perm_alternative "${PERM_ALTERNATIVE}" --perm_local_z "${PERM_LOCAL_Z}" \
          --perm_seed "${PERM_SEED}" --perm_cores "${PERM_CORES}")
    if [[ "${PERM_RANDOMIZE}" == "randomize" ]]; then
        cmd+=(--perm_genome "${PERM_GENOME}")
        [[ -n "${PERM_MASK}" ]] && cmd+=(--perm_mask "${PERM_MASK}")
    fi
fi

printf '%q ' "${cmd[@]}"; echo
"${cmd[@]}"

# ============================================================================
# Completion information
# ============================================================================

END_TIME=$(date +%s)
ELAPSED_SECONDS=$((END_TIME - START_TIME))

echo
echo "DMR overlap (${RUN_TAG}) completed successfully."
echo "Sets used: ${OUTDIR}/sets_used.tsv"
[[ "${PERM_TEST}" == "1" ]] && echo "regioneR summary: ${OUTDIR}/regioneR_permutation_summary.pdf"
echo "End time: $(date)"
echo "Elapsed time: ${ELAPSED_SECONDS} seconds"
