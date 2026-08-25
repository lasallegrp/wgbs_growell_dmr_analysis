#!/bin/bash

# ============================================================================
# run_enrichment.sbatch
#
# Run enrichment for three DMR timepoints as a SLURM array.
#
# Submit all timepoints:
#   mkdir -p 00_std_err_logs
#   sbatch run_enrichment.sbatch
#
# Override plotting mode:
#   sbatch --export=ALL,MODE=publication run_enrichment.sbatch
#
# Run one timepoint:
#   sbatch --array=0 run_enrichment.sbatch   # 36-38wk
#   sbatch --array=1 run_enrichment.sbatch   # Baseline
#   sbatch --array=2 run_enrichment.sbatch   # Postpartum
# ============================================================================

#SBATCH --job-name=DMRichR_enrichment
#SBATCH --mail-user=gekuodza@ucdavis.edu
#SBATCH --account=lasallegrp
#SBATCH --partition=high
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=50G
#SBATCH --time=1-00:00:00
#SBATCH --array=0-2
#SBATCH --output=/quobyte/lasallegrp/projects/GROWELL/WGBS/2025_DBS_dmrich_George/results/logs/00_std_err_logs/DMRichR_enrichment_%A_%a.out
#SBATCH --error=/quobyte/lasallegrp/projects/GROWELL/WGBS/2025_DBS_dmrich_George/results/logs/00_std_err_logs/DMRichR_enrichment_%A_%a.err
#SBATCH --mail-type=ALL

set -euo pipefail

# ============================================================================
# Run information
# ============================================================================

START_TIME=$(date +%s)

echo "Job ID: ${SLURM_JOB_ID}"
echo "Array task ID: ${SLURM_ARRAY_TASK_ID}"
echo "Host: $(hostname)"
echo "Working directory: $(pwd)"
echo "Start time: $(date)"

# ============================================================================
# Pixi environment
# ============================================================================

export PIXI_HOME="/quobyte/lasallegrp/home/George/.pixi"
export PATH="${PIXI_HOME}/bin:${PATH}"

export PIXI_CACHE_DIR="/tmp/${USER}/pixi-cache"
mkdir -p "${PIXI_CACHE_DIR}"

THREADS="${SLURM_CPUS_PER_TASK:-1}"

#the environment variable PIXI_PROJECT is set to the reproducibility/pixi folder of the comethyl project, because the dmrichr version is old.
# So KEGG with clusterProfiler will not work with the dmrichr version. 
# The comethyl version is up to date and works with clusterProfiler KEGG enrichment.
PIXI_PROJECT="/quobyte/lasallegrp/projects/GROWELL/WGBS/2025_DBS_comethyl_George/reproducibility/pixi"
# PIXI_PROJECT="/quobyte/lasallegrp/projects/GROWELL/WGBS/2025_DBS_dmrich_George/reproducibility/pixi"

PIXI_MANIFEST="${PIXI_PROJECT}/pixi.toml"

if [[ ! -f "${PIXI_MANIFEST}" ]]; then
    echo "ERROR: Pixi manifest not found: ${PIXI_MANIFEST}" >&2
    exit 1
fi

echo "Allocated threads: ${THREADS}"

# Optional environment check
pixi run \
    --manifest-path "${PIXI_MANIFEST}" \
    Rscript -e 'cat("R home:", R.home(), "\n"); sessionInfo()'

# ============================================================================
# Project paths
# ============================================================================

SCRIPT_DIR="/quobyte/lasallegrp/projects/GROWELL/WGBS/2025_DBS_dmrich_George/scripts/DMRichr/run_in_DMRs_folder"
ENRICHMENT_SCRIPT="${SCRIPT_DIR}/02_enrichment_plot.R"

# Select one analysis
# BASE="/quobyte/lasallegrp/projects/GROWELL/WGBS/2025_DBS_dmrich_George/results/egwg"

# Other options:
# BASE="/quobyte/lasallegrp/projects/GROWELL/WGBS/2025_DBS_dmrich_George/results/apo"
# BASE="/quobyte/lasallegrp/projects/GROWELL/WGBS/2025_DBS_dmrich_George/results/ppwr"
# BASE="/quobyte/lasallegrp/projects/GROWELL/WGBS/2025_DBS_dmrich_George/results/experimental_group"

if [[ ! -f "${ENRICHMENT_SCRIPT}" ]]; then
    echo "ERROR: enrichment script not found: ${ENRICHMENT_SCRIPT}" >&2
    exit 1
fi

# ============================================================================
# Array inputs
# ============================================================================

GENE_FILES=(
    "${BASE}/36-38wk/DMRs/DMRs_annotated.xlsx"
    "${BASE}/Baseline/DMRs/DMRs_annotated.xlsx"
    "${BASE}/Postpartum/DMRs/DMRs_annotated.xlsx"
)

LABELS=(
    "36-38wk"
    "Baseline"
    "Postpartum"
)

TASK_ID="${SLURM_ARRAY_TASK_ID}"

if (( TASK_ID < 0 || TASK_ID >= ${#GENE_FILES[@]} )); then
    echo "ERROR: invalid array task ID: ${TASK_ID}" >&2
    exit 1
fi

GENES="${GENE_FILES[$TASK_ID]}"
LABEL="${LABELS[$TASK_ID]}"
OUTDIR="${BASE}/${LABEL}/enrichment"

# File prefix = <group>_<timepoint>, e.g. "egwg_36-38wk", so output files
# read egwg_36-38wk_GO.pdf, egwg_36-38wk_full_results.xlsx, etc. Group name
# is derived from BASE's directory name; if BASE is ever unset/empty, falls
# back to "results" rather than erroring, so the naming convention always
# resolves to something (e.g. "results_36-38wk").
if [[ -n "${BASE:-}" ]]; then
    GROUP=$(basename "${BASE}")
else
    GROUP="results"
    echo "WARNING: BASE not set -- defaulting group name to '${GROUP}'" >&2
fi
PREFIX="${GROUP}_${LABEL}"

# Choose either presentation or publication.
MODE="${MODE:-publication}"

if [[ "${MODE}" != "presentation" && "${MODE}" != "publication" ]]; then
    echo "ERROR: MODE must be 'presentation' or 'publication'." >&2
    echo "Received: ${MODE}" >&2
    exit 1
fi

if [[ ! -f "${GENES}" ]]; then
    echo "ERROR: gene file not found: ${GENES}" >&2
    exit 1
fi

mkdir -p "${OUTDIR}"

# ============================================================================
# Run enrichment
# ============================================================================

echo
echo "Running enrichment"
echo "  Analysis:  $(basename "${BASE}")"
echo "  Task:      ${TASK_ID}"
echo "  Timepoint: ${LABEL}"
echo "  Prefix:    ${PREFIX}"
echo "  Input:     ${GENES}"
echo "  Output:    ${OUTDIR}"
echo "  Mode:      ${MODE}"
echo

pixi run \
    --manifest-path "${PIXI_MANIFEST}" \
    Rscript --vanilla "${ENRICHMENT_SCRIPT}" \
    --genes "${GENES}" \
    --outdir "${OUTDIR}" \
    --prefix "${PREFIX}" \
    --mode "${MODE}" \
    --kegg_method "clusterProfiler" \
    --top_n 9 \
    --go_dbs "GO_Biological_Process_2025,GO_Cellular_Component_2025,GO_Molecular_Function_2025" \
    --gwas_db "GWAS_Catalog_2025"


    # --go_method enrichr \
    # --kegg_method "clusterProfiler" \
# ============================================================================
# Completion information
# ============================================================================

END_TIME=$(date +%s)
ELAPSED_SECONDS=$((END_TIME - START_TIME))

echo
echo "Task ${TASK_ID} (${PREFIX}) completed successfully."
echo "End time: $(date)"
echo "Elapsed time: ${ELAPSED_SECONDS} seconds"

#wanted to fix folder naming but this worked prior to the fix so I am leaving it as is.
# #!/bin/bash

# # ============================================================================
# # run_enrichment.sbatch
# #
# # Run enrichment for three DMR timepoints as a SLURM array.
# #
# # Submit all timepoints:
# #   mkdir -p 00_std_err_logs
# #   sbatch run_enrichment.sbatch
# #
# # Override plotting mode:
# #   sbatch --export=ALL,MODE=publication run_enrichment.sbatch
# #
# # Run one timepoint:
# #   sbatch --array=0 run_enrichment.sbatch   # 36-38wk
# #   sbatch --array=1 run_enrichment.sbatch   # Baseline
# #   sbatch --array=2 run_enrichment.sbatch   # Postpartum
# # ============================================================================

# #SBATCH --job-name=DMRichR_enrichment
# #SBATCH --mail-user=gekuodza@ucdavis.edu
# #SBATCH --account=lasallegrp
# #SBATCH --partition=high
# #SBATCH --ntasks=1
# #SBATCH --cpus-per-task=4
# #SBATCH --mem=50G
# #SBATCH --time=1-00:00:00
# #SBATCH --array=0-2
# #SBATCH --output=00_std_err_logs/DMRichR_enrichment_%A_%a.out
# #SBATCH --error=00_std_err_logs/DMRichR_enrichment_%A_%a.err
# #SBATCH --mail-type=ALL

# set -euo pipefail

# # ============================================================================
# # Run information
# # ============================================================================

# START_TIME=$(date +%s)

# echo "Job ID: ${SLURM_JOB_ID}"
# echo "Array task ID: ${SLURM_ARRAY_TASK_ID}"
# echo "Host: $(hostname)"
# echo "Working directory: $(pwd)"
# echo "Start time: $(date)"

# # ============================================================================
# # Pixi environment
# # ============================================================================

# export PIXI_HOME="/quobyte/lasallegrp/home/George/.pixi"
# export PATH="${PIXI_HOME}/bin:${PATH}"

# export PIXI_CACHE_DIR="/tmp/${USER}/pixi-cache"
# mkdir -p "${PIXI_CACHE_DIR}"

# THREADS="${SLURM_CPUS_PER_TASK:-1}"

# #the environment variable PIXI_PROJECT is set to the reproducibility/pixi folder of the comethyl project, because the dmrichr version is old.
# # So KEGG with clusterProfiler will not work with the dmrichr version. 
# # The comethyl version is up to date and works with clusterProfiler KEGG enrichment.
# PIXI_PROJECT="/quobyte/lasallegrp/projects/GROWELL/WGBS/2025_DBS_comethyl_George/reproducibility/pixi"
# # PIXI_PROJECT="/quobyte/lasallegrp/projects/GROWELL/WGBS/2025_DBS_dmrich_George/reproducibility/pixi"

# PIXI_MANIFEST="${PIXI_PROJECT}/pixi.toml"

# if [[ ! -f "${PIXI_MANIFEST}" ]]; then
#     echo "ERROR: Pixi manifest not found: ${PIXI_MANIFEST}" >&2
#     exit 1
# fi

# echo "Allocated threads: ${THREADS}"

# # Optional environment check
# pixi run \
#     --manifest-path "${PIXI_MANIFEST}" \
#     Rscript -e 'cat("R home:", R.home(), "\n"); sessionInfo()'

# # ============================================================================
# # Project paths
# # ============================================================================

# SCRIPT_DIR="/quobyte/lasallegrp/projects/GROWELL/WGBS/2025_DBS_dmrich_George/scripts/DMRichr/run_in_DMRs_folder"
# ENRICHMENT_SCRIPT="${SCRIPT_DIR}/02_enrichment_plot.R"

# # Select one analysis
# BASE="/quobyte/lasallegrp/projects/GROWELL/WGBS/2025_DBS_dmrich_George/results/egwg"

# # Other options:
# # BASE="/quobyte/lasallegrp/projects/GROWELL/WGBS/2025_DBS_dmrich_George/results/apo"
# # BASE="/quobyte/lasallegrp/projects/GROWELL/WGBS/2025_DBS_dmrich_George/results/ppwr"
# # BASE="/quobyte/lasallegrp/projects/GROWELL/WGBS/2025_DBS_dmrich_George/results/experimental_group"

# if [[ ! -f "${ENRICHMENT_SCRIPT}" ]]; then
#     echo "ERROR: enrichment script not found: ${ENRICHMENT_SCRIPT}" >&2
#     exit 1
# fi

# # ============================================================================
# # Array inputs
# # ============================================================================

# GENE_FILES=(
#     "${BASE}/36-38wk/DMRs/DMRs_annotated.xlsx"
#     "${BASE}/Baseline/DMRs/DMRs_annotated.xlsx"
#     "${BASE}/Postpartum/DMRs/DMRs_annotated.xlsx"
# )

# LABELS=(
#     "36-38wk"
#     "Baseline"
#     "Postpartum"
# )

# TASK_ID="${SLURM_ARRAY_TASK_ID}"

# if (( TASK_ID < 0 || TASK_ID >= ${#GENE_FILES[@]} )); then
#     echo "ERROR: invalid array task ID: ${TASK_ID}" >&2
#     exit 1
# fi

# GENES="${GENE_FILES[$TASK_ID]}"
# LABEL="${LABELS[$TASK_ID]}"
# OUTDIR="${BASE}/${LABEL}/enrichment"

# # Choose either presentation or publication.
# MODE="${MODE:-publication}"

# if [[ "${MODE}" != "presentation" && "${MODE}" != "publication" ]]; then
#     echo "ERROR: MODE must be 'presentation' or 'publication'." >&2
#     echo "Received: ${MODE}" >&2
#     exit 1
# fi

# if [[ ! -f "${GENES}" ]]; then
#     echo "ERROR: gene file not found: ${GENES}" >&2
#     exit 1
# fi

# mkdir -p "${OUTDIR}"

# # ============================================================================
# # Run enrichment
# # ============================================================================

# echo
# echo "Running enrichment"
# echo "  Analysis:  $(basename "${BASE}")"
# echo "  Task:      ${TASK_ID}"
# echo "  Timepoint: ${LABEL}"
# echo "  Input:     ${GENES}"
# echo "  Output:    ${OUTDIR}"
# echo "  Mode:      ${MODE}"
# echo

# pixi run \
#     --manifest-path "${PIXI_MANIFEST}" \
#     Rscript --vanilla "${ENRICHMENT_SCRIPT}" \
#     --genes "${GENES}" \
#     --outdir "${OUTDIR}" \
#     --prefix "${LABEL}" \
#     --mode "${MODE}" \
#     --kegg_method "clusterProfiler" \
#     --top_n 7 \
#     --go_dbs "GO_Biological_Process_2025,GO_Cellular_Component_2025,GO_Molecular_Function_2025" \
#     --gwas_db "GWAS_Catalog_2025"


#     # --go_method enrichr \
#     # --kegg_method "clusterProfiler" \
# # ============================================================================
# # Completion information
# # ============================================================================

# END_TIME=$(date +%s)
# ELAPSED_SECONDS=$((END_TIME - START_TIME))

# echo
# echo "Task ${TASK_ID} (${LABEL}) completed successfully."
# echo "End time: $(date)"
# echo "Elapsed time: ${ELAPSED_SECONDS} seconds"