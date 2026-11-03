#!/usr/bin/env Rscript

# ================================================================
# SCRIPT 12B: Consensus ME-Trait Presentation Heatmaps
#
# Pipeline: comethyl WGBS consensus analysis
#
# PURPOSE
#   Generates focused, publication-ready ME-trait heatmaps by
#   subsetting the full stats from script 12a to user-defined
#   trait sets:
#     1) Load ME-trait stats per dataset from script 12a
#     2) For each trait-set file, subset to traits of interest
#     3) Generate per-dataset FULL, TOP-N, and FOCUSED heatmaps
#     4) Generate cross-dataset comparison heatmaps (same trait
#        set, different datasets side-by-side in cross_dataset/)
#   Module ordering (dendro or exposure/outcome) and trait ordering
#   (input or exposure/outcome) are computed ONCE per run, using
#   dataset1 as the reference dataset, and then applied consistently
#   across all datasets and the cross-dataset comparison so that
#   panels stay visually aligned for comparison.
#
# REQUIRED INPUTS
#   --project_root        : root directory of the analysis project
#   --regions_file         : path to Filtered_Regions.txt from script 02
#   --dataset1_label       : label for dataset 1
#   --dataset1_stats_file  : path to ME-trait stats TSV from script 12a
#
# OPTIONAL INPUTS
#   --consensus_modules_rds  : path to Consensus_Modules.rds (for dendro module ordering)
#   --dataset2_label         : label for dataset 2
#   --dataset2_stats_file    : stats TSV for dataset 2
#   --dataset3_label         : label for dataset 3
#   --dataset3_stats_file    : stats TSV for dataset 3
#   --adjustment_version     : label matching script 12a run [default = unadjusted]
#   --set_dir                : directory of .txt trait-set files (all processed)
#   --set_dir_recursive      : TRUE/FALSE; search --set_dir recursively [default = FALSE]
#   --set_file               : path to single trait-set .txt file
#   --set_file2 to set_file5 : additional trait-set files
#   --module_dendro_distance : bicor or pearson [default = bicor]
#
#   --module_order_mode       : default, dendro, or exposure_outcome [default = dendro]
#   --module_order_trait_file : TSV/CSV/TXT with columns category, trait.
#                                Required when --module_order_mode exposure_outcome.
#   --module_order_p_thresh   : p-value cutoff for classifying modules [default = same as --p_thresh]
#   --module_order_score      : p or abs_cor, within-class ordering [default = p]
#
#   --trait_order_mode        : input or exposure_outcome [default = input]
#   --trait_order_p_thresh    : p-value cutoff for trait classification [default = module_order_p_thresh]
#   --trait_order_score       : p or abs_cor, within-class trait ordering [default = p]
#   --write_reordered_trait_files : TRUE/FALSE; save ordered trait text files [default = TRUE]
#
#   --drop_grey_modules      : TRUE/FALSE; remove grey module from all datasets [default = TRUE]
#   --p_thresh               : significance threshold [default = 0.05]
#   --top_n                  : top N associations for TOP heatmap [default = 250]
#   --min_sig_traits         : minimum number of significant trait associations (p < p_thresh)
#                              a module must have (within a given dataset's trait subset) to
#                              appear in that dataset's FOCUSED heatmap. [default = 1]
#   --full_width/full_height : FULL heatmap dimensions in inches [default = 12 x 12]
#   --top_width/top_height   : TOP heatmap dimensions in inches [default = 9 x 6]
#   --focused_width          : FOCUSED heatmap width in inches [default = same as --full_width]
#   --focused_height         : FOCUSED heatmap height in inches [default = 9]
#
# MODULE ORDER TRAIT FILE FORMAT
#   A tab-delimited file with two required columns:
#
#     category    trait
#     exposure    DDE
#     exposure    op_DDT
#     outcome     breast_cancer
#
#   category must be one of: exposure, outcome
#
# OUTPUTS
#   comethyl_output/consensus/12b_me_trait_presentation/
#       module_order/<cpg_label>/<region_label>/<adjustment_version>/
#           ordered_modules.txt
#           module_order_exposure_outcome_summary.tsv   (only if module_order_mode = exposure_outcome)
#           <set_name>/
#               traits_ordered.txt
#               <set_name>_ORDERED.txt
#               trait_order_exposure_outcome_summary.tsv (only if trait_order_mode = exposure_outcome)
#       <dataset_label>/<cpg_label>/<region_label>/<adjustment_version>/
#           <set_name>/
#               <dataset_label>_<set_name>_ME_Trait_Heatmap_FULL.pdf
#               <dataset_label>_<set_name>_ME_Trait_Heatmap_TOP.pdf
#               <dataset_label>_<set_name>_ME_Trait_Heatmap_FOCUSED.pdf
#                   A filtered version of the FULL heatmap showing only modules that meet
#                   the --min_sig_traits threshold for that dataset's trait subset.
#                   Replaced by a _FOCUSED_SKIPPED.txt if no modules meet the threshold.
#               <dataset_label>_<set_name>_focused_modules.txt
#               <dataset_label>_<set_name>_subset_stats.tsv
#               <dataset_label>_<set_name>_subset_stats_significant.tsv
#               <dataset_label>_<set_name>_traits_requested/found/missing.txt
#               run_parameters.txt
#       cross_dataset/<cpg_label>/<region_label>/<adjustment_version>/
#           <set_name>/
#               <set_name>_Comparison_<dataset_label>_FULL.pdf
#               <set_name>_Comparison_<dataset_label>_TOP.pdf
#               <set_name>_Comparison_<dataset_label>_FOCUSED.pdf
#
# NOTES
#   - cpg_label and region_label are derived from --regions_file path
#   - Trait sets are plain text files with one trait name per line
#   - If --set_dir is provided, all .txt files in that directory are processed
#     (recursively if --set_dir_recursive TRUE)
#   - Module and trait ordering are derived from dataset1 only, then reused for every
#     dataset and the cross-dataset comparison so panels line up. This differs from
#     script 09b, which only ever had a single stats table to order against.
#
# EXAMPLE
#   Rscript /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_DBS_dmrich_George/scripts/comethyl/consensus/12b_me_trait_presentation_consensus.R \
#     --project_root /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_DBS_dmrich_George \
#     --regions_file /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_DBS_dmrich_George/results/comethyl_output/consensus/02_reference_region_filter/Baseline/cov3_75pct/covMin4_methSD0p08/Filtered_Regions.txt \
#     --consensus_modules_rds /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_DBS_dmrich_George/results/comethyl_output/consensus/09_consensus_modules/cov3_75pct/covMin4_methSD0p08/v1_all_pcs/shared/Consensus_Modules.rds \
#     --dataset1_label Baseline \
#     --dataset1_stats_file /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_DBS_dmrich_George/results/comethyl_output/consensus/12a_me_trait_analysis/cov3_75pct/covMin4_methSD0p08/v1_all_pcs/Baseline/ME_Trait_Correlation_Stats_Bicor.tsv \
#     --dataset2_label 36_38wks \
#     --dataset2_stats_file /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_DBS_dmrich_George/results/comethyl_output/consensus/12a_me_trait_analysis/cov3_75pct/covMin4_methSD0p08/v1_all_pcs/36_38wks/ME_Trait_Correlation_Stats_Bicor.tsv \
#     --adjustment_version v1_all_pcs \
#     --set_dir /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_DBS_dmrich_George/analysis/config/comethyl/trait_sets/
# ================================================================

message("Starting Script 12b \u2713")

suppressPackageStartupMessages({
  library(optparse)
  library(comethyl)
  library(WGCNA)
  library(readr)
  library(dplyr)
})

# ------------------------------------------------------------
# Load helper.R
# ------------------------------------------------------------
script_file_arg <- commandArgs(trailingOnly = FALSE)[
  grep("^--file=", commandArgs(trailingOnly = FALSE))
]

if (length(script_file_arg) == 0) {
  stop("Could not determine script path.")
}

script_dir  <- dirname(normalizePath(sub("^--file=", "", script_file_arg[1])))
helper_file <- file.path(script_dir, "helper.R")

if (!file.exists(helper_file)) {
  stop("helper.R not found: ", helper_file)
}

source(helper_file)

required_helpers <- c("write_vector_file", "write_log_lines", "getDendro", "plotMEtraitCor")
missing_helpers <- required_helpers[!vapply(required_helpers, exists, logical(1), mode = "function")]
if (length(missing_helpers) > 0) {
  stop("Missing helper functions in helper.R: ", paste(missing_helpers, collapse = ", "))
}

# ------------------------------------------------------------
# Helpers
# ------------------------------------------------------------
read_trait_set_file <- function(file) {
  x <- readLines(file, warn = FALSE)
  x <- trimws(x)
  x <- x[nzchar(x)]
  x <- x[!grepl("^#", x)]
  unique(x)
}

get_set_name <- function(file) {
  tools::file_path_sans_ext(basename(file))
}

collect_set_files <- function(set_dir, set_file, set_file2, set_file3, set_file4, set_file5,
                               recursive = FALSE) {
  direct <- c(set_file, set_file2, set_file3, set_file4, set_file5)
  direct <- direct[!is.na(direct) & nzchar(direct)]

  if (length(direct) > 0) {
    missing <- direct[!file.exists(direct)]
    if (length(missing) > 0) {
      stop("Trait set files not found:\n  ", paste(missing, collapse = "\n  "))
    }
  }

  from_dir <- character(0)
  if (!is.null(set_dir) && nzchar(set_dir)) {
    if (!dir.exists(set_dir)) stop("set_dir not found: ", set_dir)
    from_dir <- list.files(
      set_dir,
      pattern = "\\.txt$",
      full.names = TRUE,
      recursive = recursive,
      ignore.case = TRUE
    )
    if (length(from_dir) == 0) stop("No .txt files found in set_dir: ", set_dir)
  }

  all_files <- unique(normalizePath(c(direct, from_dir), mustWork = TRUE))

  if (length(all_files) == 0) {
    stop("No trait set files provided. Use --set_dir or --set_file.")
  }

  all_files
}

get_cor_col <- function(stats_df) {
  if ("bicor" %in% colnames(stats_df)) return("bicor")
  if ("cor" %in% colnames(stats_df)) return("cor")
  if ("correlation" %in% colnames(stats_df)) return("correlation")
  stop("Stats table must contain one of: bicor, cor, correlation")
}

standardize_p_col <- function(stats_df) {
  if ("p" %in% colnames(stats_df)) return(stats_df)
  candidates <- c("pvalue", "p.value", "p_value", "P.Value", "P.value")
  hit <- candidates[candidates %in% colnames(stats_df)][1]
  if (is.na(hit)) return(stats_df)
  colnames(stats_df)[match(hit, colnames(stats_df))] <- "p"
  stats_df
}

validate_stats_table <- function(df, label) {
  df <- standardize_p_col(df)
  required <- c("module", "trait", "p")
  missing  <- setdiff(required, colnames(df))
  if (length(missing) > 0) {
    stop(label, " missing columns: ", paste(missing, collapse = ", "))
  }
  get_cor_col(df)
  df
}

module_name_from_me_col <- function(x) {
  sub("^ME", "", as.character(x))
}

is_grey_module <- function(x) {
  tolower(module_name_from_me_col(x)) == "grey"
}

read_module_order_trait_file <- function(file) {
  if (is.null(file) || !nzchar(file)) {
    stop("--module_order_trait_file is required when --module_order_mode exposure_outcome")
  }
  if (!file.exists(file)) stop("module_order_trait_file not found: ", file)

  ext <- tolower(tools::file_ext(file))
  x <- switch(ext,
    csv = read.csv(file, stringsAsFactors = FALSE, check.names = FALSE),
    tsv = read.delim(file, stringsAsFactors = FALSE, check.names = FALSE),
    txt = read.delim(file, stringsAsFactors = FALSE, check.names = FALSE),
    read.delim(file, stringsAsFactors = FALSE, check.names = FALSE)
  )

  required <- c("category", "trait")
  missing <- setdiff(required, colnames(x))
  if (length(missing) > 0) {
    stop("module_order_trait_file missing required columns: ",
         paste(missing, collapse = ", "))
  }

  x$category <- tolower(trimws(as.character(x$category)))
  x$trait <- trimws(as.character(x$trait))

  x <- x[x$category %in% c("exposure", "outcome") &
         !is.na(x$trait) & x$trait != "", , drop = FALSE]

  if (nrow(x) == 0) {
    stop("module_order_trait_file has no valid exposure/outcome rows.")
  }

  unique(x[, c("category", "trait")])
}

# ------------------------------------------------------------
# Exposure/outcome module ordering
# (identical logic to script 09b; operates on a single reference
#  stats_df — here, dataset1 — so ordering is shared across datasets)
# ------------------------------------------------------------
build_exposure_outcome_module_order <- function(stats_df,
                                                 trait_def_file,
                                                 p_thresh = 0.05,
                                                 score_by = "p",
                                                 current_levels = NULL) {

  trait_def <- read_module_order_trait_file(trait_def_file)

  exposure_traits <- unique(trait_def$trait[trait_def$category == "exposure"])
  outcome_traits  <- unique(trait_def$trait[trait_def$category == "outcome"])

  if (length(exposure_traits) == 0) stop("No exposure traits defined in module_order_trait_file.")
  if (length(outcome_traits) == 0)  stop("No outcome traits defined in module_order_trait_file.")

  if (is.null(current_levels)) {
    current_levels <- unique(as.character(stats_df$module))
  }

  cor_col <- get_cor_col(stats_df)

  df <- stats_df %>%
    dplyr::mutate(
      module = as.character(module),
      trait = as.character(trait),
      p = suppressWarnings(as.numeric(p)),
      cor_value = suppressWarnings(as.numeric(.data[[cor_col]])),
      trait_class = dplyr::case_when(
        trait %in% exposure_traits ~ "exposure",
        trait %in% outcome_traits ~ "outcome",
        TRUE ~ "other"
      )
    ) %>%
    dplyr::filter(
      module %in% current_levels,
      trait_class %in% c("exposure", "outcome"),
      !is.na(p),
      p <= p_thresh
    )

  exposure_modules <- unique(df$module[df$trait_class == "exposure"])
  outcome_modules  <- unique(df$module[df$trait_class == "outcome"])

  module_summary <- data.frame(
    module = current_levels,
    module_class = dplyr::case_when(
      current_levels %in% exposure_modules & current_levels %in% outcome_modules ~ "exposure_and_outcome",
      current_levels %in% exposure_modules ~ "exposure_only",
      current_levels %in% outcome_modules ~ "outcome_only",
      TRUE ~ "other"
    ),
    stringsAsFactors = FALSE
  )

  score_df <- df %>%
    dplyr::group_by(module) %>%
    dplyr::summarise(
      min_p = min(p, na.rm = TRUE),
      max_abs_cor = max(abs(cor_value), na.rm = TRUE),
      n_sig_traits = dplyr::n_distinct(trait),
      n_sig_exposure_traits = dplyr::n_distinct(trait[trait_class == "exposure"]),
      n_sig_outcome_traits = dplyr::n_distinct(trait[trait_class == "outcome"]),
      sig_exposure_traits = paste(sort(unique(trait[trait_class == "exposure"])), collapse = ";"),
      sig_outcome_traits = paste(sort(unique(trait[trait_class == "outcome"])), collapse = ";"),
      .groups = "drop"
    )

  module_summary <- module_summary %>%
    dplyr::left_join(score_df, by = "module") %>%
    dplyr::mutate(
      min_p = ifelse(is.na(min_p), Inf, min_p),
      max_abs_cor = ifelse(is.na(max_abs_cor), 0, max_abs_cor),
      n_sig_traits = ifelse(is.na(n_sig_traits), 0, n_sig_traits),
      n_sig_exposure_traits = ifelse(is.na(n_sig_exposure_traits), 0, n_sig_exposure_traits),
      n_sig_outcome_traits = ifelse(is.na(n_sig_outcome_traits), 0, n_sig_outcome_traits),
      sig_exposure_traits = ifelse(is.na(sig_exposure_traits), "", sig_exposure_traits),
      sig_outcome_traits = ifelse(is.na(sig_outcome_traits), "", sig_outcome_traits),
      class_rank = dplyr::case_when(
        module_class == "exposure_and_outcome" ~ 1L,
        module_class == "exposure_only" ~ 2L,
        module_class == "outcome_only" ~ 3L,
        TRUE ~ 4L
      )
    )

  if (score_by == "abs_cor") {
    module_summary <- module_summary %>%
      dplyr::arrange(class_rank, dplyr::desc(max_abs_cor), min_p, module)
  } else {
    module_summary <- module_summary %>%
      dplyr::arrange(class_rank, min_p, dplyr::desc(max_abs_cor), module)
  }

  ordered_modules <- module_summary$module
  module_order <- match(ordered_modules, current_levels)
  module_order <- module_order[!is.na(module_order)]

  list(
    module_order = module_order,
    ordered_modules = ordered_modules,
    module_summary = module_summary,
    exposure_traits = exposure_traits,
    outcome_traits = outcome_traits
  )
}

# ------------------------------------------------------------
# Exposure/outcome trait ordering
# (identical logic to script 09b; operates on a single reference
#  stats_df — here, dataset1)
# ------------------------------------------------------------
build_exposure_outcome_trait_order <- function(stats_df,
                                                found_traits,
                                                module_order_summary,
                                                p_thresh = 0.05,
                                                score_by = "p") {

  found_traits <- unique(as.character(found_traits))

  if (length(found_traits) == 0) {
    return(list(
      trait_order = integer(0),
      ordered_traits = character(0),
      trait_summary = data.frame()
    ))
  }

  if (is.null(module_order_summary) || nrow(module_order_summary) == 0) {
    trait_summary <- data.frame(
      trait = found_traits,
      trait_class = "input_order",
      min_p = NA_real_,
      max_abs_cor = NA_real_,
      n_sig_modules = NA_integer_,
      sig_exposure_and_outcome_modules = "",
      sig_exposure_only_modules = "",
      sig_outcome_only_modules = "",
      sig_other_modules = "",
      original_rank = seq_along(found_traits),
      class_rank = seq_along(found_traits),
      stringsAsFactors = FALSE
    )
    return(list(
      trait_order = seq_along(found_traits),
      ordered_traits = found_traits,
      trait_summary = trait_summary
    ))
  }

  cor_col <- get_cor_col(stats_df)

  module_class_df <- module_order_summary %>%
    dplyr::select(module, module_class)

  sig_df <- stats_df %>%
    dplyr::mutate(
      module = as.character(module),
      trait = as.character(trait),
      p = suppressWarnings(as.numeric(p)),
      cor_value = suppressWarnings(as.numeric(.data[[cor_col]]))
    ) %>%
    dplyr::filter(
      trait %in% found_traits,
      !is.na(p),
      p <= p_thresh
    ) %>%
    dplyr::left_join(module_class_df, by = "module") %>%
    dplyr::mutate(
      module_class = ifelse(is.na(module_class), "other", module_class)
    )

  trait_summary <- data.frame(
    trait = found_traits,
    original_rank = seq_along(found_traits),
    stringsAsFactors = FALSE
  )

  score_df <- sig_df %>%
    dplyr::group_by(trait) %>%
    dplyr::summarise(
      has_exposure_and_outcome = any(module_class == "exposure_and_outcome"),
      has_exposure_only = any(module_class == "exposure_only"),
      has_outcome_only = any(module_class == "outcome_only"),
      has_other = any(module_class == "other"),
      min_p = min(p, na.rm = TRUE),
      max_abs_cor = max(abs(cor_value), na.rm = TRUE),
      n_sig_modules = dplyr::n_distinct(module),
      sig_exposure_and_outcome_modules = paste(sort(unique(module[module_class == "exposure_and_outcome"])), collapse = ";"),
      sig_exposure_only_modules = paste(sort(unique(module[module_class == "exposure_only"])), collapse = ";"),
      sig_outcome_only_modules = paste(sort(unique(module[module_class == "outcome_only"])), collapse = ";"),
      sig_other_modules = paste(sort(unique(module[module_class == "other"])), collapse = ";"),
      .groups = "drop"
    )

  trait_summary <- trait_summary %>%
    dplyr::left_join(score_df, by = "trait") %>%
    dplyr::mutate(
      has_exposure_and_outcome = ifelse(is.na(has_exposure_and_outcome), FALSE, has_exposure_and_outcome),
      has_exposure_only = ifelse(is.na(has_exposure_only), FALSE, has_exposure_only),
      has_outcome_only = ifelse(is.na(has_outcome_only), FALSE, has_outcome_only),
      has_other = ifelse(is.na(has_other), FALSE, has_other),
      trait_class = dplyr::case_when(
        has_exposure_and_outcome ~ "exposure_and_outcome",
        has_exposure_only ~ "exposure_only",
        has_outcome_only ~ "outcome_only",
        TRUE ~ "other"
      ),
      min_p = ifelse(is.na(min_p), Inf, min_p),
      max_abs_cor = ifelse(is.na(max_abs_cor), 0, max_abs_cor),
      n_sig_modules = ifelse(is.na(n_sig_modules), 0, n_sig_modules),
      sig_exposure_and_outcome_modules = ifelse(is.na(sig_exposure_and_outcome_modules), "", sig_exposure_and_outcome_modules),
      sig_exposure_only_modules = ifelse(is.na(sig_exposure_only_modules), "", sig_exposure_only_modules),
      sig_outcome_only_modules = ifelse(is.na(sig_outcome_only_modules), "", sig_outcome_only_modules),
      sig_other_modules = ifelse(is.na(sig_other_modules), "", sig_other_modules),
      class_rank = dplyr::case_when(
        trait_class == "exposure_and_outcome" ~ 1L,
        trait_class == "exposure_only" ~ 2L,
        trait_class == "outcome_only" ~ 3L,
        TRUE ~ 4L
      )
    )

  if (score_by == "abs_cor") {
    trait_summary <- trait_summary %>%
      dplyr::arrange(class_rank, dplyr::desc(max_abs_cor), min_p, original_rank)
  } else {
    trait_summary <- trait_summary %>%
      dplyr::arrange(class_rank, min_p, dplyr::desc(max_abs_cor), original_rank)
  }

  ordered_traits <- trait_summary$trait
  trait_order <- match(ordered_traits, found_traits)
  trait_order <- trait_order[!is.na(trait_order)]

  list(
    trait_order = trait_order,
    ordered_traits = ordered_traits,
    trait_summary = trait_summary
  )
}

# ------------------------------------------------------------
# Dendro module ordering from consensus MEs (dataset1 reference)
# ------------------------------------------------------------
build_dendro_module_order <- function(MEs_ref, current_levels, module_dendro_distance = "bicor") {
  module_order <- seq_along(current_levels)

  MEs_ref <- as.matrix(MEs_ref)
  me_modules <- module_name_from_me_col(colnames(MEs_ref))
  keep <- !is_grey_module(colnames(MEs_ref)) & me_modules %in% current_levels
  MEs_ref <- MEs_ref[, keep, drop = FALSE]

  if (ncol(MEs_ref) < 2) return(module_order)

  moduleDendro <- getDendro(MEs_ref, distance = module_dendro_distance)
  dendro_modules <- module_name_from_me_col(colnames(MEs_ref)[moduleDendro$order])

  module_order <- match(intersect(dendro_modules, current_levels), current_levels)
  module_order <- module_order[!is.na(module_order)]

  remaining <- setdiff(seq_along(current_levels), module_order)
  module_order <- c(module_order, remaining)

  if (length(module_order) == 0) module_order <- seq_along(current_levels)
  module_order
}

# ------------------------------------------------------------
# Parse args
# ------------------------------------------------------------
option_list <- list(
  make_option("--project_root",           type = "character"),
  make_option("--consensus_modules_rds",  type = "character", default = NULL),
  make_option("--regions_file",           type = "character"),
  make_option("--adjustment_version",     type = "character", default = "unadjusted"),

  make_option("--dataset1_label",         type = "character"),
  make_option("--dataset1_stats_file",    type = "character"),
  make_option("--dataset2_label",         type = "character", default = NULL),
  make_option("--dataset2_stats_file",    type = "character", default = NULL),
  make_option("--dataset3_label",         type = "character", default = NULL),
  make_option("--dataset3_stats_file",    type = "character", default = NULL),

  make_option("--set_dir",                type = "character", default = NULL),
  make_option("--set_dir_recursive",      type = "logical",   default = FALSE,
              help = "If TRUE, search --set_dir recursively for .txt files [default = FALSE]"),
  make_option("--set_file",               type = "character", default = NULL),
  make_option("--set_file2",              type = "character", default = NULL),
  make_option("--set_file3",              type = "character", default = NULL),
  make_option("--set_file4",              type = "character", default = NULL),
  make_option("--set_file5",              type = "character", default = NULL),

  make_option("--module_dendro_distance", type = "character", default = "bicor",
              help = "bicor or pearson [default = bicor]"),

  make_option("--module_order_mode", type = "character", default = "dendro",
              help = "Module order mode: default, dendro, or exposure_outcome [default = dendro]"),

  make_option("--module_order_trait_file", type = "character", default = NULL,
              help = "TSV/CSV/TXT with columns category and trait. Required when --module_order_mode exposure_outcome."),

  make_option("--module_order_p_thresh", type = "double", default = NA,
              help = "P-value cutoff for exposure/outcome module classification. If not set, uses --p_thresh."),

  make_option("--module_order_score", type = "character", default = "p",
              help = "Within-class ordering score: p or abs_cor [default = p]"),

  make_option("--trait_order_mode", type = "character", default = "input",
              help = "Trait order mode inside each heatmap: input or exposure_outcome [default = input]"),

  make_option("--trait_order_p_thresh", type = "double", default = NA,
              help = "P-value cutoff for exposure/outcome trait ordering. If not set, uses --module_order_p_thresh."),

  make_option("--trait_order_score", type = "character", default = "p",
              help = "Within-class trait ordering score: p or abs_cor [default = p]"),

  make_option("--write_reordered_trait_files", type = "logical", default = TRUE,
              help = "Write a reordered trait text file per set under module_order/ [default = TRUE]"),

  make_option("--drop_grey_modules", type = "logical", default = TRUE,
              help = "Remove grey module from all datasets before plotting [default = TRUE]"),

  make_option("--p_thresh",               type = "double",    default = 0.05),
  make_option("--top_n",                  type = "integer",   default = 250),

  make_option("--min_sig_traits", type = "integer", default = 1,
              help = paste("Minimum number of significant trait associations (p < p_thresh)",
                           "a module must have (per dataset, within the current trait set)",
                           "to appear in that dataset's FOCUSED heatmap. [default = 1]")),

  make_option("--full_width",             type = "double",    default = 12),
  make_option("--full_height",            type = "double",    default = 12),
  make_option("--top_width",              type = "double",    default = 9),
  make_option("--top_height",             type = "double",    default = 6),

  make_option("--focused_width", type = "double", default = NULL,
              help = "FOCUSED heatmap width in inches. Defaults to --full_width if not set."),
  make_option("--focused_height", type = "double", default = 9,
              help = "FOCUSED heatmap height in inches [default = 9].")
)

opt <- parse_args(OptionParser(option_list = option_list))

# ------------------------------------------------------------
# Validate
# ------------------------------------------------------------
if (is.null(opt$project_root))        stop("--project_root is required")
if (is.null(opt$regions_file))        stop("--regions_file is required")
if (is.null(opt$dataset1_label))      stop("--dataset1_label is required")
if (is.null(opt$dataset1_stats_file)) stop("--dataset1_stats_file is required")

if (!dir.exists(opt$project_root)) stop("project_root not found: ", opt$project_root)
if (!file.exists(opt$regions_file)) stop("regions_file not found: ", opt$regions_file)
if (!file.exists(opt$dataset1_stats_file)) stop("dataset1_stats_file not found: ", opt$dataset1_stats_file)

module_dendro_distance <- tolower(opt$module_dendro_distance)
if (!module_dendro_distance %in% c("bicor", "pearson")) {
  stop("--module_dendro_distance must be bicor or pearson")
}

module_order_mode <- tolower(opt$module_order_mode)
if (!module_order_mode %in% c("default", "dendro", "exposure_outcome")) {
  stop("--module_order_mode must be one of: default, dendro, exposure_outcome")
}

module_order_score <- tolower(opt$module_order_score)
if (!module_order_score %in% c("p", "abs_cor")) {
  stop("--module_order_score must be one of: p, abs_cor")
}

trait_order_mode <- tolower(opt$trait_order_mode)
if (!trait_order_mode %in% c("input", "default", "exposure_outcome")) {
  stop("--trait_order_mode must be one of: input, exposure_outcome")
}
if (trait_order_mode == "default") trait_order_mode <- "input"

trait_order_score <- tolower(opt$trait_order_score)
if (!trait_order_score %in% c("p", "abs_cor")) {
  stop("--trait_order_score must be one of: p, abs_cor")
}

if (trait_order_mode == "exposure_outcome" && module_order_mode != "exposure_outcome") {
  stop("--trait_order_mode exposure_outcome requires --module_order_mode exposure_outcome")
}

if (module_order_mode == "exposure_outcome") {
  if (is.null(opt$module_order_trait_file) || !file.exists(opt$module_order_trait_file)) {
    stop("--module_order_trait_file is required and must exist when --module_order_mode exposure_outcome")
  }
}

if (opt$p_thresh <= 0 || opt$p_thresh >= 1) stop("--p_thresh must be > 0 and < 1")
if (opt$top_n < 1) stop("--top_n must be >= 1")
if (opt$min_sig_traits < 1) stop("--min_sig_traits must be >= 1")
if (!is.null(opt$focused_width) && opt$focused_width <= 0) stop("--focused_width must be > 0")
if (opt$focused_height <= 0) stop("--focused_height must be > 0")

module_order_p_thresh <- opt$module_order_p_thresh
if (is.na(module_order_p_thresh)) module_order_p_thresh <- opt$p_thresh
if (module_order_p_thresh <= 0 || module_order_p_thresh >= 1) {
  stop("--module_order_p_thresh must be > 0 and < 1")
}

trait_order_p_thresh <- opt$trait_order_p_thresh
if (is.na(trait_order_p_thresh)) trait_order_p_thresh <- module_order_p_thresh
if (trait_order_p_thresh <= 0 || trait_order_p_thresh >= 1) {
  stop("--trait_order_p_thresh must be > 0 and < 1")
}

dataset_inputs <- list(
  list(label = opt$dataset1_label, stats_file = opt$dataset1_stats_file)
)

for (i in 2:3) {
  lbl <- opt[[paste0("dataset", i, "_label")]]
  sf  <- opt[[paste0("dataset", i, "_stats_file")]]

  if (!is.null(lbl) || !is.null(sf)) {
    if (is.null(lbl) || is.null(sf)) {
      stop("Provide both --dataset", i, "_label and --dataset", i, "_stats_file")
    }
    if (!file.exists(sf)) {
      stop("dataset", i, "_stats_file not found: ", sf)
    }
    dataset_inputs[[length(dataset_inputs) + 1]] <- list(label = lbl, stats_file = sf)
  }
}

set_files <- collect_set_files(
  set_dir   = opt$set_dir,
  set_file  = opt$set_file,
  set_file2 = opt$set_file2,
  set_file3 = opt$set_file3,
  set_file4 = opt$set_file4,
  set_file5 = opt$set_file5,
  recursive = isTRUE(opt$set_dir_recursive)
)

WGCNA::enableWGCNAThreads()

# ------------------------------------------------------------
# Output dirs
# ------------------------------------------------------------

# ----------------------------------------------------------------
# Derive cpg_label / region_label from --regions_file path:
#   .../02_reference_region_filter/<ds>/<cpg_label>/<region_label>/Filtered_Regions.txt
# ----------------------------------------------------------------
{
  rfile_dir    <- dirname(opt$regions_file)
  region_label <- basename(rfile_dir)
  cpg_label    <- basename(dirname(rfile_dir))
  message("Derived cpg_label    : ", cpg_label)
  message("Derived region_label : ", region_label)
}
pipeline_root <- file.path(opt$project_root, "comethyl_output", "consensus")
step_dir      <- file.path(pipeline_root, "12b_me_trait_presentation")
cross_dir     <- file.path(step_dir, "cross_dataset", cpg_label, region_label, opt$adjustment_version)
order_dir     <- file.path(step_dir, "module_order", cpg_label, region_label, opt$adjustment_version)

dir.create(cross_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(order_dir, recursive = TRUE, showWarnings = FALSE)
message("Output root: ", step_dir)

# ------------------------------------------------------------
# Load all stats tables (drop grey modules here, before any ordering)
# ------------------------------------------------------------
stats_list <- list()

for (ds in dataset_inputs) {
  df <- readr::read_tsv(ds$stats_file, show_col_types = FALSE)
  df <- validate_stats_table(df, paste0(ds$label, " stats_file"))

  df$module <- as.character(df$module)
  df$trait  <- as.character(df$trait)
  df$p      <- suppressWarnings(as.numeric(df$p))

  if (isTRUE(opt$drop_grey_modules)) {
    n_before <- dplyr::n_distinct(df$module)
    df <- df %>% dplyr::filter(!is_grey_module(module))
    n_after <- dplyr::n_distinct(df$module)
    message("[", ds$label, "] Dropped grey module if present: ", n_before, " -> ", n_after, " module(s)")
  }

  stats_list[[ds$label]] <- df

  message(
    "Loaded stats for: ", ds$label,
    " (", nrow(df), " rows, ",
    length(unique(df$module)), " modules, ",
    length(unique(df$trait)), " traits)"
  )
}

# ------------------------------------------------------------
# Module ordering
#
# NOTE: Ordering is derived once, from dataset1 (the "reference"
# dataset), then reused for every dataset and the cross-dataset
# comparison. This keeps panels visually aligned. If a module
# present in another dataset is absent from dataset1's module set,
# it is appended at the end of that dataset's factor levels.
# ------------------------------------------------------------
first_label     <- dataset_inputs[[1]]$label
current_levels  <- unique(stats_list[[first_label]]$module)

module_order          <- NULL
ordered_modules        <- current_levels
module_order_summary  <- NULL

if (module_order_mode == "exposure_outcome") {

  order_res <- build_exposure_outcome_module_order(
    stats_df        = stats_list[[first_label]],
    trait_def_file  = opt$module_order_trait_file,
    p_thresh        = module_order_p_thresh,
    score_by        = module_order_score,
    current_levels  = current_levels
  )

  module_order         <- order_res$module_order
  module_order_summary <- order_res$module_summary
  ordered_modules       <- order_res$ordered_modules

  if (length(module_order) == 0) {
    warning("exposure_outcome module ordering returned zero modules; using default order.")
    module_order   <- NULL
    ordered_modules <- current_levels
  }

} else if (module_order_mode == "dendro") {

  if (!is.null(opt$consensus_modules_rds)) {
    if (!file.exists(opt$consensus_modules_rds)) {
      stop("consensus_modules_rds not found: ", opt$consensus_modules_rds)
    }

    consensusMods <- readRDS(opt$consensus_modules_rds)

    if (is.null(consensusMods$multiMEs)) {
      stop("Consensus object does not contain $multiMEs")
    }
    if (!first_label %in% names(consensusMods$multiMEs)) {
      stop(
        "dataset1 label '", first_label,
        "' not found in consensusMods$multiMEs. Available: ",
        paste(names(consensusMods$multiMEs), collapse = ", ")
      )
    }

    MEs_ref <- as.data.frame(consensusMods$multiMEs[[first_label]]$data)

    module_order <- build_dendro_module_order(
      MEs_ref = MEs_ref,
      current_levels = current_levels,
      module_dendro_distance = module_dendro_distance
    )
    ordered_modules <- current_levels[module_order]
    message("Module ordering derived from consensus MEs using dataset: ", first_label)
  } else {
    message("module_order_mode = dendro but --consensus_modules_rds not supplied; ",
            "falling back to default module order.")
    module_order   <- NULL
    ordered_modules <- current_levels
  }

} else {
  # "default": preserve original behavior of leaving module_order NULL
  module_order   <- NULL
  ordered_modules <- current_levels
}

write_vector_file(ordered_modules, file.path(order_dir, "ordered_modules.txt"))
if (!is.null(module_order_summary)) {
  readr::write_tsv(module_order_summary, file.path(order_dir, "module_order_exposure_outcome_summary.tsv"))
}

message("Module order mode: ", module_order_mode)
message("Trait order mode: ", trait_order_mode)

# ------------------------------------------------------------
# Per-dataset presentation heatmaps + comparison-ready outputs
# ------------------------------------------------------------
for (set_file in set_files) {

  set_name <- get_set_name(set_file)
  requested_traits <- read_trait_set_file(set_file)

  message("\n--- Trait set: ", set_name, " ---")

  set_order_dir <- file.path(order_dir, set_name)
  dir.create(set_order_dir, recursive = TRUE, showWarnings = FALSE)

  # ----------------------------------------------------------
  # Trait ordering (derived once from dataset1, reused everywhere)
  # ----------------------------------------------------------
  found_traits_ref <- intersect(requested_traits, unique(stats_list[[first_label]]$trait))
  all_available_traits <- unique(unlist(lapply(stats_list, function(df) unique(df$trait))))

  trait_order_summary <- NULL

  if (trait_order_mode == "exposure_outcome") {
    if (is.null(module_order_summary)) {
      warning("[", set_name, "] trait_order_mode=exposure_outcome requires module_order_mode=exposure_outcome. ",
              "Using input trait order.")
      ordered_traits_global <- requested_traits
    } else {
      trait_order_res <- build_exposure_outcome_trait_order(
        stats_df              = stats_list[[first_label]],
        found_traits           = found_traits_ref,
        module_order_summary  = module_order_summary,
        p_thresh               = trait_order_p_thresh,
        score_by                = trait_order_score
      )
      ordered_traits_ref   <- trait_order_res$ordered_traits
      trait_order_summary <- trait_order_res$trait_summary

      # Append any requested traits not found in dataset1 but available in
      # another dataset, preserving their original requested-file order.
      remainder <- intersect(
        requested_traits[!requested_traits %in% ordered_traits_ref],
        all_available_traits
      )
      ordered_traits_global <- c(ordered_traits_ref, remainder)

      message("[", set_name, "] trait_order_mode=exposure_outcome: ",
              "ordered ", length(ordered_traits_ref), " trait(s) from dataset1 (", first_label, "), ",
              length(remainder), " appended from other datasets.")
    }
  } else {
    # "input": keep requested order as-is
    ordered_traits_global <- requested_traits
  }

  write_vector_file(ordered_traits_global, file.path(set_order_dir, "traits_ordered.txt"))
  if (!is.null(trait_order_summary)) {
    readr::write_tsv(trait_order_summary, file.path(set_order_dir, "trait_order_exposure_outcome_summary.tsv"))
  }
  if (isTRUE(opt$write_reordered_trait_files)) {
    write_vector_file(ordered_traits_global, file.path(set_order_dir, paste0(set_name, "_ORDERED.txt")))
  }

  # ----------------------------------------------------------
  # Per-dataset outputs
  # ----------------------------------------------------------
  for (ds in dataset_inputs) {

    ds_label   <- ds$label
    ds_set_dir <- file.path(step_dir, ds_label, cpg_label, region_label, opt$adjustment_version, set_name)
    dir.create(ds_set_dir, recursive = TRUE, showWarnings = FALSE)

    stats_df <- stats_list[[ds_label]]
    available_traits <- unique(stats_df$trait)
    found_traits <- intersect(requested_traits, available_traits)
    missing_traits <- setdiff(requested_traits, available_traits)

    write_vector_file(
      requested_traits,
      file.path(ds_set_dir, paste0(ds_label, "_", set_name, "_traits_requested.txt"))
    )
    write_vector_file(
      found_traits,
      file.path(ds_set_dir, paste0(ds_label, "_", set_name, "_traits_found.txt"))
    )
    write_vector_file(
      missing_traits,
      file.path(ds_set_dir, paste0(ds_label, "_", set_name, "_traits_missing.txt"))
    )

    if (length(found_traits) == 0) {
      message(ds_label, "/", set_name, ": no traits found — skipping")
      next
    }

    # Order this dataset's found traits using the shared global trait order,
    # dropping any global-order traits absent from this dataset.
    ordered_traits_ds <- ordered_traits_global[ordered_traits_global %in% found_traits]
    ordered_traits_ds <- c(ordered_traits_ds, setdiff(found_traits, ordered_traits_ds))

    # Order this dataset's modules using the shared module order, appending
    # any modules present here but absent from dataset1's reference levels.
    ds_levels <- unique(as.character(stats_df$module))
    ds_module_levels <- ordered_modules[ordered_modules %in% ds_levels]
    ds_module_levels <- c(ds_module_levels, setdiff(ds_levels, ds_module_levels))
    ds_module_order  <- seq_along(ds_module_levels)

    subset_df <- dplyr::filter(stats_df, trait %in% found_traits)

    # Re-factor explicitly for plotMEtraitCor()
    subset_df$module <- factor(
      as.character(subset_df$module),
      levels = ds_module_levels
    )

    subset_df$trait <- factor(
      as.character(subset_df$trait),
      levels = ordered_traits_ds
    )

    subset_sig <- dplyr::filter(subset_df, !is.na(p), p < opt$p_thresh)

    # ----------------------------------------------------------
    # Focused module filtering (min_sig_traits)
    # ----------------------------------------------------------
    focused_modules <- subset_sig %>%
      dplyr::group_by(module) %>%
      dplyr::summarise(n_sig = dplyr::n_distinct(trait), .groups = "drop") %>%
      dplyr::filter(n_sig >= opt$min_sig_traits) %>%
      dplyr::pull(module) %>%
      as.character()

    focused_df <- subset_df %>%
      dplyr::filter(as.character(module) %in% focused_modules)

    focused_df$module <- factor(
      as.character(focused_df$module),
      levels = ds_module_levels[ds_module_levels %in% focused_modules]
    )
    focused_df$trait <- factor(
      as.character(focused_df$trait),
      levels = ordered_traits_ds
    )
    focused_module_order <- seq_along(levels(focused_df$module))

    message(ds_label, "/", set_name, ": FOCUSED heatmap — ",
            length(focused_modules), " of ", length(ds_module_levels),
            " modules have >= ", opt$min_sig_traits, " significant trait association(s)")

    write_vector_file(
      focused_modules,
      file.path(ds_set_dir, paste0(ds_label, "_", set_name, "_focused_modules.txt"))
    )

    readr::write_tsv(
      subset_df,
      file.path(ds_set_dir, paste0(ds_label, "_", set_name, "_subset_stats.tsv"))
    )
    readr::write_tsv(
      subset_sig,
      file.path(ds_set_dir, paste0(ds_label, "_", set_name, "_subset_stats_significant.tsv"))
    )

    full_pdf    <- file.path(ds_set_dir, paste0(ds_label, "_", set_name, "_ME_Trait_Heatmap_FULL.pdf"))
    top_pdf     <- file.path(ds_set_dir, paste0(ds_label, "_", set_name, "_ME_Trait_Heatmap_TOP.pdf"))
    focused_pdf <- file.path(ds_set_dir, paste0(ds_label, "_", set_name, "_ME_Trait_Heatmap_FOCUSED.pdf"))

    tryCatch({
      plotMEtraitCor(
        subset_df,
        moduleOrder     = ds_module_order,
        p               = opt$p_thresh,
        topOnly         = FALSE,
        file            = full_pdf,
        width           = opt$full_width,
        height          = opt$full_height,
        colColorMargins = c(-2.5, 4.21, 3.0, 12.07)
      )
      message(ds_label, "/", set_name, ": saved FULL heatmap")
    }, error = function(e) {
      message(ds_label, "/", set_name, ": FULL heatmap failed — ", conditionMessage(e))
    })

    tryCatch({
      plotMEtraitCor(
        subset_df,
        moduleOrder     = ds_module_order,
        p               = opt$p_thresh,
        topOnly         = TRUE,
        nTop            = opt$top_n,
        label.type      = "p",
        label.size      = 4,
        label.nudge_y   = 0,
        legend.position = c(1.11, 0.795),
        colColorMargins = c(-1, 4.75, 0.5, 10.1),
        file            = top_pdf,
        width           = opt$top_width,
        height          = opt$top_height
      )
      message(ds_label, "/", set_name, ": saved TOP heatmap")
    }, error = function(e) {
      message(ds_label, "/", set_name, ": TOP heatmap failed — ", conditionMessage(e))
    })

    if (length(focused_modules) == 0) {
      message(ds_label, "/", set_name, ": Skipping FOCUSED heatmap — no modules meet min_sig_traits >= ", opt$min_sig_traits)
      writeLines(
        paste("Skipped: no modules with >= ", opt$min_sig_traits, " significant trait associations"),
        con = file.path(ds_set_dir, paste0(ds_label, "_", set_name, "_ME_Trait_Heatmap_FOCUSED_SKIPPED.txt"))
      )
    } else {
      tryCatch({
        plotMEtraitCor(
          focused_df,
          moduleOrder     = focused_module_order,
          p               = opt$p_thresh,
          topOnly         = FALSE,
          file            = focused_pdf,
          width           = if (is.null(opt$focused_width)) opt$full_width else opt$focused_width,
          height          = opt$focused_height,
          colColorMargins = c(-2.5, 4.21, 3.0, 12.07)
        )
        message(ds_label, "/", set_name, ": saved FOCUSED heatmap")
      }, error = function(e) {
        message(ds_label, "/", set_name, ": FOCUSED heatmap failed — ", conditionMessage(e))
      })
    }

    write_log_lines(
      c(
        paste("dataset_label:", ds_label),
        paste("set_name:", set_name),
        paste("set_file:", set_file),
        paste("stats_file:", ds$stats_file),
        paste("drop_grey_modules:", opt$drop_grey_modules),
        paste("module_order_mode:", module_order_mode),
        paste("module_order_trait_file:", ifelse(is.null(opt$module_order_trait_file), "NULL", opt$module_order_trait_file)),
        paste("module_order_p_thresh:", module_order_p_thresh),
        paste("module_order_score:", module_order_score),
        paste("trait_order_mode:", trait_order_mode),
        paste("trait_order_p_thresh:", trait_order_p_thresh),
        paste("trait_order_score:", trait_order_score),
        paste("p_thresh:", opt$p_thresh),
        paste("top_n:", opt$top_n),
        paste("min_sig_traits:", opt$min_sig_traits),
        paste("n_focused_modules:", length(focused_modules)),
        paste("n_traits_requested:", length(requested_traits)),
        paste("n_traits_found:", length(found_traits)),
        paste("n_traits_missing:", length(missing_traits)),
        paste("n_rows_subset:", nrow(subset_df)),
        paste("n_rows_significant:", nrow(subset_sig)),
        paste("date:", as.character(Sys.time()))
      ),
      file.path(ds_set_dir, "run_parameters.txt")
    )
  }

  # ----------------------------------------------------------
  # Comparison-ready outputs in shared folder
  # ----------------------------------------------------------
  cross_set_dir <- file.path(cross_dir, set_name)
  dir.create(cross_set_dir, recursive = TRUE, showWarnings = FALSE)

  n_cross_saved <- 0

  for (ds in dataset_inputs) {
    ds_label <- ds$label
    stats_df <- stats_list[[ds_label]]

    found_traits <- intersect(requested_traits, unique(stats_df$trait))
    if (length(found_traits) == 0) {
      message("Cross/", set_name, "/", ds_label, ": no traits found — skipping")
      next
    }

    ordered_traits_ds <- ordered_traits_global[ordered_traits_global %in% found_traits]
    ordered_traits_ds <- c(ordered_traits_ds, setdiff(found_traits, ordered_traits_ds))

    ds_levels <- unique(as.character(stats_df$module))
    ds_module_levels <- ordered_modules[ordered_modules %in% ds_levels]
    ds_module_levels <- c(ds_module_levels, setdiff(ds_levels, ds_module_levels))
    ds_module_order  <- seq_along(ds_module_levels)

    subset_df <- dplyr::filter(stats_df, trait %in% found_traits)

    # Re-factor explicitly for plotMEtraitCor()
    subset_df$module <- factor(
      as.character(subset_df$module),
      levels = ds_module_levels
    )

    subset_df$trait <- factor(
      as.character(subset_df$trait),
      levels = ordered_traits_ds
    )

    subset_sig <- dplyr::filter(subset_df, !is.na(p), p < opt$p_thresh)

    focused_modules <- subset_sig %>%
      dplyr::group_by(module) %>%
      dplyr::summarise(n_sig = dplyr::n_distinct(trait), .groups = "drop") %>%
      dplyr::filter(n_sig >= opt$min_sig_traits) %>%
      dplyr::pull(module) %>%
      as.character()

    focused_df <- subset_df %>%
      dplyr::filter(as.character(module) %in% focused_modules)

    focused_df$module <- factor(
      as.character(focused_df$module),
      levels = ds_module_levels[ds_module_levels %in% focused_modules]
    )
    focused_df$trait <- factor(
      as.character(focused_df$trait),
      levels = ordered_traits_ds
    )
    focused_module_order <- seq_along(levels(focused_df$module))

    full_pdf    <- file.path(cross_set_dir, paste0(set_name, "_Comparison_", ds_label, "_FULL.pdf"))
    top_pdf     <- file.path(cross_set_dir, paste0(set_name, "_Comparison_", ds_label, "_TOP.pdf"))
    focused_pdf <- file.path(cross_set_dir, paste0(set_name, "_Comparison_", ds_label, "_FOCUSED.pdf"))

    tryCatch({
      plotMEtraitCor(
        subset_df,
        moduleOrder     = ds_module_order,
        p               = opt$p_thresh,
        topOnly         = FALSE,
        file            = full_pdf,
        width           = opt$full_width,
        height          = opt$full_height,
        colColorMargins = c(-2.5, 4.21, 3.0, 12.07)
      )
      message("Cross/", set_name, "/", ds_label, ": saved FULL")
    }, error = function(e) {
      message("Cross/", set_name, "/", ds_label, ": FULL failed — ", conditionMessage(e))
    })

    tryCatch({
      plotMEtraitCor(
        subset_df,
        moduleOrder     = ds_module_order,
        p               = opt$p_thresh,
        topOnly         = TRUE,
        nTop            = opt$top_n,
        label.type      = "p",
        label.size      = 4,
        label.nudge_y   = 0,
        legend.position = c(1.11, 0.795),
        colColorMargins = c(-1, 4.75, 0.5, 10.1),
        file            = top_pdf,
        width           = opt$top_width,
        height          = opt$top_height
      )
      message("Cross/", set_name, "/", ds_label, ": saved TOP")
    }, error = function(e) {
      message("Cross/", set_name, "/", ds_label, ": TOP failed — ", conditionMessage(e))
    })

    if (length(focused_modules) == 0) {
      message("Cross/", set_name, "/", ds_label, ": Skipping FOCUSED — no modules meet min_sig_traits >= ", opt$min_sig_traits)
      writeLines(
        paste("Skipped: no modules with >= ", opt$min_sig_traits, " significant trait associations"),
        con = file.path(cross_set_dir, paste0(set_name, "_Comparison_", ds_label, "_FOCUSED_SKIPPED.txt"))
      )
    } else {
      tryCatch({
        plotMEtraitCor(
          focused_df,
          moduleOrder     = focused_module_order,
          p               = opt$p_thresh,
          topOnly         = FALSE,
          file            = focused_pdf,
          width           = if (is.null(opt$focused_width)) opt$full_width else opt$focused_width,
          height          = opt$focused_height,
          colColorMargins = c(-2.5, 4.21, 3.0, 12.07)
        )
        message("Cross/", set_name, "/", ds_label, ": saved FOCUSED")
      }, error = function(e) {
        message("Cross/", set_name, "/", ds_label, ": FOCUSED failed — ", conditionMessage(e))
      })
    }

    n_cross_saved <- n_cross_saved + 1
  }

  write_log_lines(
    c(
      paste("set_name:", set_name),
      paste("set_file:", set_file),
      paste("datasets:", paste(sapply(dataset_inputs, `[[`, "label"), collapse = ", ")),
      paste("n_traits_requested:", length(requested_traits)),
      paste("n_datasets_saved:", n_cross_saved),
      paste("drop_grey_modules:", opt$drop_grey_modules),
      paste("module_order_mode:", module_order_mode),
      paste("trait_order_mode:", trait_order_mode),
      paste("p_thresh:", opt$p_thresh),
      paste("top_n:", opt$top_n),
      paste("min_sig_traits:", opt$min_sig_traits),
      paste("date:", as.character(Sys.time()))
    ),
    file.path(cross_set_dir, "run_parameters.txt")
  )
}

# ------------------------------------------------------------
# Top-level log
# ------------------------------------------------------------
write_log_lines(
  c(
    paste("project_root:", opt$project_root),
    paste("consensus_modules_rds:", ifelse(is.null(opt$consensus_modules_rds), "NULL", opt$consensus_modules_rds)),
    paste("adjustment_version:", opt$adjustment_version),
    paste("datasets:", paste(sapply(dataset_inputs, `[[`, "label"), collapse = ", ")),
    paste("reference_dataset_for_ordering:", first_label),
    paste("set_files:", paste(set_files, collapse = ", ")),
    paste("set_dir_recursive:", opt$set_dir_recursive),
    paste("module_dendro_distance:", module_dendro_distance),
    paste("module_order_mode:", module_order_mode),
    paste("module_order_trait_file:", ifelse(is.null(opt$module_order_trait_file), "NULL", opt$module_order_trait_file)),
    paste("module_order_p_thresh:", module_order_p_thresh),
    paste("module_order_score:", module_order_score),
    paste("trait_order_mode:", trait_order_mode),
    paste("trait_order_p_thresh:", trait_order_p_thresh),
    paste("trait_order_score:", trait_order_score),
    paste("write_reordered_trait_files:", opt$write_reordered_trait_files),
    paste("drop_grey_modules:", opt$drop_grey_modules),
    paste("p_thresh:", opt$p_thresh),
    paste("top_n:", opt$top_n),
    paste("min_sig_traits:", opt$min_sig_traits),
    paste("date:", as.character(Sys.time()))
  ),
  file.path(step_dir, "run_parameters.txt")
)

writeLines(
  capture.output(sessionInfo()),
  con = file.path(step_dir, "sessionInfo.txt")
)

message("\n\u2713 Script 12B complete: consensus ME-trait presentation heatmaps finished")
message("Outputs saved under: ", step_dir)