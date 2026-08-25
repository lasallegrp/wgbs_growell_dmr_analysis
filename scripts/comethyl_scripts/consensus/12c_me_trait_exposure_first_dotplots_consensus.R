#!/usr/bin/env Rscript

# ================================================================
# SCRIPT 12C: Consensus ME-Trait Category-First Presentation Dot Plots
#
# Pipeline: comethyl WGBS consensus analysis
#
# PURPOSE
#   Builds a single "<category>-robust" consensus module list (outcome
#   by default; can be switched to exposure via --selection_trait_category),
#   then renders three publication-style dot plots (bicor = fill color,
#   -log10(p) = dot size, dot drawn only where nominally significant)
#   sharing that SAME module list and row order:
#     1) EXPOSURES only
#     2) OUTCOMES only
#     3) COMBINED (outcomes block + exposures block, divided by a
#        vertical line, matching the presentation figure style)
#
#   Module selection logic:
#     A module is kept if it has >= --min_exposure_traits trait(s) from
#     the --selection_trait_category category that are nominally
#     significant (p < --p_thresh) in >= --min_timepoints_replicated of
#     the datasets provided, with (by default) a consistent bicor sign
#     across those significant timepoints. This mirrors the robustness
#     convention already used elsewhere in this project. The resulting
#     module list and row order are then reused, unchanged, for the
#     other two plots -- so a module can appear with a blank row in
#     whichever plot corresponds to the category it was NOT selected on.
#
# REQUIRED INPUTS
#   --project_root               : root directory of the analysis project
#   --regions_file                : path to Filtered_Regions.txt from script 02
#   --dataset1_label              : label for dataset 1
#   --dataset1_stats_file         : ME-trait stats TSV from script 12a
#   --exposure_outcome_trait_file : TSV/CSV/TXT with columns category, trait
#                                   (category must be "exposure" or "outcome"),
#                                   plus an OPTIONAL third column "timepoint".
#                                   Base 2-column format is the same file used
#                                   by 12b's --module_order_trait_file.
#                                   timepoint restricts a trait to a single
#                                   dataset's panel -- its value must exactly
#                                   match one of --dataset1_label/2/3_label
#                                   (e.g. "Baseline", "36-38wk", "Postpartum").
#                                   Leave blank/omit the column entirely for a
#                                   trait that should appear in every panel
#                                   (this is the old, still-default, behavior).
#                                   In the cross_dataset/ combined figure, a
#                                   restricted trait's column stays in every
#                                   panel (for alignment) but only shows a dot
#                                   in its own timepoint's panel. In the
#                                   per-timepoint plots, a restricted trait's
#                                   column is dropped entirely from any panel
#                                   it doesn't belong to.
#   --selection_trait_category    : exposure or outcome [default = outcome].
#                                   Controls which trait category drives module
#                                   selection AND row order. The other category
#                                   is still plotted in its own dot plot and in
#                                   the combined plot -- switching this only
#                                   changes which category decides which modules
#                                   are kept and how rows are ranked.
#
# OPTIONAL INPUTS
#   --dataset2_label / --dataset2_stats_file : dataset 2
#   --dataset3_label / --dataset3_stats_file : dataset 3
#   --adjustment_version         : label matching script 12a run [default = unadjusted]
#
#   --drop_grey_modules          : TRUE/FALSE [default = TRUE]
#   --p_thresh                   : nominal significance threshold [default = 0.05]
#   --min_timepoints_replicated  : # of datasets the SAME exposure trait must be
#                                   nominally significant in for a module to count
#                                   as exposure-robust [default = 2]
#   --require_consistent_direction : TRUE/FALSE; require the same bicor sign across
#                                   the replicated timepoints [default = TRUE]
#   --min_exposure_traits        : minimum number of robust exposure traits a module
#                                   must have to be kept [default = 1]
#   --drop_empty_columns         : TRUE/FALSE; drop trait columns with zero
#                                   significant dots anywhere in that plot [default = FALSE]
#
#   --color_limit                 : symmetric bicor color-scale limit. If not set,
#                                    auto-computed from the max |bicor| actually plotted.
#   --dot_size_min / --dot_size_max : point size range in mm [default = 1.5 / 6]
#
#   --font_family                 : font family passed to the plotting device and
#                                    theme [default = Helvetica]. If Helvetica is not
#                                    resolvable by your Cairo/font setup, rerun with
#                                    "sans" -- this script does not auto-detect fonts.
#   --font_size                   : base point size [default = 8, per npj Women's Health spec]
#
#   --full_width_mm               : figure width in mm [default = 170, npj full-width]
#   --row_height_in               : inches of height per module row, used to
#                                    auto-compute figure height [default = 0.155]
#   --header_margin_in            : inches reserved above the plot area for
#                                    title/legend/group labels [default = 1.0]
#   --max_height_in               : height (in) above which a warning is printed
#                                    suggesting the figure won't fit a single
#                                    journal page [default = 9.5]
#   --exposure_plot_width/height/--outcome_plot_width/height/
#   --combined_plot_width/height  : inches; override the auto-computed size for
#                                    that specific plot. [default = NULL = auto]
#
#   --exposure_group_label        : label text for the exposures column block
#                                    [default = "Exposures"]
#   --outcome_group_label         : label text for the outcomes column block
#                                    [default = "Outcomes"]
#   --add_panel_letters           : TRUE/FALSE; prefix facet strips with a), b), c)
#                                    [default = TRUE]
#
#   --save_preview_png            : TRUE/FALSE; also save a PNG for slides
#                                    [default = TRUE]
#   --preview_dpi                 : DPI for the PNG preview [default = 200]
#
# OUTPUTS
#   comethyl_output/consensus/12c_me_trait_exposure_first_dotplots/
#       <cpg_label>/<region_label>/<adjustment_version>/
#           <selection_trait_category>_robust_modules.txt
#           <selection_trait_category>_robust_module_summary.tsv
#           exposure_outcome_pair_summary.tsv   (all module x trait pairs, robust or not)
#           run_parameters.txt
#           sessionInfo.txt
#
#           cross_dataset/   <- SHARED/COMBINED: all timepoints in one figure
#               Dotplot_Exposures.pdf  (+ .png if --save_preview_png)
#               Dotplot_Outcomes.pdf   (+ .png)
#               Dotplot_Combined.pdf   (+ .png)
#               Dotplot_Exposures_data.tsv   (rows actually plotted, i.e. nominally significant)
#               Dotplot_Outcomes_data.tsv
#               Dotplot_Combined_data.tsv
#
#           <dataset1_label>/   <- PER-TIMEPOINT: single-panel versions, one folder per dataset
#               Dotplot_Exposures_<dataset1_label>.pdf  (+ .png, + _data.tsv)
#               Dotplot_Outcomes_<dataset1_label>.pdf
#               Dotplot_Combined_<dataset1_label>.pdf
#           <dataset2_label>/  (same three files, if dataset2 provided)
#           <dataset3_label>/  (same three files, if dataset3 provided)
#
#   Skipped plots (no nominally significant rows for that module/trait/
#   dataset combination) get a "..._SKIPPED.txt" note instead of a PDF.
#
# NOTES
#   - cpg_label and region_label are derived from --regions_file path.
#   - Vector PDFs are written with grDevices::cairo_pdf, matching the
#     project's Illustrator-assembly convention (panels placed as
#     separate vector PDFs). Convert to 600 DPI TIFF in Illustrator/
#     Photoshop at final placement time, per npj Women's Health specs.
#   - Figure dimensions are auto-computed from the module count, trait
#     column count, and --panel_layout. This script prints the computed
#     size and warns when it won't fit journal full-width/a single page --
#     check that message before assuming a figure is presentation- or
#     print-ready.
#   - --panel_layout vertical stacks timepoint panels top-to-bottom instead
#     of side-by-side, trading width for height -- useful when a wide
#     trait-column count won't fit n_timepoints panels at journal
#     full-width. Per-timepoint plots (single dataset each) are unaffected
#     by --panel_layout since they only ever have one panel.
#   - --exposure_plot_width/height, --outcome_plot_width/height, and
#     --combined_plot_width/height override the auto-computed size for the
#     cross_dataset/ (combined) version ONLY. Per-timepoint plots are
#     always auto-sized.
#
# EXAMPLE
#   Rscript $SCRIPTS/consensus/12c_me_trait_exposure_first_dotplots_consensus.R \
#     --project_root $PROJECT \
#     --regions_file $OUT/02_reference_region_filter/Baseline/$CPG_LABEL/$REGION_LABEL/Filtered_Regions.txt \
#     --adjustment_version v2_exclude_protected_pcs \
#     --dataset1_label Baseline \
#     --dataset1_stats_file $OUT/12a_me_trait_analysis/Baseline/$CPG_LABEL/$REGION_LABEL/v2_exclude_protected_pcs/ME_Trait_Correlation_Stats_Bicor.tsv \
#     --dataset2_label 36-38wk \
#     --dataset2_stats_file $OUT/12a_me_trait_analysis/36-38wk/$CPG_LABEL/$REGION_LABEL/v2_exclude_protected_pcs/ME_Trait_Correlation_Stats_Bicor.tsv \
#     --dataset3_label Postpartum \
#     --dataset3_stats_file $OUT/12a_me_trait_analysis/Postpartum/$CPG_LABEL/$REGION_LABEL/v2_exclude_protected_pcs/ME_Trait_Correlation_Stats_Bicor.tsv \
#     --exposure_outcome_trait_file $CONFIG/exposure_outcome_traits.tsv
# ================================================================

message("Starting Script 12c \u2713")

suppressPackageStartupMessages({
  library(optparse)
  library(readr)
  library(dplyr)
  library(ggplot2)
  library(scales)
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

required_helpers <- c("write_vector_file", "write_log_lines")
missing_helpers <- required_helpers[!vapply(required_helpers, exists, logical(1), mode = "function")]
if (length(missing_helpers) > 0) {
  stop("Missing helper functions in helper.R: ", paste(missing_helpers, collapse = ", "))
}

# ggplot2's internal mm<->pt conversion constant, reimplemented locally
# so this script does not depend on an unexported ggplot2 internal.
PT_TO_MM <- 25.4 / 72.27

# ------------------------------------------------------------
# Helpers
# ------------------------------------------------------------
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

is_grey_module <- function(x) {
  tolower(sub("^ME", "", as.character(x))) == "grey"
}

read_exposure_outcome_file <- function(file) {
  if (is.null(file) || !nzchar(file)) {
    stop("--exposure_outcome_trait_file is required")
  }
  if (!file.exists(file)) stop("exposure_outcome_trait_file not found: ", file)

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
    stop("exposure_outcome_trait_file missing required columns: ",
         paste(missing, collapse = ", "))
  }

  x$category <- tolower(trimws(as.character(x$category)))
  x$trait <- trimws(as.character(x$trait))

  # Optional third column: restricts a trait to a single dataset/timepoint's
  # panel (value must exactly match one of --dataset1_label/2/3_label). Blank
  # or absent = trait applies to every panel, which preserves old 2-column
  # file behavior unchanged.
  if ("timepoint" %in% colnames(x)) {
    x$timepoint <- trimws(as.character(x$timepoint))
    x$timepoint[x$timepoint == ""] <- NA_character_
  } else {
    x$timepoint <- NA_character_
  }

  x <- x[x$category %in% c("exposure", "outcome") &
         !is.na(x$trait) & x$trait != "", , drop = FALSE]

  if (nrow(x) == 0) {
    stop("exposure_outcome_trait_file has no valid exposure/outcome rows.")
  }

  x <- x[!duplicated(x[, c("category", "trait")]), ]
  x[, c("category", "trait", "timepoint")]
}

# ------------------------------------------------------------
# Exposure-first robust module/trait pair computation
# ------------------------------------------------------------
compute_pair_summary <- function(stats_list, dataset_labels, candidate_traits,
                                  p_thresh, min_timepoints, require_consistent_direction) {

  long <- dplyr::bind_rows(lapply(dataset_labels, function(ds) {
    df <- stats_list[[ds]]
    df <- df[df$trait %in% candidate_traits, c("module", "trait", "p", "r_value")]
    if (nrow(df) == 0) return(df)
    df$dataset <- ds
    df
  }))

  if (nrow(long) == 0) {
    return(data.frame(
      module = character(0), trait = character(0), n_sig = integer(0),
      consistent_direction = logical(0), timepoints_sig = character(0),
      min_p = numeric(0), max_abs_r = numeric(0), robust = logical(0),
      stringsAsFactors = FALSE
    ))
  }

  long$sig <- !is.na(long$p) & long$p < p_thresh

  pair_summary <- long %>%
    dplyr::group_by(module, trait) %>%
    dplyr::summarise(
      n_sig = sum(sig, na.rm = TRUE),
      n_unique_signs = dplyr::n_distinct(sign(r_value[sig])),
      timepoints_sig = paste(sort(dataset[sig]), collapse = ";"),
      min_p = suppressWarnings(min(p[sig], na.rm = TRUE)),
      max_abs_r = suppressWarnings(max(abs(r_value[sig]), na.rm = TRUE)),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      consistent_direction = n_unique_signs <= 1,
      min_p = ifelse(is.infinite(min_p), NA_real_, min_p),
      max_abs_r = ifelse(is.infinite(max_abs_r) | is.nan(max_abs_r), NA_real_, max_abs_r),
      robust = n_sig >= min_timepoints & (!require_consistent_direction | consistent_direction)
    )

  as.data.frame(pair_summary)
}

# ------------------------------------------------------------
# Build the long data frame actually plotted (significant rows only,
# with full factor levels retained for consistent column/row grids)
# ------------------------------------------------------------
# ------------------------------------------------------------
# Build the long data frame actually plotted (significant rows only,
# with full factor levels retained for consistent column/row grids).
# trait_timepoint_restrict, if given, is a named character vector
# (trait -> dataset label) that blanks out rows for a trait in any
# dataset other than its assigned one, WITHOUT removing that trait's
# column from the grid -- used for the multi-panel combined view so
# columns stay aligned across panels while cross-timepoint spillover
# (e.g. a wk36-only exposure showing a dot in the Baseline panel) is
# suppressed. Traits absent from the vector, or with an NA entry, are
# unrestricted and behave exactly as before.
# ------------------------------------------------------------
build_plot_df <- function(stats_list, dataset_labels, modules_ordered, traits_ordered,
                           p_thresh, trait_timepoint_restrict = NULL) {
  long <- dplyr::bind_rows(lapply(dataset_labels, function(ds) {
    df <- stats_list[[ds]]
    df <- df[df$module %in% modules_ordered & df$trait %in% traits_ordered, ]
    if (nrow(df) == 0) return(df)
    df$dataset <- ds
    df
  }))

  if (!is.null(trait_timepoint_restrict) && nrow(long) > 0) {
    restrict_val <- trait_timepoint_restrict[long$trait]
    keep <- is.na(restrict_val) | (restrict_val == long$dataset)
    keep[is.na(keep)] <- TRUE # trait not present in the restriction map at all
    long <- long[keep, ]
  }

  long <- long[!is.na(long$p) & long$p < p_thresh, ]
  long$neglogp <- -log10(pmax(long$p, 1e-300))
  long$module <- factor(long$module, levels = modules_ordered)
  long$trait <- factor(long$trait, levels = traits_ordered)
  long$dataset <- factor(long$dataset, levels = dataset_labels)
  long
}

# ------------------------------------------------------------
# Per-panel variant for "compact" combined plots: each dataset gets its OWN
# trait column set (per_panel_traits, a named list keyed by dataset label),
# rather than one column set shared/blanked across every panel. Returns
# `trait` as plain character (not a shared factor) since each panel's
# columns differ -- ggplot facet_wrap(scales="free_x") derives each panel's
# own axis from what's actually present in its own data subset.
# ------------------------------------------------------------
build_plot_df_per_panel <- function(stats_list, dataset_labels, modules_ordered,
                                     per_panel_traits, p_thresh) {
  long <- dplyr::bind_rows(lapply(dataset_labels, function(ds) {
    traits_here <- per_panel_traits[[ds]]
    df <- stats_list[[ds]]
    df <- df[df$module %in% modules_ordered & df$trait %in% traits_here, ]
    df <- df[!is.na(df$p) & df$p < p_thresh, ]
    if (nrow(df) == 0) return(df)
    df$dataset <- ds
    df$neglogp <- -log10(pmax(df$p, 1e-300))
    df$module <- factor(df$module, levels = modules_ordered)
    # Explicit factor levels in the intended (outcomes-block, then
    # exposures-block) order -- without this, ggplot's discrete x-scale
    # defaults to ALPHABETICAL ordering of the raw trait strings, which
    # interleaves outcomes and exposures instead of keeping them in their
    # intended blocks, and desyncs the group-label/divider positions (which
    # assume the intended order) from where the columns actually land.
    df$trait <- factor(df$trait, levels = traits_here)
    df
  }))
  # Explicit factor levels in dataset-argument order -- without this,
  # facet_wrap sorts panels alphabetically by the "dataset" column's raw
  # values (e.g. "36-38wk" before "Baseline"), silently reordering panels
  # away from the intended a/b/c = dataset1/2/3 sequence.
  if (nrow(long) > 0) long$dataset <- factor(long$dataset, levels = dataset_labels)
  long
}

# ------------------------------------------------------------
# Dimension recommendation
# ------------------------------------------------------------
GROUP_LABEL_EXTRA_ROWS <- 3.2 # must match the y-scale expansion in make_dot_plot()
                               # (add = c(0.6, 2.6)); budgeted here as real
                               # height so it doesn't compress the module rows.

# Empirically-determined floor: ggplot's fixed-size chrome (title, subtitle,
# top margin) costs roughly the same ABSOLUTE amount of space regardless of
# how few module rows there are, AND the rotated x-axis trait labels at the
# BOTTOM of EACH PANEL need vertical space that scales with the longest
# trait name's character count -- a plot with names like
# "HEI2015_TOTAL_SCORE_base" (25 chars) needs meaningfully more headroom
# than one with names like "apo" or "SPL_SVM". Critically, that bottom-label
# footer is needed ONCE for horizontal layout (panels share one row, so one
# shared bottom edge) but ONCE PER STACKED PANEL for vertical layout (each
# panel has its own bottom edge) -- confirmed by direct visual testing this
# mattered: a formula that only budgeted one shared footer left vertical
# layouts with badly overlapping row labels even though the same footer
# budget was clean for horizontal. TOP_CHROME_IN/BOTTOM_LABEL_*_IN below are
# calibrated against three data points: short names/n_modules=2/horizontal
# (2.3in clean), 25-char names/n_modules=2/horizontal (~2.8-3.3in clean),
# and 25-char names/n_modules=2/vertical/3 panels (~9in clean) -- each with
# a safety margin above the tested-clean threshold.
MIN_FIGURE_HEIGHT_IN <- 2.3
TOP_CHROME_IN <- 1.0
BOTTOM_LABEL_BASE_IN <- 0.35
BOTTOM_LABEL_PER_CHAR_IN <- 0.065

bottom_label_footer_in <- function(max_trait_nchar, font_size) {
  BOTTOM_LABEL_BASE_IN + BOTTOM_LABEL_PER_CHAR_IN * (font_size / 8) * max_trait_nchar
}

# horizontal: panels sit side by side -> height is ONE panel's worth of rows
#             + ONE shared bottom-label footer (panels share a bottom edge).
# vertical:   panels stack top to bottom -> height is N panels' worth of
#             rows + N panels' worth of bottom-label footers (each stacked
#             panel has its own bottom edge), plus spacing between panels.
compute_recommended_height <- function(n_modules, row_height_in, header_margin_in,
                                        max_height_in, n_facets, panel_layout,
                                        panel_spacing_in, max_trait_nchar = 0, font_size = 8) {
  effective_top_in <- max(header_margin_in, TOP_CHROME_IN)
  bottom_in <- bottom_label_footer_in(max_trait_nchar, font_size)

  panel_multiplier <- if (panel_layout == "vertical") max(n_facets, 1) else 1
  spacing_total_in <- if (panel_layout == "vertical" && n_facets > 1) (n_facets - 1) * panel_spacing_in else 0

  height_in <- effective_top_in +
    (panel_multiplier * bottom_in) +
    (panel_multiplier * (n_modules + GROUP_LABEL_EXTRA_ROWS) * row_height_in) +
    spacing_total_in

  height_in <- max(height_in, MIN_FIGURE_HEIGHT_IN)

  list(
    height_in = height_in,
    exceeds_page = height_in > max_height_in
  )
}

# horizontal: panels sit side by side -> width is the SUM of each panel's
#             own column count (in "aligned" mode every panel has the same
#             count, so this equals n_traits * n_facets as before; in
#             "compact" mode panels can have different counts).
# vertical:   panels stack top to bottom -> width is the WIDEST single
#             panel's column count (this is the whole point of switching to
#             vertical -- trading width for height when there are too many
#             trait columns to fit n_timepoints panels side by side at
#             journal full-width).
#
# Each panel is floored at MIN_PANEL_WIDTH_IN regardless of how few columns
# it has -- a panel with only 2-3 applicable columns (e.g. after
# --drop_empty_columns trims a timepoint down small in "compact" mode)
# still needs enough width for its own "Outcomes"/"Exposures"
# group-block label text, or those labels overlap each other. This mirrors
# MIN_FIGURE_HEIGHT_IN's role for the row-label-overlap case.
MIN_PANEL_WIDTH_IN <- 1.8

compute_recommended_width <- function(n_traits_per_panel, col_width_in,
                                       width_margin_in, journal_width_in, panel_layout) {
  n_traits_per_panel <- pmax(n_traits_per_panel, 0)
  panel_widths_in <- pmax(n_traits_per_panel * col_width_in, MIN_PANEL_WIDTH_IN)
  panel_total <- if (panel_layout == "vertical") {
    if (length(panel_widths_in) == 0) 0 else max(panel_widths_in)
  } else {
    sum(panel_widths_in)
  }
  comfortable_width_in <- width_margin_in + panel_total
  list(
    comfortable_width_in = comfortable_width_in,
    journal_width_in     = journal_width_in,
    width_in             = max(comfortable_width_in, journal_width_in),
    exceeds_journal_width = comfortable_width_in > journal_width_in
  )
}

# ------------------------------------------------------------
# Dot plot builder
# ------------------------------------------------------------
make_dot_plot <- function(plot_df, modules_ordered, traits_ordered, dataset_labels,
                           group_blocks = NULL, color_limit = NA,
                           size_range = c(1.5, 6), font_family = "Helvetica",
                           font_size = 8, title = NULL, subtitle = NULL,
                           add_panel_letters = TRUE, panel_layout = "horizontal",
                           panel_spacing_in = 0.15, column_mode = "aligned",
                           grid_registration_df = NULL, per_panel_group_blocks = NULL) {

  n_mod <- length(modules_ordered)

  if (is.na(color_limit) || is.null(color_limit)) {
    color_limit <- suppressWarnings(max(abs(plot_df$r_value), na.rm = TRUE))
    if (!is.finite(color_limit) || color_limit <= 0) color_limit <- 0.6
  }

  facet_labels <- dataset_labels
  if (isTRUE(add_panel_letters)) {
    facet_labels <- paste0(letters[seq_along(dataset_labels)], ") ", dataset_labels)
  }
  names(facet_labels) <- dataset_labels

  # "compact" mode: each panel gets its own independent x-axis (free_x),
  # showing only the columns that actually apply to that timepoint --
  # avoids the "aligned" mode's wall of empty columns when many traits are
  # timepoint-restricted. "aligned" mode keeps one shared column layout
  # across every panel (useful when you want a trait's position to be
  # directly comparable panel to panel).
  facet_scales <- if (column_mode == "compact") "free_x" else "fixed"
  facet_layer <- if (panel_layout == "vertical") {
    facet_wrap(~dataset, ncol = 1, labeller = as_labeller(facet_labels), scales = facet_scales)
  } else {
    facet_wrap(~dataset, nrow = 1, labeller = as_labeller(facet_labels), scales = facet_scales)
  }

  p <- ggplot(plot_df, aes(x = trait, y = module)) +
    geom_point(aes(fill = r_value, size = neglogp), shape = 21, color = "black", stroke = 0.3) +
    scale_fill_gradient2(
      low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0,
      limits = c(-color_limit, color_limit), oob = scales::squish,
      name = "bicor (r)"
    ) +
    scale_size_continuous(
      range = size_range,
      name = expression(-log[10](p)),
      guide = guide_legend(override.aes = list(fill = "grey40"))
    )

  if (column_mode == "compact") {
    # No shared x-scale limits -- each panel derives its own columns from
    # whatever is present in its own data (real dots + the geom_blank
    # registration rows below, which register empty columns without
    # drawing anything visible).
    p <- p + scale_x_discrete(drop = TRUE)
    if (!is.null(grid_registration_df) && nrow(grid_registration_df) > 0) {
      p <- p + geom_blank(data = grid_registration_df, aes(x = trait, y = module), inherit.aes = FALSE)
    }
  } else {
    p <- p + scale_x_discrete(limits = traits_ordered, drop = FALSE)
  }

  p <- p +
    scale_y_discrete(
      limits = rev(modules_ordered), drop = FALSE,
      # Reserve fixed headroom (in row-units, independent of n_modules) above
      # the top row so the group-block labels never collide with the facet
      # strip text, regardless of how many modules are plotted.
      expand = expansion(add = c(0.6, 2.6))
    ) +
    facet_layer +
    labs(x = NULL, y = NULL, title = title, subtitle = subtitle) +
    theme_minimal(base_family = font_family, base_size = font_size) +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = font_size - 1),
      axis.text.y = element_text(size = font_size - 1),
      strip.text = element_text(face = "bold", size = font_size + 1),
      panel.grid.major = element_line(color = "grey85", linewidth = 0.2),
      panel.grid.minor = element_blank(),
      plot.title = element_text(size = font_size + 2, hjust = 0),
      plot.subtitle = element_text(size = font_size - 1),
      legend.title = element_text(size = font_size),
      legend.text = element_text(size = font_size - 1),
      plot.margin = margin(t = 30, r = 5, b = 5, l = 5, unit = "pt"),
      panel.spacing = grid::unit(panel_spacing_in, "in")
    ) +
    coord_cartesian(clip = "off")

  if (column_mode == "compact") {
    # Group-block labels/dividers differ per panel (each panel has a
    # different exposure column count), so they're built as per-panel data
    # frames and drawn as faceted layers instead of one annotate() call
    # replicated identically across every panel.
    if (!is.null(per_panel_group_blocks) && length(per_panel_group_blocks) > 0) {
      label_rows <- list()
      vline_rows <- list()
      for (ds in names(per_panel_group_blocks)) {
        blocks <- per_panel_group_blocks[[ds]]
        for (i in seq_along(blocks)) {
          blk <- blocks[[i]]
          label_rows[[length(label_rows) + 1]] <- data.frame(
            dataset = ds, x = mean(c(blk$start, blk$end)), y = n_mod + 1.4,
            label = blk$label, stringsAsFactors = FALSE
          )
          if (i < length(blocks)) {
            vline_rows[[length(vline_rows) + 1]] <- data.frame(
              dataset = ds, x = blk$end + 0.5, stringsAsFactors = FALSE
            )
          }
        }
      }
      label_df <- dplyr::bind_rows(label_rows)
      label_df$dataset <- factor(label_df$dataset, levels = dataset_labels)
      p <- p + geom_text(
        data = label_df, aes(x = x, y = y, label = label), inherit.aes = FALSE,
        fontface = "bold", size = (font_size - 0.5) * PT_TO_MM, hjust = 0.5
      )
      if (length(vline_rows) > 0) {
        vline_df <- dplyr::bind_rows(vline_rows)
        vline_df$dataset <- factor(vline_df$dataset, levels = dataset_labels)
        p <- p + geom_vline(data = vline_df, aes(xintercept = x),
                             color = "black", linewidth = 0.6)
      }
    }
  } else if (!is.null(group_blocks) && length(group_blocks) > 0) {
    for (i in seq_along(group_blocks)) {
      blk <- group_blocks[[i]]
      p <- p + annotate(
        "text",
        x = mean(c(blk$start, blk$end)),
        y = n_mod + 1.4,
        label = blk$label,
        fontface = "bold",
        size = (font_size - 0.5) * PT_TO_MM,
        hjust = 0.5
      )
      if (i < length(group_blocks)) {
        p <- p + geom_vline(xintercept = blk$end + 0.5, color = "black", linewidth = 0.6)
      }
    }
  }

  p
}

save_plot <- function(p, out_prefix, width_in, height_in, font_family,
                       save_png, preview_dpi) {
  pdf_file <- paste0(out_prefix, ".pdf")
  grDevices::cairo_pdf(pdf_file, width = width_in, height = height_in, family = font_family)
  print(p)
  grDevices::dev.off()

  if (isTRUE(save_png)) {
    png_file <- paste0(out_prefix, ".png")
    ggplot2::ggsave(
      png_file, plot = p, width = width_in, height = height_in,
      dpi = preview_dpi, units = "in", bg = "white"
    )
  }
  invisible(pdf_file)
}

# ------------------------------------------------------------
# Parse args
# ------------------------------------------------------------
option_list <- list(
  make_option("--project_root",           type = "character"),
  make_option("--regions_file",           type = "character"),
  make_option("--adjustment_version",     type = "character", default = "unadjusted"),

  make_option("--dataset1_label",         type = "character"),
  make_option("--dataset1_stats_file",    type = "character"),
  make_option("--dataset2_label",         type = "character", default = NULL),
  make_option("--dataset2_stats_file",    type = "character", default = NULL),
  make_option("--dataset3_label",         type = "character", default = NULL),
  make_option("--dataset3_stats_file",    type = "character", default = NULL),

  make_option("--exposure_outcome_trait_file", type = "character"),

  make_option("--selection_trait_category", type = "character", default = "outcome",
              help = paste("Which trait category drives module selection AND row order:",
                           "exposure or outcome [default = outcome]. The other category's",
                           "traits are still plotted in their own dot plot and in the combined",
                           "plot -- this only controls which category's significance decides",
                           "which modules are kept and how rows are ranked.")),

  make_option("--drop_grey_modules", type = "logical", default = TRUE),
  make_option("--p_thresh", type = "double", default = 0.05),
  make_option("--min_timepoints_replicated", type = "integer", default = 2),
  make_option("--require_consistent_direction", type = "logical", default = TRUE),
  make_option("--min_exposure_traits", type = "integer", default = 1),
  make_option("--drop_empty_columns", type = "logical", default = FALSE),

  make_option("--color_limit", type = "double", default = NA),
  make_option("--dot_size_min", type = "double", default = 1.5),
  make_option("--dot_size_max", type = "double", default = 6),

  make_option("--font_family", type = "character", default = "Helvetica"),
  make_option("--font_size", type = "double", default = 8),

  make_option("--full_width_mm", type = "double", default = 170),
  make_option("--row_height_in", type = "double", default = 0.155),
  make_option("--header_margin_in", type = "double", default = 1.0),
  make_option("--max_height_in", type = "double", default = 9.5),
  make_option("--col_width_in", type = "double", default = 0.16,
              help = paste("Inches of width reserved per trait column per facet panel,",
                           "used to auto-widen the figure beyond --full_width_mm when",
                           "there are enough columns that the journal full-width would",
                           "make columns/group-labels illegible or overlapping. [default = 0.16]")),
  make_option("--width_margin_in", type = "double", default = 1.3,
              help = "Inches reserved for y-axis module labels + legend, outside the facet panels [default = 1.3]"),

  make_option("--exposure_plot_width", type = "double", default = NA),
  make_option("--exposure_plot_height", type = "double", default = NA),
  make_option("--outcome_plot_width", type = "double", default = NA),
  make_option("--outcome_plot_height", type = "double", default = NA),
  make_option("--combined_plot_width", type = "double", default = NA),
  make_option("--combined_plot_height", type = "double", default = NA),

  make_option("--exposure_group_label", type = "character", default = "Exposures"),
  make_option("--outcome_group_label", type = "character", default = "Outcomes"),
  make_option("--add_panel_letters", type = "logical", default = TRUE),

  make_option("--panel_layout", type = "character", default = "horizontal",
              help = paste("horizontal or vertical [default = horizontal]. horizontal places",
                           "timepoint panels side by side (a, b, c across); vertical stacks",
                           "them (a, b, c down). vertical trades width for height -- useful",
                           "when a wide trait-column count won't fit journal full-width, since",
                           "each panel only needs ONE panel's worth of width instead of",
                           "n_timepoints worth.")),
  make_option("--per_timepoint_plots", type = "logical", default = TRUE,
              help = paste("TRUE/FALSE; also render single-panel (one timepoint each) versions",
                           "of all three plots into per-timepoint folders, in addition to the",
                           "multi-panel combined version in cross_dataset/ [default = TRUE]")),
  make_option("--panel_spacing_in", type = "double", default = NA,
              help = paste("Inches of spacing between stacked/side-by-side panels.",
                           "If not set, defaults to 0.15in for horizontal layout and",
                           "0.55in for vertical layout (vertical needs more room for",
                           "the next panel's facet strip + group labels).")),
  make_option("--combined_column_mode", type = "character", default = "compact",
              help = paste("aligned or compact [default = compact]. Affects ONLY the",
                           "cross_dataset/ (combined, multi-panel) figures. 'aligned' shows",
                           "the same trait columns in every panel, blanking (not removing)",
                           "any dot outside a trait's assigned timepoint -- keeps a trait's",
                           "position directly comparable across panels, but gets cluttered",
                           "fast when many traits are timepoint-restricted. 'compact' gives",
                           "each panel its own independent set of columns, showing only the",
                           "traits that actually apply to that timepoint -- matches what the",
                           "per-timepoint plots already do. With no timepoint-restricted",
                           "traits, aligned and compact produce identical output.")),

  make_option("--save_preview_png", type = "logical", default = TRUE),
  make_option("--preview_dpi", type = "double", default = 200)
)

opt <- parse_args(OptionParser(option_list = option_list))

# ------------------------------------------------------------
# Validate
# ------------------------------------------------------------
if (is.null(opt$project_root))        stop("--project_root is required")
if (is.null(opt$regions_file))        stop("--regions_file is required")
if (is.null(opt$dataset1_label))      stop("--dataset1_label is required")
if (is.null(opt$dataset1_stats_file)) stop("--dataset1_stats_file is required")
if (is.null(opt$exposure_outcome_trait_file)) stop("--exposure_outcome_trait_file is required")

if (!dir.exists(opt$project_root)) stop("project_root not found: ", opt$project_root)
if (!file.exists(opt$regions_file)) stop("regions_file not found: ", opt$regions_file)
if (!file.exists(opt$dataset1_stats_file)) stop("dataset1_stats_file not found: ", opt$dataset1_stats_file)
if (!file.exists(opt$exposure_outcome_trait_file)) stop("exposure_outcome_trait_file not found: ", opt$exposure_outcome_trait_file)

if (opt$p_thresh <= 0 || opt$p_thresh >= 1) stop("--p_thresh must be > 0 and < 1")
if (opt$min_timepoints_replicated < 1) stop("--min_timepoints_replicated must be >= 1")
selection_trait_category <- tolower(opt$selection_trait_category)
if (!selection_trait_category %in% c("exposure", "outcome")) {
  stop("--selection_trait_category must be one of: exposure, outcome")
}
panel_layout <- tolower(opt$panel_layout)
if (!panel_layout %in% c("horizontal", "vertical")) {
  stop("--panel_layout must be one of: horizontal, vertical")
}
panel_spacing_in <- opt$panel_spacing_in
if (is.na(panel_spacing_in)) {
  panel_spacing_in <- if (panel_layout == "vertical") 0.55 else 0.15
}
if (panel_spacing_in < 0) stop("--panel_spacing_in must be >= 0")
combined_column_mode <- tolower(opt$combined_column_mode)
if (!combined_column_mode %in% c("aligned", "compact")) {
  stop("--combined_column_mode must be one of: aligned, compact")
}
if (opt$min_exposure_traits < 1) stop("--min_exposure_traits must be >= 1")
if (opt$dot_size_min <= 0 || opt$dot_size_max <= opt$dot_size_min) {
  stop("--dot_size_max must be > --dot_size_min > 0")
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

dataset_labels <- sapply(dataset_inputs, `[[`, "label")

if (opt$min_timepoints_replicated > length(dataset_labels)) {
  stop("--min_timepoints_replicated (", opt$min_timepoints_replicated,
       ") cannot exceed the number of datasets provided (", length(dataset_labels), ")")
}

# ------------------------------------------------------------
# Output dirs
# ------------------------------------------------------------
{
  rfile_dir    <- dirname(opt$regions_file)
  region_label <- basename(rfile_dir)
  cpg_label    <- basename(dirname(rfile_dir))
  message("Derived cpg_label    : ", cpg_label)
  message("Derived region_label : ", region_label)
}

pipeline_root <- file.path(opt$project_root, "comethyl_output", "consensus")
step_dir      <- file.path(pipeline_root, "12c_me_trait_exposure_first_dotplots",
                            cpg_label, region_label, opt$adjustment_version)
cross_dir     <- file.path(step_dir, "cross_dataset")
dir.create(step_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(cross_dir, recursive = TRUE, showWarnings = FALSE)
message("Output root: ", step_dir)
message("  Shared/combined (all timepoints together): ", cross_dir)
for (lbl in dataset_labels) {
  dir.create(file.path(step_dir, lbl), recursive = TRUE, showWarnings = FALSE)
  message("  Per-timepoint: ", file.path(step_dir, lbl))
}

# ------------------------------------------------------------
# Load exposure/outcome trait definitions
# ------------------------------------------------------------
trait_def <- read_exposure_outcome_file(opt$exposure_outcome_trait_file)
exposure_traits <- unique(trait_def$trait[trait_def$category == "exposure"])
outcome_traits  <- unique(trait_def$trait[trait_def$category == "outcome"])

if (length(exposure_traits) == 0) stop("No exposure traits defined in exposure_outcome_trait_file.")
if (length(outcome_traits) == 0)  stop("No outcome traits defined in exposure_outcome_trait_file.")

message("Exposure traits: ", length(exposure_traits))
message("Outcome traits : ", length(outcome_traits))

# ------------------------------------------------------------
# Timepoint-restricted traits (optional third "timepoint" column in the
# trait file). trait_timepoint is a named vector: trait -> dataset label
# it's restricted to, or NA if unrestricted (applies to every panel).
# ------------------------------------------------------------
trait_timepoint <- setNames(trait_def$timepoint, trait_def$trait)

restricted_traits <- trait_def[!is.na(trait_def$timepoint), , drop = FALSE]
if (nrow(restricted_traits) > 0) {
  bad_timepoints <- setdiff(unique(restricted_traits$timepoint), dataset_labels)
  if (length(bad_timepoints) > 0) {
    warning(
      "exposure_outcome_trait_file has timepoint value(s) that don't match any ",
      "--dataset*_label provided to this run: ", paste(bad_timepoints, collapse = ", "),
      ". Traits restricted to a timepoint not in this run will never be plotted. ",
      "(dataset labels this run knows about: ", paste(dataset_labels, collapse = ", "), ")"
    )
  }
  message("Timepoint-restricted traits: ", nrow(restricted_traits), " (",
          paste(paste0(restricted_traits$trait, "->", restricted_traits$timepoint), collapse = ", "), ")")
}

# Returns the subset of `traits` valid for dataset ds_label -- i.e. traits
# with no timepoint restriction, or restricted specifically to ds_label.
# Used for per-timepoint plots, where an out-of-timepoint trait's column
# should be dropped entirely rather than just shown empty.
filter_traits_for_dataset <- function(traits, ds_label, trait_timepoint) {
  restrict_val <- trait_timepoint[traits]
  traits[is.na(restrict_val) | restrict_val == ds_label]
}

# ------------------------------------------------------------
# Load stats tables
# ------------------------------------------------------------
stats_list <- list()

for (ds in dataset_inputs) {
  df <- readr::read_tsv(ds$stats_file, show_col_types = FALSE)
  df <- validate_stats_table(df, paste0(ds$label, " stats_file"))

  cor_col <- get_cor_col(df)
  df$r_value <- suppressWarnings(as.numeric(df[[cor_col]]))
  df$p       <- suppressWarnings(as.numeric(df$p))
  df$module  <- as.character(df$module)
  df$trait   <- as.character(df$trait)

  if (isTRUE(opt$drop_grey_modules)) {
    n_before <- dplyr::n_distinct(df$module)
    df <- df[!is_grey_module(df$module), ]
    n_after <- dplyr::n_distinct(df$module)
    message("[", ds$label, "] Dropped grey module if present: ", n_before, " -> ", n_after, " module(s)")
  }

  stats_list[[ds$label]] <- df[, c("module", "trait", "p", "r_value")]

  message("Loaded stats for: ", ds$label, " (", nrow(df), " rows, ",
          length(unique(df$module)), " modules, ", length(unique(df$trait)), " traits)")
}

# ------------------------------------------------------------
# Compute robust module list (the core selection step).
# Category is configurable via --selection_trait_category:
# defaults to "exposure" but can be switched to "outcome".
# ------------------------------------------------------------
candidate_traits <- if (selection_trait_category == "exposure") exposure_traits else outcome_traits
message("\nSelection category: ", selection_trait_category,
        " (", length(candidate_traits), " candidate trait(s))")

pair_summary <- compute_pair_summary(
  stats_list             = stats_list,
  dataset_labels          = dataset_labels,
  candidate_traits        = candidate_traits,
  p_thresh                = opt$p_thresh,
  min_timepoints          = opt$min_timepoints_replicated,
  require_consistent_direction = isTRUE(opt$require_consistent_direction)
)

readr::write_tsv(pair_summary, file.path(step_dir, "exposure_outcome_pair_summary.tsv"))

robust_pairs <- pair_summary[pair_summary$robust %in% TRUE, ]

if (nrow(robust_pairs) == 0) {
  stop(
    "No module x ", selection_trait_category, "-trait pairs met the robustness criterion ",
    "(p < ", opt$p_thresh, " in >= ", opt$min_timepoints_replicated, " timepoint(s)",
    if (isTRUE(opt$require_consistent_direction)) ", consistent direction" else "",
    "). Try lowering --min_timepoints_replicated, raising --p_thresh, setting ",
    "--require_consistent_direction FALSE, or switching --selection_trait_category."
  )
}

module_summary <- robust_pairs %>%
  dplyr::group_by(module) %>%
  dplyr::summarise(
    n_robust_traits = dplyr::n(),
    robust_traits = paste(sort(unique(trait)), collapse = ";"),
    primary_trait = trait[which.min(match(trait, candidate_traits))],
    primary_trait_min_p = min_p[which.min(match(trait, candidate_traits))],
    .groups = "drop"
  ) %>%
  dplyr::filter(n_robust_traits >= opt$min_exposure_traits) %>%
  dplyr::mutate(primary_trait_rank = match(primary_trait, candidate_traits)) %>%
  dplyr::arrange(primary_trait_rank, module)

if (nrow(module_summary) == 0) {
  stop("No modules met --min_exposure_traits (", opt$min_exposure_traits, ") for ",
       "selection category '", selection_trait_category, "'. Lower this threshold and rerun.")
}

modules_ordered <- module_summary$module
n_modules <- length(modules_ordered)

message("\n", tools::toTitleCase(selection_trait_category), "-robust modules selected: ", n_modules)
message("(criterion: p < ", opt$p_thresh, " in >= ", opt$min_timepoints_replicated,
        " of ", length(dataset_labels), " timepoints",
        if (isTRUE(opt$require_consistent_direction)) ", consistent direction)" else ")")

write_vector_file(modules_ordered, file.path(step_dir, paste0(selection_trait_category, "_robust_modules.txt")))
readr::write_tsv(module_summary, file.path(step_dir, paste0(selection_trait_category, "_robust_module_summary.tsv")))


# ------------------------------------------------------------
# Dimension planning
# ------------------------------------------------------------
journal_width_in <- opt$full_width_mm / 25.4

message("\nPanel layout: ", panel_layout, " (panel_spacing = ", panel_spacing_in, " in)")
message("Journal full-width constraint: ", round(journal_width_in, 2),
        " in (", opt$full_width_mm, " mm, npj full-width column)")

resolve_dim <- function(user_val, auto_val) if (is.na(user_val)) auto_val else user_val

# ------------------------------------------------------------
# Generic renderer, used for both:
#   (a) the multi-panel "cross_dataset" version (all timepoints together)
#   (b) each single-panel per-timepoint version
# so the sizing/skip/save logic only needs to exist once.
# ------------------------------------------------------------
emit_dot_plot <- function(dataset_labels_subset, traits_ordered, group_blocks,
                           title, out_dir, file_stub, user_width, user_height,
                           label_for_msg, trait_timepoint_restrict = NULL,
                           column_mode = "aligned", per_panel_traits = NULL,
                           per_panel_group_blocks = NULL) {

  n_facets_here <- length(dataset_labels_subset)

  if (column_mode == "compact") {
    plot_df <- build_plot_df_per_panel(stats_list, dataset_labels_subset, modules_ordered,
                                        per_panel_traits, opt$p_thresh)
    n_traits_per_panel <- vapply(dataset_labels_subset, function(ds) length(per_panel_traits[[ds]]), numeric(1))
    grid_registration_df <- dplyr::bind_rows(lapply(dataset_labels_subset, function(ds) {
      traits_here <- per_panel_traits[[ds]]
      if (length(traits_here) == 0) return(NULL)
      data.frame(
        module = modules_ordered[1],
        # Explicit factor levels matching plot_df's per-dataset block order
        # (outcomes, then exposures) -- this layer's "trait" values feed
        # into the SAME shared x-scale as the main geom_point layer, so it
        # must carry the same intended order or ggplot's scale training
        # can revert to alphabetical for whichever values only appear here.
        trait = factor(traits_here, levels = traits_here),
        dataset = ds, stringsAsFactors = FALSE
      )
    }))
    if (!is.null(grid_registration_df) && nrow(grid_registration_df) > 0) {
      grid_registration_df$dataset <- factor(grid_registration_df$dataset, levels = dataset_labels_subset)
    }
  } else {
    plot_df <- build_plot_df(stats_list, dataset_labels_subset, modules_ordered, traits_ordered,
                              opt$p_thresh, trait_timepoint_restrict = trait_timepoint_restrict)
    n_traits_per_panel <- rep(length(traits_ordered), n_facets_here)
    grid_registration_df <- NULL
  }

  readr::write_tsv(plot_df, file.path(out_dir, paste0(file_stub, "_data.tsv")))

  if (nrow(plot_df) == 0) {
    message(label_for_msg, ": skipped -- no nominally significant associations (p < ",
            opt$p_thresh, ") for the ", n_modules, " selected module(s) in ",
            paste(dataset_labels_subset, collapse = ", "), ".")
    writeLines(
      paste0("Skipped: no significant associations for the ", n_modules,
             " module(s) selected via --selection_trait_category ", selection_trait_category,
             " in ", paste(dataset_labels_subset, collapse = ", ")),
      con = file.path(out_dir, paste0(file_stub, "_SKIPPED.txt"))
    )
    return(invisible(NULL))
  }

  max_trait_nchar <- if (column_mode == "compact") {
    all_here <- unique(unlist(per_panel_traits))
    if (length(all_here) == 0) 0 else max(nchar(all_here))
  } else {
    if (length(traits_ordered) == 0) 0 else max(nchar(traits_ordered))
  }

  h <- compute_recommended_height(
    n_modules        = n_modules,
    row_height_in    = opt$row_height_in,
    header_margin_in = opt$header_margin_in,
    max_height_in    = opt$max_height_in,
    n_facets         = n_facets_here,
    panel_layout     = panel_layout,
    panel_spacing_in = panel_spacing_in,
    max_trait_nchar  = max_trait_nchar,
    font_size        = opt$font_size
  )
  w <- compute_recommended_width(
    n_traits_per_panel = n_traits_per_panel,
    col_width_in        = opt$col_width_in,
    width_margin_in     = opt$width_margin_in,
    journal_width_in    = journal_width_in,
    panel_layout         = panel_layout
  )

  width_in  <- resolve_dim(user_width, w$width_in)
  height_in <- resolve_dim(user_height, h$height_in)

  message(label_for_msg, ": ", paste(n_traits_per_panel, collapse = "/"), " trait column(s) x ", n_facets_here,
          " panel(s) [", panel_layout, ", ", column_mode, "] -> ", round(width_in, 2), " x ", round(height_in, 2), " in",
          if (isTRUE(w$exceeds_journal_width) && panel_layout == "horizontal") paste0(
            " (EXCEEDS the ", round(journal_width_in, 2),
            " in journal full-width constraint -- consider --panel_layout vertical, ",
            "--drop_empty_columns TRUE, or a smaller trait set)"
          ) else if (isTRUE(h$exceeds_page)) paste0(
            " (exceeds --max_height_in ", opt$max_height_in, " in)"
          ) else "")

  p <- make_dot_plot(
    plot_df                = plot_df,
    modules_ordered         = modules_ordered,
    traits_ordered          = traits_ordered,
    dataset_labels           = dataset_labels_subset,
    group_blocks             = group_blocks,
    color_limit              = opt$color_limit,
    size_range                = size_range,
    font_family               = opt$font_family,
    font_size                  = opt$font_size,
    title                      = title,
    subtitle                   = paste0("(modules shown: ", selection_trait_category, "-robust hits [>= ",
                                opt$min_timepoints_replicated, " timepoints]; dot plotted only ",
                                "where nominally significant, p<", opt$p_thresh, ")"),
    add_panel_letters         = isTRUE(opt$add_panel_letters) && n_facets_here > 1,
    panel_layout               = panel_layout,
    panel_spacing_in           = panel_spacing_in,
    column_mode                = column_mode,
    grid_registration_df       = grid_registration_df,
    per_panel_group_blocks     = per_panel_group_blocks
  )

  save_plot(p, file.path(out_dir, file_stub), width_in, height_in, opt$font_family,
            isTRUE(opt$save_preview_png), opt$preview_dpi)
  message("Saved: ", file_stub, ".pdf (", round(width_in, 2), " x ", round(height_in, 2),
          " in) -> ", out_dir)

  invisible(list(width_in = width_in, height_in = height_in))
}

size_range <- c(opt$dot_size_min, opt$dot_size_max)

# ------------------------------------------------------------
# Resolve trait columns per category (optionally dropping empty ones)
# ------------------------------------------------------------
exposure_traits_plot <- exposure_traits
if (isTRUE(opt$drop_empty_columns)) {
  tmp <- build_plot_df(stats_list, dataset_labels, modules_ordered, exposure_traits, opt$p_thresh)
  exposure_traits_plot <- exposure_traits[exposure_traits %in% unique(as.character(tmp$trait))]
  if (length(exposure_traits_plot) == 0) exposure_traits_plot <- exposure_traits
}

outcome_traits_plot <- outcome_traits
if (isTRUE(opt$drop_empty_columns)) {
  tmp <- build_plot_df(stats_list, dataset_labels, modules_ordered, outcome_traits, opt$p_thresh)
  outcome_traits_plot <- outcome_traits[outcome_traits %in% unique(as.character(tmp$trait))]
  if (length(outcome_traits_plot) == 0) outcome_traits_plot <- outcome_traits
}

combined_traits <- c(outcome_traits_plot, exposure_traits_plot)

plot_specs <- list(
  Exposures = list(
    traits       = exposure_traits_plot,
    group_blocks = list(list(label = opt$exposure_group_label, start = 1, end = length(exposure_traits_plot))),
    title        = "Consensus methylation module associations with exposures",
    user_width   = opt$exposure_plot_width,
    user_height  = opt$exposure_plot_height
  ),
  Outcomes = list(
    traits       = outcome_traits_plot,
    group_blocks = list(list(label = opt$outcome_group_label, start = 1, end = length(outcome_traits_plot))),
    title        = "Consensus methylation module associations with pregnancy/PPWR outcomes",
    user_width   = opt$outcome_plot_width,
    user_height  = opt$outcome_plot_height
  ),
  Combined = list(
    traits       = combined_traits,
    group_blocks = list(
      list(label = opt$outcome_group_label,  start = 1, end = length(outcome_traits_plot)),
      list(label = opt$exposure_group_label, start = length(outcome_traits_plot) + 1, end = length(combined_traits))
    ),
    title        = "Consensus methylation module associations with outcomes and exposures",
    user_width   = opt$combined_plot_width,
    user_height  = opt$combined_plot_height
  )
)

# ------------------------------------------------------------
# Render: shared/combined (cross_dataset/) version of each plot --
# all timepoints together, layout controlled by --panel_layout.
# Per-plot user width/height overrides (--exposure_plot_width etc.) apply
# ONLY to this combined version; per-timepoint versions below are always
# auto-sized single panels.
# ------------------------------------------------------------
# ------------------------------------------------------------
# For a given dataset, returns the per-category trait columns and
# group_blocks applicable to it (dropping, not blanking, out-of-timepoint
# columns). Shared by per-timepoint plots and by the combined plot's
# "compact" column mode -- both need exactly the same per-panel filtering.
# ------------------------------------------------------------
build_ds_plot_specs <- function(ds_label) {
  exposure_here <- filter_traits_for_dataset(exposure_traits_plot, ds_label, trait_timepoint)
  outcome_here  <- filter_traits_for_dataset(outcome_traits_plot, ds_label, trait_timepoint)
  combined_here <- c(outcome_here, exposure_here)

  list(
    Exposures = list(
      traits       = exposure_here,
      group_blocks = list(list(label = opt$exposure_group_label, start = 1, end = length(exposure_here))),
      title        = plot_specs$Exposures$title
    ),
    Outcomes = list(
      traits       = outcome_here,
      group_blocks = list(list(label = opt$outcome_group_label, start = 1, end = length(outcome_here))),
      title        = plot_specs$Outcomes$title
    ),
    Combined = list(
      traits       = combined_here,
      group_blocks = list(
        list(label = opt$outcome_group_label,  start = 1, end = length(outcome_here)),
        list(label = opt$exposure_group_label, start = length(outcome_here) + 1, end = length(combined_here))
      ),
      title        = plot_specs$Combined$title
    )
  )
}

message("\n--- Shared/combined plots (cross_dataset/) ---")
message("  combined_column_mode = ", combined_column_mode)
if (combined_column_mode == "compact") {
  per_ds_specs <- setNames(lapply(dataset_labels, build_ds_plot_specs), dataset_labels)

  for (plot_name in names(plot_specs)) {
    per_panel_traits       <- setNames(lapply(dataset_labels, function(ds) per_ds_specs[[ds]][[plot_name]]$traits), dataset_labels)
    per_panel_group_blocks <- setNames(lapply(dataset_labels, function(ds) per_ds_specs[[ds]][[plot_name]]$group_blocks), dataset_labels)

    emit_dot_plot(
      dataset_labels_subset   = dataset_labels,
      traits_ordered          = plot_specs[[plot_name]]$traits, # unused in compact mode, kept for messaging fallback
      group_blocks            = NULL,
      title                   = plot_specs[[plot_name]]$title,
      out_dir                 = cross_dir,
      file_stub               = paste0("Dotplot_", plot_name),
      user_width              = plot_specs[[plot_name]]$user_width,
      user_height             = plot_specs[[plot_name]]$user_height,
      label_for_msg           = paste0(plot_name, " (combined)"),
      column_mode             = "compact",
      per_panel_traits        = per_panel_traits,
      per_panel_group_blocks  = per_panel_group_blocks
    )
  }
} else {
  for (plot_name in names(plot_specs)) {
    spec <- plot_specs[[plot_name]]
    emit_dot_plot(
      dataset_labels_subset     = dataset_labels,
      traits_ordered            = spec$traits,
      group_blocks              = spec$group_blocks,
      title                     = spec$title,
      out_dir                   = cross_dir,
      file_stub                 = paste0("Dotplot_", plot_name),
      user_width                = spec$user_width,
      user_height                = spec$user_height,
      label_for_msg             = paste0(plot_name, " (combined)"),
      trait_timepoint_restrict  = trait_timepoint,
      column_mode                = "aligned"
    )
  }
}

# ------------------------------------------------------------
# Render: per-timepoint (single-panel) version of each plot, one folder
# per dataset label -- matches the per-dataset folder convention used
# elsewhere in this project (e.g. script 12b).
# ------------------------------------------------------------
if (isTRUE(opt$per_timepoint_plots)) {
  message("\n--- Per-timepoint plots ---")
  for (ds_label in dataset_labels) {
    ds_dir <- file.path(step_dir, ds_label)

    ds_plot_specs <- build_ds_plot_specs(ds_label)

    n_dropped <- length(exposure_traits_plot) - length(ds_plot_specs$Exposures$traits)
    if (n_dropped > 0) {
      message("[", ds_label, "] dropped ", n_dropped, " exposure column(s) restricted to a different timepoint")
    }

    for (plot_name in names(ds_plot_specs)) {
      spec <- ds_plot_specs[[plot_name]]
      emit_dot_plot(
        dataset_labels_subset    = ds_label,
        traits_ordered           = spec$traits,
        group_blocks             = spec$group_blocks,
        title                    = spec$title,
        out_dir                  = ds_dir,
        file_stub                = paste0("Dotplot_", plot_name, "_", ds_label),
        user_width               = NA,
        user_height              = NA,
        label_for_msg            = paste0(plot_name, " (", ds_label, ")"),
        trait_timepoint_restrict = NULL, # already column-filtered above; no row-level restriction needed
        column_mode              = "aligned" # single panel -- aligned/compact are equivalent here
      )
    }
  }
} else {
  message("\n--per_timepoint_plots FALSE -- skipping per-timepoint folders.")
}

# ------------------------------------------------------------
# Logs
# ------------------------------------------------------------
write_log_lines(
  c(
    paste("project_root:", opt$project_root),
    paste("adjustment_version:", opt$adjustment_version),
    paste("datasets:", paste(dataset_labels, collapse = ", ")),
    paste("exposure_outcome_trait_file:", opt$exposure_outcome_trait_file),
    paste("n_exposure_traits_defined:", length(exposure_traits)),
    paste("n_outcome_traits_defined:", length(outcome_traits)),
    paste("drop_grey_modules:", opt$drop_grey_modules),
    paste("p_thresh:", opt$p_thresh),
    paste("min_timepoints_replicated:", opt$min_timepoints_replicated),
    paste("require_consistent_direction:", opt$require_consistent_direction),
    paste("selection_trait_category:", selection_trait_category),
    paste("min_exposure_traits:", opt$min_exposure_traits),
    paste("drop_empty_columns:", opt$drop_empty_columns),
    paste("n_modules_selected:", n_modules),
    paste("panel_layout:", panel_layout),
    paste("panel_spacing_in:", panel_spacing_in),
    paste("combined_column_mode:", combined_column_mode),
    paste("per_timepoint_plots:", opt$per_timepoint_plots),
    paste("color_limit:", ifelse(is.na(opt$color_limit), "auto", opt$color_limit)),
    paste("font_family:", opt$font_family),
    paste("font_size:", opt$font_size),
    paste("full_width_mm:", opt$full_width_mm),
    paste("row_height_in:", opt$row_height_in),
    paste("date:", as.character(Sys.time()))
  ),
  file.path(step_dir, "run_parameters.txt")
)

writeLines(
  capture.output(sessionInfo()),
  con = file.path(step_dir, "sessionInfo.txt")
)

message("\n\u2713 Script 12C complete: exposure-first presentation dot plots finished")
message("Outputs saved under: ", step_dir)
message("  Shared/combined: ", cross_dir)
if (isTRUE(opt$per_timepoint_plots)) {
  message("  Per-timepoint: ", paste(file.path(step_dir, dataset_labels), collapse = ", "))
}

#This worked but too many timepoint labels
# #!/usr/bin/env Rscript

# # ================================================================
# # SCRIPT 12C: Consensus ME-Trait Category-First Presentation Dot Plots
# #
# # Pipeline: comethyl WGBS consensus analysis
# #
# # PURPOSE
# #   Builds a single "<category>-robust" consensus module list (outcome
# #   by default; can be switched to exposure via --selection_trait_category),
# #   then renders three publication-style dot plots (bicor = fill color,
# #   -log10(p) = dot size, dot drawn only where nominally significant)
# #   sharing that SAME module list and row order:
# #     1) EXPOSURES only
# #     2) OUTCOMES only
# #     3) COMBINED (outcomes block + exposures block, divided by a
# #        vertical line, matching the presentation figure style)
# #
# #   Module selection logic:
# #     A module is kept if it has >= --min_exposure_traits trait(s) from
# #     the --selection_trait_category category that are nominally
# #     significant (p < --p_thresh) in >= --min_timepoints_replicated of
# #     the datasets provided, with (by default) a consistent bicor sign
# #     across those significant timepoints. This mirrors the robustness
# #     convention already used elsewhere in this project. The resulting
# #     module list and row order are then reused, unchanged, for the
# #     other two plots -- so a module can appear with a blank row in
# #     whichever plot corresponds to the category it was NOT selected on.
# #
# # REQUIRED INPUTS
# #   --project_root               : root directory of the analysis project
# #   --regions_file                : path to Filtered_Regions.txt from script 02
# #   --dataset1_label              : label for dataset 1
# #   --dataset1_stats_file         : ME-trait stats TSV from script 12a
# #   --exposure_outcome_trait_file : TSV/CSV/TXT with columns category, trait
# #                                   (category must be "exposure" or "outcome"),
# #                                   plus an OPTIONAL third column "timepoint".
# #                                   Base 2-column format is the same file used
# #                                   by 12b's --module_order_trait_file.
# #                                   timepoint restricts a trait to a single
# #                                   dataset's panel -- its value must exactly
# #                                   match one of --dataset1_label/2/3_label
# #                                   (e.g. "Baseline", "36-38wk", "Postpartum").
# #                                   Leave blank/omit the column entirely for a
# #                                   trait that should appear in every panel
# #                                   (this is the old, still-default, behavior).
# #                                   In the cross_dataset/ combined figure, a
# #                                   restricted trait's column stays in every
# #                                   panel (for alignment) but only shows a dot
# #                                   in its own timepoint's panel. In the
# #                                   per-timepoint plots, a restricted trait's
# #                                   column is dropped entirely from any panel
# #                                   it doesn't belong to.
# #   --selection_trait_category    : exposure or outcome [default = outcome].
# #                                   Controls which trait category drives module
# #                                   selection AND row order. The other category
# #                                   is still plotted in its own dot plot and in
# #                                   the combined plot -- switching this only
# #                                   changes which category decides which modules
# #                                   are kept and how rows are ranked.
# #
# # OPTIONAL INPUTS
# #   --dataset2_label / --dataset2_stats_file : dataset 2
# #   --dataset3_label / --dataset3_stats_file : dataset 3
# #   --adjustment_version         : label matching script 12a run [default = unadjusted]
# #
# #   --drop_grey_modules          : TRUE/FALSE [default = TRUE]
# #   --p_thresh                   : nominal significance threshold [default = 0.05]
# #   --min_timepoints_replicated  : # of datasets the SAME exposure trait must be
# #                                   nominally significant in for a module to count
# #                                   as exposure-robust [default = 2]
# #   --require_consistent_direction : TRUE/FALSE; require the same bicor sign across
# #                                   the replicated timepoints [default = TRUE]
# #   --min_exposure_traits        : minimum number of robust exposure traits a module
# #                                   must have to be kept [default = 1]
# #   --drop_empty_columns         : TRUE/FALSE; drop trait columns with zero
# #                                   significant dots anywhere in that plot [default = FALSE]
# #
# #   --color_limit                 : symmetric bicor color-scale limit. If not set,
# #                                    auto-computed from the max |bicor| actually plotted.
# #   --dot_size_min / --dot_size_max : point size range in mm [default = 1.5 / 6]
# #
# #   --font_family                 : font family passed to the plotting device and
# #                                    theme [default = Helvetica]. If Helvetica is not
# #                                    resolvable by your Cairo/font setup, rerun with
# #                                    "sans" -- this script does not auto-detect fonts.
# #   --font_size                   : base point size [default = 8, per npj Women's Health spec]
# #
# #   --full_width_mm               : figure width in mm [default = 170, npj full-width]
# #   --row_height_in               : inches of height per module row, used to
# #                                    auto-compute figure height [default = 0.155]
# #   --header_margin_in            : inches reserved above the plot area for
# #                                    title/legend/group labels [default = 1.0]
# #   --max_height_in               : height (in) above which a warning is printed
# #                                    suggesting the figure won't fit a single
# #                                    journal page [default = 9.5]
# #   --exposure_plot_width/height/--outcome_plot_width/height/
# #   --combined_plot_width/height  : inches; override the auto-computed size for
# #                                    that specific plot. [default = NULL = auto]
# #
# #   --exposure_group_label        : label text for the exposures column block
# #                                    [default = "Exposures"]
# #   --outcome_group_label         : label text for the outcomes column block
# #                                    [default = "Outcomes"]
# #   --add_panel_letters           : TRUE/FALSE; prefix facet strips with a), b), c)
# #                                    [default = TRUE]
# #
# #   --save_preview_png            : TRUE/FALSE; also save a PNG for slides
# #                                    [default = TRUE]
# #   --preview_dpi                 : DPI for the PNG preview [default = 200]
# #
# # OUTPUTS
# #   comethyl_output/consensus/12c_me_trait_exposure_first_dotplots/
# #       <cpg_label>/<region_label>/<adjustment_version>/
# #           <selection_trait_category>_robust_modules.txt
# #           <selection_trait_category>_robust_module_summary.tsv
# #           exposure_outcome_pair_summary.tsv   (all module x trait pairs, robust or not)
# #           run_parameters.txt
# #           sessionInfo.txt
# #
# #           cross_dataset/   <- SHARED/COMBINED: all timepoints in one figure
# #               Dotplot_Exposures.pdf  (+ .png if --save_preview_png)
# #               Dotplot_Outcomes.pdf   (+ .png)
# #               Dotplot_Combined.pdf   (+ .png)
# #               Dotplot_Exposures_data.tsv   (rows actually plotted, i.e. nominally significant)
# #               Dotplot_Outcomes_data.tsv
# #               Dotplot_Combined_data.tsv
# #
# #           <dataset1_label>/   <- PER-TIMEPOINT: single-panel versions, one folder per dataset
# #               Dotplot_Exposures_<dataset1_label>.pdf  (+ .png, + _data.tsv)
# #               Dotplot_Outcomes_<dataset1_label>.pdf
# #               Dotplot_Combined_<dataset1_label>.pdf
# #           <dataset2_label>/  (same three files, if dataset2 provided)
# #           <dataset3_label>/  (same three files, if dataset3 provided)
# #
# #   Skipped plots (no nominally significant rows for that module/trait/
# #   dataset combination) get a "..._SKIPPED.txt" note instead of a PDF.
# #
# # NOTES
# #   - cpg_label and region_label are derived from --regions_file path.
# #   - Vector PDFs are written with grDevices::cairo_pdf, matching the
# #     project's Illustrator-assembly convention (panels placed as
# #     separate vector PDFs). Convert to 600 DPI TIFF in Illustrator/
# #     Photoshop at final placement time, per npj Women's Health specs.
# #   - Figure dimensions are auto-computed from the module count, trait
# #     column count, and --panel_layout. This script prints the computed
# #     size and warns when it won't fit journal full-width/a single page --
# #     check that message before assuming a figure is presentation- or
# #     print-ready.
# #   - --panel_layout vertical stacks timepoint panels top-to-bottom instead
# #     of side-by-side, trading width for height -- useful when a wide
# #     trait-column count won't fit n_timepoints panels at journal
# #     full-width. Per-timepoint plots (single dataset each) are unaffected
# #     by --panel_layout since they only ever have one panel.
# #   - --exposure_plot_width/height, --outcome_plot_width/height, and
# #     --combined_plot_width/height override the auto-computed size for the
# #     cross_dataset/ (combined) version ONLY. Per-timepoint plots are
# #     always auto-sized.
# #
# # EXAMPLE
# #   Rscript $SCRIPTS/consensus/12c_me_trait_exposure_first_dotplots_consensus.R \
# #     --project_root $PROJECT \
# #     --regions_file $OUT/02_reference_region_filter/Baseline/$CPG_LABEL/$REGION_LABEL/Filtered_Regions.txt \
# #     --adjustment_version v2_exclude_protected_pcs \
# #     --dataset1_label Baseline \
# #     --dataset1_stats_file $OUT/12a_me_trait_analysis/Baseline/$CPG_LABEL/$REGION_LABEL/v2_exclude_protected_pcs/ME_Trait_Correlation_Stats_Bicor.tsv \
# #     --dataset2_label 36-38wk \
# #     --dataset2_stats_file $OUT/12a_me_trait_analysis/36-38wk/$CPG_LABEL/$REGION_LABEL/v2_exclude_protected_pcs/ME_Trait_Correlation_Stats_Bicor.tsv \
# #     --dataset3_label Postpartum \
# #     --dataset3_stats_file $OUT/12a_me_trait_analysis/Postpartum/$CPG_LABEL/$REGION_LABEL/v2_exclude_protected_pcs/ME_Trait_Correlation_Stats_Bicor.tsv \
# #     --exposure_outcome_trait_file $CONFIG/exposure_outcome_traits.tsv
# # ================================================================

# message("Starting Script 12c \u2713")

# suppressPackageStartupMessages({
#   library(optparse)
#   library(readr)
#   library(dplyr)
#   library(ggplot2)
#   library(scales)
# })

# # ------------------------------------------------------------
# # Load helper.R
# # ------------------------------------------------------------
# script_file_arg <- commandArgs(trailingOnly = FALSE)[
#   grep("^--file=", commandArgs(trailingOnly = FALSE))
# ]

# if (length(script_file_arg) == 0) {
#   stop("Could not determine script path.")
# }

# script_dir  <- dirname(normalizePath(sub("^--file=", "", script_file_arg[1])))
# helper_file <- file.path(script_dir, "helper.R")

# if (!file.exists(helper_file)) {
#   stop("helper.R not found: ", helper_file)
# }

# source(helper_file)

# required_helpers <- c("write_vector_file", "write_log_lines")
# missing_helpers <- required_helpers[!vapply(required_helpers, exists, logical(1), mode = "function")]
# if (length(missing_helpers) > 0) {
#   stop("Missing helper functions in helper.R: ", paste(missing_helpers, collapse = ", "))
# }

# # ggplot2's internal mm<->pt conversion constant, reimplemented locally
# # so this script does not depend on an unexported ggplot2 internal.
# PT_TO_MM <- 25.4 / 72.27

# # ------------------------------------------------------------
# # Helpers
# # ------------------------------------------------------------
# get_cor_col <- function(stats_df) {
#   if ("bicor" %in% colnames(stats_df)) return("bicor")
#   if ("cor" %in% colnames(stats_df)) return("cor")
#   if ("correlation" %in% colnames(stats_df)) return("correlation")
#   stop("Stats table must contain one of: bicor, cor, correlation")
# }

# standardize_p_col <- function(stats_df) {
#   if ("p" %in% colnames(stats_df)) return(stats_df)
#   candidates <- c("pvalue", "p.value", "p_value", "P.Value", "P.value")
#   hit <- candidates[candidates %in% colnames(stats_df)][1]
#   if (is.na(hit)) return(stats_df)
#   colnames(stats_df)[match(hit, colnames(stats_df))] <- "p"
#   stats_df
# }

# validate_stats_table <- function(df, label) {
#   df <- standardize_p_col(df)
#   required <- c("module", "trait", "p")
#   missing  <- setdiff(required, colnames(df))
#   if (length(missing) > 0) {
#     stop(label, " missing columns: ", paste(missing, collapse = ", "))
#   }
#   get_cor_col(df)
#   df
# }

# is_grey_module <- function(x) {
#   tolower(sub("^ME", "", as.character(x))) == "grey"
# }

# read_exposure_outcome_file <- function(file) {
#   if (is.null(file) || !nzchar(file)) {
#     stop("--exposure_outcome_trait_file is required")
#   }
#   if (!file.exists(file)) stop("exposure_outcome_trait_file not found: ", file)

#   ext <- tolower(tools::file_ext(file))
#   x <- switch(ext,
#     csv = read.csv(file, stringsAsFactors = FALSE, check.names = FALSE),
#     tsv = read.delim(file, stringsAsFactors = FALSE, check.names = FALSE),
#     txt = read.delim(file, stringsAsFactors = FALSE, check.names = FALSE),
#     read.delim(file, stringsAsFactors = FALSE, check.names = FALSE)
#   )

#   required <- c("category", "trait")
#   missing <- setdiff(required, colnames(x))
#   if (length(missing) > 0) {
#     stop("exposure_outcome_trait_file missing required columns: ",
#          paste(missing, collapse = ", "))
#   }

#   x$category <- tolower(trimws(as.character(x$category)))
#   x$trait <- trimws(as.character(x$trait))

#   # Optional third column: restricts a trait to a single dataset/timepoint's
#   # panel (value must exactly match one of --dataset1_label/2/3_label). Blank
#   # or absent = trait applies to every panel, which preserves old 2-column
#   # file behavior unchanged.
#   if ("timepoint" %in% colnames(x)) {
#     x$timepoint <- trimws(as.character(x$timepoint))
#     x$timepoint[x$timepoint == ""] <- NA_character_
#   } else {
#     x$timepoint <- NA_character_
#   }

#   x <- x[x$category %in% c("exposure", "outcome") &
#          !is.na(x$trait) & x$trait != "", , drop = FALSE]

#   if (nrow(x) == 0) {
#     stop("exposure_outcome_trait_file has no valid exposure/outcome rows.")
#   }

#   x <- x[!duplicated(x[, c("category", "trait")]), ]
#   x[, c("category", "trait", "timepoint")]
# }

# # ------------------------------------------------------------
# # Exposure-first robust module/trait pair computation
# # ------------------------------------------------------------
# compute_pair_summary <- function(stats_list, dataset_labels, candidate_traits,
#                                   p_thresh, min_timepoints, require_consistent_direction) {

#   long <- dplyr::bind_rows(lapply(dataset_labels, function(ds) {
#     df <- stats_list[[ds]]
#     df <- df[df$trait %in% candidate_traits, c("module", "trait", "p", "r_value")]
#     if (nrow(df) == 0) return(df)
#     df$dataset <- ds
#     df
#   }))

#   if (nrow(long) == 0) {
#     return(data.frame(
#       module = character(0), trait = character(0), n_sig = integer(0),
#       consistent_direction = logical(0), timepoints_sig = character(0),
#       min_p = numeric(0), max_abs_r = numeric(0), robust = logical(0),
#       stringsAsFactors = FALSE
#     ))
#   }

#   long$sig <- !is.na(long$p) & long$p < p_thresh

#   pair_summary <- long %>%
#     dplyr::group_by(module, trait) %>%
#     dplyr::summarise(
#       n_sig = sum(sig, na.rm = TRUE),
#       n_unique_signs = dplyr::n_distinct(sign(r_value[sig])),
#       timepoints_sig = paste(sort(dataset[sig]), collapse = ";"),
#       min_p = suppressWarnings(min(p[sig], na.rm = TRUE)),
#       max_abs_r = suppressWarnings(max(abs(r_value[sig]), na.rm = TRUE)),
#       .groups = "drop"
#     ) %>%
#     dplyr::mutate(
#       consistent_direction = n_unique_signs <= 1,
#       min_p = ifelse(is.infinite(min_p), NA_real_, min_p),
#       max_abs_r = ifelse(is.infinite(max_abs_r) | is.nan(max_abs_r), NA_real_, max_abs_r),
#       robust = n_sig >= min_timepoints & (!require_consistent_direction | consistent_direction)
#     )

#   as.data.frame(pair_summary)
# }

# # ------------------------------------------------------------
# # Build the long data frame actually plotted (significant rows only,
# # with full factor levels retained for consistent column/row grids)
# # ------------------------------------------------------------
# # ------------------------------------------------------------
# # Build the long data frame actually plotted (significant rows only,
# # with full factor levels retained for consistent column/row grids).
# # trait_timepoint_restrict, if given, is a named character vector
# # (trait -> dataset label) that blanks out rows for a trait in any
# # dataset other than its assigned one, WITHOUT removing that trait's
# # column from the grid -- used for the multi-panel combined view so
# # columns stay aligned across panels while cross-timepoint spillover
# # (e.g. a wk36-only exposure showing a dot in the Baseline panel) is
# # suppressed. Traits absent from the vector, or with an NA entry, are
# # unrestricted and behave exactly as before.
# # ------------------------------------------------------------
# build_plot_df <- function(stats_list, dataset_labels, modules_ordered, traits_ordered,
#                            p_thresh, trait_timepoint_restrict = NULL) {
#   long <- dplyr::bind_rows(lapply(dataset_labels, function(ds) {
#     df <- stats_list[[ds]]
#     df <- df[df$module %in% modules_ordered & df$trait %in% traits_ordered, ]
#     if (nrow(df) == 0) return(df)
#     df$dataset <- ds
#     df
#   }))

#   if (!is.null(trait_timepoint_restrict) && nrow(long) > 0) {
#     restrict_val <- trait_timepoint_restrict[long$trait]
#     keep <- is.na(restrict_val) | (restrict_val == long$dataset)
#     keep[is.na(keep)] <- TRUE # trait not present in the restriction map at all
#     long <- long[keep, ]
#   }

#   long <- long[!is.na(long$p) & long$p < p_thresh, ]
#   long$neglogp <- -log10(pmax(long$p, 1e-300))
#   long$module <- factor(long$module, levels = modules_ordered)
#   long$trait <- factor(long$trait, levels = traits_ordered)
#   long$dataset <- factor(long$dataset, levels = dataset_labels)
#   long
# }

# # ------------------------------------------------------------
# # Dimension recommendation
# # ------------------------------------------------------------
# GROUP_LABEL_EXTRA_ROWS <- 3.2 # must match the y-scale expansion in make_dot_plot()
#                                # (add = c(0.6, 2.6)); budgeted here as real
#                                # height so it doesn't compress the module rows.

# # Empirically-determined floor: ggplot's fixed-size chrome (title, subtitle,
# # top margin, facet strip) costs roughly the same ABSOLUTE amount of space
# # regardless of how few module rows there are. With very few modules
# # (n_modules ~1-3), that fixed cost eats a large enough fraction of a small
# # total figure height that the remaining panel space compresses module row
# # labels into unreadable overlap -- confirmed by direct visual testing at
# # n_modules=2, where anything below ~2.1in total height overlapped and 2.3in
# # was clean. This floor is intentionally a bit above that tested threshold.
# MIN_FIGURE_HEIGHT_IN <- 2.3

# # horizontal: panels sit side by side -> height is ONE panel's worth of rows.
# # vertical:   panels stack top to bottom -> height is N panels' worth of rows,
# #             plus spacing between each stacked panel for its own facet strip
# #             + group labels.
# compute_recommended_height <- function(n_modules, row_height_in, header_margin_in,
#                                         max_height_in, n_facets, panel_layout,
#                                         panel_spacing_in) {
#   panel_multiplier <- if (panel_layout == "vertical") max(n_facets, 1) else 1
#   spacing_total_in <- if (panel_layout == "vertical" && n_facets > 1) (n_facets - 1) * panel_spacing_in else 0

#   height_in <- header_margin_in +
#     (panel_multiplier * (n_modules + GROUP_LABEL_EXTRA_ROWS) * row_height_in) +
#     spacing_total_in

#   height_in <- max(height_in, MIN_FIGURE_HEIGHT_IN)

#   list(
#     height_in = height_in,
#     exceeds_page = height_in > max_height_in
#   )
# }

# # horizontal: panels sit side by side -> width is N panels' worth of columns.
# # vertical:   panels stack top to bottom -> width is ONE panel's worth of
# #             columns (this is the whole point of switching to vertical --
# #             trading width for height when there are too many trait columns
# #             to fit n_timepoints panels side by side at journal full-width).
# compute_recommended_width <- function(n_traits, n_facets, col_width_in,
#                                        width_margin_in, journal_width_in, panel_layout) {
#   panel_multiplier <- if (panel_layout == "vertical") 1 else max(n_facets, 1)
#   comfortable_width_in <- width_margin_in + (panel_multiplier * n_traits * col_width_in)
#   list(
#     comfortable_width_in = comfortable_width_in,
#     journal_width_in     = journal_width_in,
#     width_in             = max(comfortable_width_in, journal_width_in),
#     exceeds_journal_width = comfortable_width_in > journal_width_in
#   )
# }

# # ------------------------------------------------------------
# # Dot plot builder
# # ------------------------------------------------------------
# make_dot_plot <- function(plot_df, modules_ordered, traits_ordered, dataset_labels,
#                            group_blocks = NULL, color_limit = NA,
#                            size_range = c(1.5, 6), font_family = "Helvetica",
#                            font_size = 8, title = NULL, subtitle = NULL,
#                            add_panel_letters = TRUE, panel_layout = "horizontal",
#                            panel_spacing_in = 0.15) {

#   n_mod <- length(modules_ordered)

#   if (is.na(color_limit) || is.null(color_limit)) {
#     color_limit <- suppressWarnings(max(abs(plot_df$r_value), na.rm = TRUE))
#     if (!is.finite(color_limit) || color_limit <= 0) color_limit <- 0.6
#   }

#   facet_labels <- dataset_labels
#   if (isTRUE(add_panel_letters)) {
#     facet_labels <- paste0(letters[seq_along(dataset_labels)], ") ", dataset_labels)
#   }
#   names(facet_labels) <- dataset_labels

#   facet_layer <- if (panel_layout == "vertical") {
#     facet_wrap(~dataset, ncol = 1, labeller = as_labeller(facet_labels))
#   } else {
#     facet_wrap(~dataset, nrow = 1, labeller = as_labeller(facet_labels))
#   }

#   p <- ggplot(plot_df, aes(x = trait, y = module)) +
#     geom_point(aes(fill = r_value, size = neglogp), shape = 21, color = "black", stroke = 0.3) +
#     scale_fill_gradient2(
#       low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0,
#       limits = c(-color_limit, color_limit), oob = scales::squish,
#       name = "bicor (r)"
#     ) +
#     scale_size_continuous(
#       range = size_range,
#       name = expression(-log[10](p)),
#       guide = guide_legend(override.aes = list(fill = "grey40"))
#     ) +
#     scale_x_discrete(limits = traits_ordered, drop = FALSE) +
#     scale_y_discrete(
#       limits = rev(modules_ordered), drop = FALSE,
#       # Reserve fixed headroom (in row-units, independent of n_modules) above
#       # the top row so the group-block labels never collide with the facet
#       # strip text, regardless of how many modules are plotted.
#       expand = expansion(add = c(0.6, 2.6))
#     ) +
#     facet_layer +
#     labs(x = NULL, y = NULL, title = title, subtitle = subtitle) +
#     theme_minimal(base_family = font_family, base_size = font_size) +
#     theme(
#       axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = font_size - 1),
#       axis.text.y = element_text(size = font_size - 1),
#       strip.text = element_text(face = "bold", size = font_size + 1),
#       panel.grid.major = element_line(color = "grey85", linewidth = 0.2),
#       panel.grid.minor = element_blank(),
#       plot.title = element_text(size = font_size + 2, hjust = 0),
#       plot.subtitle = element_text(size = font_size - 1),
#       legend.title = element_text(size = font_size),
#       legend.text = element_text(size = font_size - 1),
#       plot.margin = margin(t = 30, r = 5, b = 5, l = 5, unit = "pt"),
#       panel.spacing = grid::unit(panel_spacing_in, "in")
#     ) +
#     coord_cartesian(clip = "off")

#   if (!is.null(group_blocks) && length(group_blocks) > 0) {
#     for (i in seq_along(group_blocks)) {
#       blk <- group_blocks[[i]]
#       p <- p + annotate(
#         "text",
#         x = mean(c(blk$start, blk$end)),
#         y = n_mod + 1.4,
#         label = blk$label,
#         fontface = "bold",
#         size = (font_size - 0.5) * PT_TO_MM,
#         hjust = 0.5
#       )
#       if (i < length(group_blocks)) {
#         p <- p + geom_vline(xintercept = blk$end + 0.5, color = "black", linewidth = 0.6)
#       }
#     }
#   }

#   p
# }

# save_plot <- function(p, out_prefix, width_in, height_in, font_family,
#                        save_png, preview_dpi) {
#   pdf_file <- paste0(out_prefix, ".pdf")
#   grDevices::cairo_pdf(pdf_file, width = width_in, height = height_in, family = font_family)
#   print(p)
#   grDevices::dev.off()

#   if (isTRUE(save_png)) {
#     png_file <- paste0(out_prefix, ".png")
#     ggplot2::ggsave(
#       png_file, plot = p, width = width_in, height = height_in,
#       dpi = preview_dpi, units = "in", bg = "white"
#     )
#   }
#   invisible(pdf_file)
# }

# # ------------------------------------------------------------
# # Parse args
# # ------------------------------------------------------------
# option_list <- list(
#   make_option("--project_root",           type = "character"),
#   make_option("--regions_file",           type = "character"),
#   make_option("--adjustment_version",     type = "character", default = "unadjusted"),

#   make_option("--dataset1_label",         type = "character"),
#   make_option("--dataset1_stats_file",    type = "character"),
#   make_option("--dataset2_label",         type = "character", default = NULL),
#   make_option("--dataset2_stats_file",    type = "character", default = NULL),
#   make_option("--dataset3_label",         type = "character", default = NULL),
#   make_option("--dataset3_stats_file",    type = "character", default = NULL),

#   make_option("--exposure_outcome_trait_file", type = "character"),

#   make_option("--selection_trait_category", type = "character", default = "outcome",
#               help = paste("Which trait category drives module selection AND row order:",
#                            "exposure or outcome [default = outcome]. The other category's",
#                            "traits are still plotted in their own dot plot and in the combined",
#                            "plot -- this only controls which category's significance decides",
#                            "which modules are kept and how rows are ranked.")),

#   make_option("--drop_grey_modules", type = "logical", default = TRUE),
#   make_option("--p_thresh", type = "double", default = 0.05),
#   make_option("--min_timepoints_replicated", type = "integer", default = 2),
#   make_option("--require_consistent_direction", type = "logical", default = TRUE),
#   make_option("--min_exposure_traits", type = "integer", default = 1),
#   make_option("--drop_empty_columns", type = "logical", default = FALSE),

#   make_option("--color_limit", type = "double", default = NA),
#   make_option("--dot_size_min", type = "double", default = 1.5),
#   make_option("--dot_size_max", type = "double", default = 6),

#   make_option("--font_family", type = "character", default = "Helvetica"),
#   make_option("--font_size", type = "double", default = 8),

#   make_option("--full_width_mm", type = "double", default = 170),
#   make_option("--row_height_in", type = "double", default = 0.155),
#   make_option("--header_margin_in", type = "double", default = 1.0),
#   make_option("--max_height_in", type = "double", default = 9.5),
#   make_option("--col_width_in", type = "double", default = 0.16,
#               help = paste("Inches of width reserved per trait column per facet panel,",
#                            "used to auto-widen the figure beyond --full_width_mm when",
#                            "there are enough columns that the journal full-width would",
#                            "make columns/group-labels illegible or overlapping. [default = 0.16]")),
#   make_option("--width_margin_in", type = "double", default = 1.3,
#               help = "Inches reserved for y-axis module labels + legend, outside the facet panels [default = 1.3]"),

#   make_option("--exposure_plot_width", type = "double", default = NA),
#   make_option("--exposure_plot_height", type = "double", default = NA),
#   make_option("--outcome_plot_width", type = "double", default = NA),
#   make_option("--outcome_plot_height", type = "double", default = NA),
#   make_option("--combined_plot_width", type = "double", default = NA),
#   make_option("--combined_plot_height", type = "double", default = NA),

#   make_option("--exposure_group_label", type = "character", default = "Exposures"),
#   make_option("--outcome_group_label", type = "character", default = "Outcomes"),
#   make_option("--add_panel_letters", type = "logical", default = TRUE),

#   make_option("--panel_layout", type = "character", default = "horizontal",
#               help = paste("horizontal or vertical [default = horizontal]. horizontal places",
#                            "timepoint panels side by side (a, b, c across); vertical stacks",
#                            "them (a, b, c down). vertical trades width for height -- useful",
#                            "when a wide trait-column count won't fit journal full-width, since",
#                            "each panel only needs ONE panel's worth of width instead of",
#                            "n_timepoints worth.")),
#   make_option("--per_timepoint_plots", type = "logical", default = TRUE,
#               help = paste("TRUE/FALSE; also render single-panel (one timepoint each) versions",
#                            "of all three plots into per-timepoint folders, in addition to the",
#                            "multi-panel combined version in cross_dataset/ [default = TRUE]")),
#   make_option("--panel_spacing_in", type = "double", default = NA,
#               help = paste("Inches of spacing between stacked/side-by-side panels.",
#                            "If not set, defaults to 0.15in for horizontal layout and",
#                            "0.55in for vertical layout (vertical needs more room for",
#                            "the next panel's facet strip + group labels).")),

#   make_option("--save_preview_png", type = "logical", default = TRUE),
#   make_option("--preview_dpi", type = "double", default = 200)
# )

# opt <- parse_args(OptionParser(option_list = option_list))

# # ------------------------------------------------------------
# # Validate
# # ------------------------------------------------------------
# if (is.null(opt$project_root))        stop("--project_root is required")
# if (is.null(opt$regions_file))        stop("--regions_file is required")
# if (is.null(opt$dataset1_label))      stop("--dataset1_label is required")
# if (is.null(opt$dataset1_stats_file)) stop("--dataset1_stats_file is required")
# if (is.null(opt$exposure_outcome_trait_file)) stop("--exposure_outcome_trait_file is required")

# if (!dir.exists(opt$project_root)) stop("project_root not found: ", opt$project_root)
# if (!file.exists(opt$regions_file)) stop("regions_file not found: ", opt$regions_file)
# if (!file.exists(opt$dataset1_stats_file)) stop("dataset1_stats_file not found: ", opt$dataset1_stats_file)
# if (!file.exists(opt$exposure_outcome_trait_file)) stop("exposure_outcome_trait_file not found: ", opt$exposure_outcome_trait_file)

# if (opt$p_thresh <= 0 || opt$p_thresh >= 1) stop("--p_thresh must be > 0 and < 1")
# if (opt$min_timepoints_replicated < 1) stop("--min_timepoints_replicated must be >= 1")
# selection_trait_category <- tolower(opt$selection_trait_category)
# if (!selection_trait_category %in% c("exposure", "outcome")) {
#   stop("--selection_trait_category must be one of: exposure, outcome")
# }
# panel_layout <- tolower(opt$panel_layout)
# if (!panel_layout %in% c("horizontal", "vertical")) {
#   stop("--panel_layout must be one of: horizontal, vertical")
# }
# panel_spacing_in <- opt$panel_spacing_in
# if (is.na(panel_spacing_in)) {
#   panel_spacing_in <- if (panel_layout == "vertical") 0.55 else 0.15
# }
# if (panel_spacing_in < 0) stop("--panel_spacing_in must be >= 0")
# if (opt$min_exposure_traits < 1) stop("--min_exposure_traits must be >= 1")
# if (opt$dot_size_min <= 0 || opt$dot_size_max <= opt$dot_size_min) {
#   stop("--dot_size_max must be > --dot_size_min > 0")
# }

# dataset_inputs <- list(
#   list(label = opt$dataset1_label, stats_file = opt$dataset1_stats_file)
# )

# for (i in 2:3) {
#   lbl <- opt[[paste0("dataset", i, "_label")]]
#   sf  <- opt[[paste0("dataset", i, "_stats_file")]]

#   if (!is.null(lbl) || !is.null(sf)) {
#     if (is.null(lbl) || is.null(sf)) {
#       stop("Provide both --dataset", i, "_label and --dataset", i, "_stats_file")
#     }
#     if (!file.exists(sf)) {
#       stop("dataset", i, "_stats_file not found: ", sf)
#     }
#     dataset_inputs[[length(dataset_inputs) + 1]] <- list(label = lbl, stats_file = sf)
#   }
# }

# dataset_labels <- sapply(dataset_inputs, `[[`, "label")

# if (opt$min_timepoints_replicated > length(dataset_labels)) {
#   stop("--min_timepoints_replicated (", opt$min_timepoints_replicated,
#        ") cannot exceed the number of datasets provided (", length(dataset_labels), ")")
# }

# # ------------------------------------------------------------
# # Output dirs
# # ------------------------------------------------------------
# {
#   rfile_dir    <- dirname(opt$regions_file)
#   region_label <- basename(rfile_dir)
#   cpg_label    <- basename(dirname(rfile_dir))
#   message("Derived cpg_label    : ", cpg_label)
#   message("Derived region_label : ", region_label)
# }

# pipeline_root <- file.path(opt$project_root, "comethyl_output", "consensus")
# step_dir      <- file.path(pipeline_root, "12c_me_trait_exposure_first_dotplots",
#                             cpg_label, region_label, opt$adjustment_version)
# cross_dir     <- file.path(step_dir, "cross_dataset")
# dir.create(step_dir, recursive = TRUE, showWarnings = FALSE)
# dir.create(cross_dir, recursive = TRUE, showWarnings = FALSE)
# message("Output root: ", step_dir)
# message("  Shared/combined (all timepoints together): ", cross_dir)
# for (lbl in dataset_labels) {
#   dir.create(file.path(step_dir, lbl), recursive = TRUE, showWarnings = FALSE)
#   message("  Per-timepoint: ", file.path(step_dir, lbl))
# }

# # ------------------------------------------------------------
# # Load exposure/outcome trait definitions
# # ------------------------------------------------------------
# trait_def <- read_exposure_outcome_file(opt$exposure_outcome_trait_file)
# exposure_traits <- unique(trait_def$trait[trait_def$category == "exposure"])
# outcome_traits  <- unique(trait_def$trait[trait_def$category == "outcome"])

# if (length(exposure_traits) == 0) stop("No exposure traits defined in exposure_outcome_trait_file.")
# if (length(outcome_traits) == 0)  stop("No outcome traits defined in exposure_outcome_trait_file.")

# message("Exposure traits: ", length(exposure_traits))
# message("Outcome traits : ", length(outcome_traits))

# # ------------------------------------------------------------
# # Timepoint-restricted traits (optional third "timepoint" column in the
# # trait file). trait_timepoint is a named vector: trait -> dataset label
# # it's restricted to, or NA if unrestricted (applies to every panel).
# # ------------------------------------------------------------
# trait_timepoint <- setNames(trait_def$timepoint, trait_def$trait)

# restricted_traits <- trait_def[!is.na(trait_def$timepoint), , drop = FALSE]
# if (nrow(restricted_traits) > 0) {
#   bad_timepoints <- setdiff(unique(restricted_traits$timepoint), dataset_labels)
#   if (length(bad_timepoints) > 0) {
#     warning(
#       "exposure_outcome_trait_file has timepoint value(s) that don't match any ",
#       "--dataset*_label provided to this run: ", paste(bad_timepoints, collapse = ", "),
#       ". Traits restricted to a timepoint not in this run will never be plotted. ",
#       "(dataset labels this run knows about: ", paste(dataset_labels, collapse = ", "), ")"
#     )
#   }
#   message("Timepoint-restricted traits: ", nrow(restricted_traits), " (",
#           paste(paste0(restricted_traits$trait, "->", restricted_traits$timepoint), collapse = ", "), ")")
# }

# # Returns the subset of `traits` valid for dataset ds_label -- i.e. traits
# # with no timepoint restriction, or restricted specifically to ds_label.
# # Used for per-timepoint plots, where an out-of-timepoint trait's column
# # should be dropped entirely rather than just shown empty.
# filter_traits_for_dataset <- function(traits, ds_label, trait_timepoint) {
#   restrict_val <- trait_timepoint[traits]
#   traits[is.na(restrict_val) | restrict_val == ds_label]
# }

# # ------------------------------------------------------------
# # Load stats tables
# # ------------------------------------------------------------
# stats_list <- list()

# for (ds in dataset_inputs) {
#   df <- readr::read_tsv(ds$stats_file, show_col_types = FALSE)
#   df <- validate_stats_table(df, paste0(ds$label, " stats_file"))

#   cor_col <- get_cor_col(df)
#   df$r_value <- suppressWarnings(as.numeric(df[[cor_col]]))
#   df$p       <- suppressWarnings(as.numeric(df$p))
#   df$module  <- as.character(df$module)
#   df$trait   <- as.character(df$trait)

#   if (isTRUE(opt$drop_grey_modules)) {
#     n_before <- dplyr::n_distinct(df$module)
#     df <- df[!is_grey_module(df$module), ]
#     n_after <- dplyr::n_distinct(df$module)
#     message("[", ds$label, "] Dropped grey module if present: ", n_before, " -> ", n_after, " module(s)")
#   }

#   stats_list[[ds$label]] <- df[, c("module", "trait", "p", "r_value")]

#   message("Loaded stats for: ", ds$label, " (", nrow(df), " rows, ",
#           length(unique(df$module)), " modules, ", length(unique(df$trait)), " traits)")
# }

# # ------------------------------------------------------------
# # Compute robust module list (the core selection step).
# # Category is configurable via --selection_trait_category:
# # defaults to "exposure" but can be switched to "outcome".
# # ------------------------------------------------------------
# candidate_traits <- if (selection_trait_category == "exposure") exposure_traits else outcome_traits
# message("\nSelection category: ", selection_trait_category,
#         " (", length(candidate_traits), " candidate trait(s))")

# pair_summary <- compute_pair_summary(
#   stats_list             = stats_list,
#   dataset_labels          = dataset_labels,
#   candidate_traits        = candidate_traits,
#   p_thresh                = opt$p_thresh,
#   min_timepoints          = opt$min_timepoints_replicated,
#   require_consistent_direction = isTRUE(opt$require_consistent_direction)
# )

# readr::write_tsv(pair_summary, file.path(step_dir, "exposure_outcome_pair_summary.tsv"))

# robust_pairs <- pair_summary[pair_summary$robust %in% TRUE, ]

# if (nrow(robust_pairs) == 0) {
#   stop(
#     "No module x ", selection_trait_category, "-trait pairs met the robustness criterion ",
#     "(p < ", opt$p_thresh, " in >= ", opt$min_timepoints_replicated, " timepoint(s)",
#     if (isTRUE(opt$require_consistent_direction)) ", consistent direction" else "",
#     "). Try lowering --min_timepoints_replicated, raising --p_thresh, setting ",
#     "--require_consistent_direction FALSE, or switching --selection_trait_category."
#   )
# }

# module_summary <- robust_pairs %>%
#   dplyr::group_by(module) %>%
#   dplyr::summarise(
#     n_robust_traits = dplyr::n(),
#     robust_traits = paste(sort(unique(trait)), collapse = ";"),
#     primary_trait = trait[which.min(match(trait, candidate_traits))],
#     primary_trait_min_p = min_p[which.min(match(trait, candidate_traits))],
#     .groups = "drop"
#   ) %>%
#   dplyr::filter(n_robust_traits >= opt$min_exposure_traits) %>%
#   dplyr::mutate(primary_trait_rank = match(primary_trait, candidate_traits)) %>%
#   dplyr::arrange(primary_trait_rank, module)

# if (nrow(module_summary) == 0) {
#   stop("No modules met --min_exposure_traits (", opt$min_exposure_traits, ") for ",
#        "selection category '", selection_trait_category, "'. Lower this threshold and rerun.")
# }

# modules_ordered <- module_summary$module
# n_modules <- length(modules_ordered)

# message("\n", tools::toTitleCase(selection_trait_category), "-robust modules selected: ", n_modules)
# message("(criterion: p < ", opt$p_thresh, " in >= ", opt$min_timepoints_replicated,
#         " of ", length(dataset_labels), " timepoints",
#         if (isTRUE(opt$require_consistent_direction)) ", consistent direction)" else ")")

# write_vector_file(modules_ordered, file.path(step_dir, paste0(selection_trait_category, "_robust_modules.txt")))
# readr::write_tsv(module_summary, file.path(step_dir, paste0(selection_trait_category, "_robust_module_summary.tsv")))


# # ------------------------------------------------------------
# # Dimension planning
# # ------------------------------------------------------------
# journal_width_in <- opt$full_width_mm / 25.4

# message("\nPanel layout: ", panel_layout, " (panel_spacing = ", panel_spacing_in, " in)")
# message("Journal full-width constraint: ", round(journal_width_in, 2),
#         " in (", opt$full_width_mm, " mm, npj full-width column)")

# resolve_dim <- function(user_val, auto_val) if (is.na(user_val)) auto_val else user_val

# # ------------------------------------------------------------
# # Generic renderer, used for both:
# #   (a) the multi-panel "cross_dataset" version (all timepoints together)
# #   (b) each single-panel per-timepoint version
# # so the sizing/skip/save logic only needs to exist once.
# # ------------------------------------------------------------
# emit_dot_plot <- function(dataset_labels_subset, traits_ordered, group_blocks,
#                            title, out_dir, file_stub, user_width, user_height,
#                            label_for_msg, trait_timepoint_restrict = NULL) {

#   plot_df <- build_plot_df(stats_list, dataset_labels_subset, modules_ordered, traits_ordered,
#                             opt$p_thresh, trait_timepoint_restrict = trait_timepoint_restrict)
#   readr::write_tsv(plot_df, file.path(out_dir, paste0(file_stub, "_data.tsv")))

#   if (nrow(plot_df) == 0) {
#     message(label_for_msg, ": skipped -- no nominally significant associations (p < ",
#             opt$p_thresh, ") for the ", n_modules, " selected module(s) in ",
#             paste(dataset_labels_subset, collapse = ", "), ".")
#     writeLines(
#       paste0("Skipped: no significant associations for the ", n_modules,
#              " module(s) selected via --selection_trait_category ", selection_trait_category,
#              " in ", paste(dataset_labels_subset, collapse = ", ")),
#       con = file.path(out_dir, paste0(file_stub, "_SKIPPED.txt"))
#     )
#     return(invisible(NULL))
#   }

#   n_facets_here <- length(dataset_labels_subset)

#   h <- compute_recommended_height(
#     n_modules        = n_modules,
#     row_height_in    = opt$row_height_in,
#     header_margin_in = opt$header_margin_in,
#     max_height_in    = opt$max_height_in,
#     n_facets         = n_facets_here,
#     panel_layout     = panel_layout,
#     panel_spacing_in = panel_spacing_in
#   )
#   w <- compute_recommended_width(
#     n_traits          = length(traits_ordered),
#     n_facets          = n_facets_here,
#     col_width_in      = opt$col_width_in,
#     width_margin_in   = opt$width_margin_in,
#     journal_width_in  = journal_width_in,
#     panel_layout      = panel_layout
#   )

#   width_in  <- resolve_dim(user_width, w$width_in)
#   height_in <- resolve_dim(user_height, h$height_in)

#   message(label_for_msg, ": ", length(traits_ordered), " trait column(s) x ", n_facets_here,
#           " panel(s) [", panel_layout, "] -> ", round(width_in, 2), " x ", round(height_in, 2), " in",
#           if (isTRUE(w$exceeds_journal_width) && panel_layout == "horizontal") paste0(
#             " (EXCEEDS the ", round(journal_width_in, 2),
#             " in journal full-width constraint -- consider --panel_layout vertical, ",
#             "--drop_empty_columns TRUE, or a smaller trait set)"
#           ) else if (isTRUE(h$exceeds_page)) paste0(
#             " (exceeds --max_height_in ", opt$max_height_in, " in)"
#           ) else "")

#   p <- make_dot_plot(
#     plot_df           = plot_df,
#     modules_ordered   = modules_ordered,
#     traits_ordered    = traits_ordered,
#     dataset_labels    = dataset_labels_subset,
#     group_blocks      = group_blocks,
#     color_limit       = opt$color_limit,
#     size_range        = size_range,
#     font_family       = opt$font_family,
#     font_size         = opt$font_size,
#     title             = title,
#     subtitle          = paste0("(modules shown: ", selection_trait_category, "-robust hits [>= ",
#                                 opt$min_timepoints_replicated, " timepoints]; dot plotted only ",
#                                 "where nominally significant, p<", opt$p_thresh, ")"),
#     add_panel_letters = isTRUE(opt$add_panel_letters) && n_facets_here > 1,
#     panel_layout      = panel_layout,
#     panel_spacing_in  = panel_spacing_in
#   )

#   save_plot(p, file.path(out_dir, file_stub), width_in, height_in, opt$font_family,
#             isTRUE(opt$save_preview_png), opt$preview_dpi)
#   message("Saved: ", file_stub, ".pdf (", round(width_in, 2), " x ", round(height_in, 2),
#           " in) -> ", out_dir)

#   invisible(list(width_in = width_in, height_in = height_in))
# }

# size_range <- c(opt$dot_size_min, opt$dot_size_max)

# # ------------------------------------------------------------
# # Resolve trait columns per category (optionally dropping empty ones)
# # ------------------------------------------------------------
# exposure_traits_plot <- exposure_traits
# if (isTRUE(opt$drop_empty_columns)) {
#   tmp <- build_plot_df(stats_list, dataset_labels, modules_ordered, exposure_traits, opt$p_thresh)
#   exposure_traits_plot <- exposure_traits[exposure_traits %in% unique(as.character(tmp$trait))]
#   if (length(exposure_traits_plot) == 0) exposure_traits_plot <- exposure_traits
# }

# outcome_traits_plot <- outcome_traits
# if (isTRUE(opt$drop_empty_columns)) {
#   tmp <- build_plot_df(stats_list, dataset_labels, modules_ordered, outcome_traits, opt$p_thresh)
#   outcome_traits_plot <- outcome_traits[outcome_traits %in% unique(as.character(tmp$trait))]
#   if (length(outcome_traits_plot) == 0) outcome_traits_plot <- outcome_traits
# }

# combined_traits <- c(outcome_traits_plot, exposure_traits_plot)

# plot_specs <- list(
#   Exposures = list(
#     traits       = exposure_traits_plot,
#     group_blocks = list(list(label = opt$exposure_group_label, start = 1, end = length(exposure_traits_plot))),
#     title        = "Consensus methylation module associations with exposures",
#     user_width   = opt$exposure_plot_width,
#     user_height  = opt$exposure_plot_height
#   ),
#   Outcomes = list(
#     traits       = outcome_traits_plot,
#     group_blocks = list(list(label = opt$outcome_group_label, start = 1, end = length(outcome_traits_plot))),
#     title        = "Consensus methylation module associations with pregnancy outcomes",
#     user_width   = opt$outcome_plot_width,
#     user_height  = opt$outcome_plot_height
#   ),
#   Combined = list(
#     traits       = combined_traits,
#     group_blocks = list(
#       list(label = opt$outcome_group_label,  start = 1, end = length(outcome_traits_plot)),
#       list(label = opt$exposure_group_label, start = length(outcome_traits_plot) + 1, end = length(combined_traits))
#     ),
#     title        = "Consensus methylation module associations with outcomes and exposures",
#     user_width   = opt$combined_plot_width,
#     user_height  = opt$combined_plot_height
#   )
# )

# # ------------------------------------------------------------
# # Render: shared/combined (cross_dataset/) version of each plot --
# # all timepoints together, layout controlled by --panel_layout.
# # Per-plot user width/height overrides (--exposure_plot_width etc.) apply
# # ONLY to this combined version; per-timepoint versions below are always
# # auto-sized single panels.
# # ------------------------------------------------------------
# message("\n--- Shared/combined plots (cross_dataset/) ---")
# for (plot_name in names(plot_specs)) {
#   spec <- plot_specs[[plot_name]]
#   emit_dot_plot(
#     dataset_labels_subset     = dataset_labels,
#     traits_ordered            = spec$traits,
#     group_blocks              = spec$group_blocks,
#     title                     = spec$title,
#     out_dir                   = cross_dir,
#     file_stub                 = paste0("Dotplot_", plot_name),
#     user_width                = spec$user_width,
#     user_height                = spec$user_height,
#     label_for_msg             = paste0(plot_name, " (combined)"),
#     trait_timepoint_restrict  = trait_timepoint
#   )
# }

# # ------------------------------------------------------------
# # Render: per-timepoint (single-panel) version of each plot, one folder
# # per dataset label -- matches the per-dataset folder convention used
# # elsewhere in this project (e.g. script 12b).
# # ------------------------------------------------------------
# if (isTRUE(opt$per_timepoint_plots)) {
#   message("\n--- Per-timepoint plots ---")
#   for (ds_label in dataset_labels) {
#     ds_dir <- file.path(step_dir, ds_label)

#     # Drop (not just blank) columns that don't apply to this timepoint --
#     # a single-panel plot has no cross-panel alignment to preserve, so
#     # there's no reason to keep an always-empty column around.
#     exposure_here <- filter_traits_for_dataset(exposure_traits_plot, ds_label, trait_timepoint)
#     outcome_here  <- filter_traits_for_dataset(outcome_traits_plot, ds_label, trait_timepoint)
#     combined_here <- c(outcome_here, exposure_here)

#     n_dropped <- length(exposure_traits_plot) - length(exposure_here)
#     if (n_dropped > 0) {
#       message("[", ds_label, "] dropped ", n_dropped, " exposure column(s) restricted to a different timepoint")
#     }

#     ds_plot_specs <- list(
#       Exposures = list(
#         traits       = exposure_here,
#         group_blocks = list(list(label = opt$exposure_group_label, start = 1, end = length(exposure_here))),
#         title        = plot_specs$Exposures$title
#       ),
#       Outcomes = list(
#         traits       = outcome_here,
#         group_blocks = list(list(label = opt$outcome_group_label, start = 1, end = length(outcome_here))),
#         title        = plot_specs$Outcomes$title
#       ),
#       Combined = list(
#         traits       = combined_here,
#         group_blocks = list(
#           list(label = opt$outcome_group_label,  start = 1, end = length(outcome_here)),
#           list(label = opt$exposure_group_label, start = length(outcome_here) + 1, end = length(combined_here))
#         ),
#         title        = plot_specs$Combined$title
#       )
#     )

#     for (plot_name in names(ds_plot_specs)) {
#       spec <- ds_plot_specs[[plot_name]]
#       emit_dot_plot(
#         dataset_labels_subset    = ds_label,
#         traits_ordered           = spec$traits,
#         group_blocks             = spec$group_blocks,
#         title                    = spec$title,
#         out_dir                  = ds_dir,
#         file_stub                = paste0("Dotplot_", plot_name, "_", ds_label),
#         user_width               = NA,
#         user_height              = NA,
#         label_for_msg            = paste0(plot_name, " (", ds_label, ")"),
#         trait_timepoint_restrict = NULL # already column-filtered above; no row-level restriction needed
#       )
#     }
#   }
# } else {
#   message("\n--per_timepoint_plots FALSE -- skipping per-timepoint folders.")
# }

# # ------------------------------------------------------------
# # Logs
# # ------------------------------------------------------------
# write_log_lines(
#   c(
#     paste("project_root:", opt$project_root),
#     paste("adjustment_version:", opt$adjustment_version),
#     paste("datasets:", paste(dataset_labels, collapse = ", ")),
#     paste("exposure_outcome_trait_file:", opt$exposure_outcome_trait_file),
#     paste("n_exposure_traits_defined:", length(exposure_traits)),
#     paste("n_outcome_traits_defined:", length(outcome_traits)),
#     paste("drop_grey_modules:", opt$drop_grey_modules),
#     paste("p_thresh:", opt$p_thresh),
#     paste("min_timepoints_replicated:", opt$min_timepoints_replicated),
#     paste("require_consistent_direction:", opt$require_consistent_direction),
#     paste("selection_trait_category:", selection_trait_category),
#     paste("min_exposure_traits:", opt$min_exposure_traits),
#     paste("drop_empty_columns:", opt$drop_empty_columns),
#     paste("n_modules_selected:", n_modules),
#     paste("panel_layout:", panel_layout),
#     paste("panel_spacing_in:", panel_spacing_in),
#     paste("per_timepoint_plots:", opt$per_timepoint_plots),
#     paste("color_limit:", ifelse(is.na(opt$color_limit), "auto", opt$color_limit)),
#     paste("font_family:", opt$font_family),
#     paste("font_size:", opt$font_size),
#     paste("full_width_mm:", opt$full_width_mm),
#     paste("row_height_in:", opt$row_height_in),
#     paste("date:", as.character(Sys.time()))
#   ),
#   file.path(step_dir, "run_parameters.txt")
# )

# writeLines(
#   capture.output(sessionInfo()),
#   con = file.path(step_dir, "sessionInfo.txt")
# )

# message("\n\u2713 Script 12C complete: exposure-first presentation dot plots finished")
# message("Outputs saved under: ", step_dir)
# message("  Shared/combined: ", cross_dir)
# if (isTRUE(opt$per_timepoint_plots)) {
#   message("  Per-timepoint: ", paste(file.path(step_dir, dataset_labels), collapse = ", "))
# }
