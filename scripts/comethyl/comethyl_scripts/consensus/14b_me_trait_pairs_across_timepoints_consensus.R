#!/usr/bin/env Rscript

# ================================================================
# SCRIPT 14B: Consensus ME-trait pair plots across timepoints
#
# Pipeline: comethyl WGBS consensus analysis
#
# PURPOSE
#   One figure per module x trait pair with one panel per timepoint
#   (e.g. lightcyan1 vs PPWR 6 mo at 10-16 wk | 36-38 wk | 3 mo PP).
#   Continuous traits: scatter + linear fit. Binary/categorical traits:
#   violin + box + points (No/Yes). Each panel is labelled with that
#   timepoint's bicor and p from the 12a stats file.
#   Script 14 makes the same plots one timepoint at a time.
#
# WHICH PAIRS (combine as needed; duplicates removed)
#   --pairs "lightcyan1:ppwr_6m,lightcyan1:ppwr_e_6m"   explicit pairs
#   --pair_file pairs.tsv                    columns module, trait
#   --modules lightcyan1,coral3 --traits ppwr_6m,ppwr_3m   all combinations
#   --auto TRUE  every pair with p < --p_thresh in >= --min_timepoints
#                datasets, restricted to --traits / --modules if given
#
# INPUTS (as script 14)
#   --project_root, --consensus_modules_rds, --regions_file, --adjustment_version
#   --datasetN_label / --datasetN_sample_info / --datasetN_stats_file (N = 1..3)
#   --sample_id_col (optional), --trait_label_file (optional TSV trait, label;
#   e.g. ppwr_trait_groups.tsv), --dataset_long_labels "10-16 wk,36-38 wk,3 mo PP"
#   gene labels in titles (optional): --gene_assoc_file, --kme_source, --kme_file
#
# OPTIONS
#   --shared_y TRUE        same eigengene axis in every panel [TRUE]
#   --cat_plot_style violin|dot  [violin]
#   --panel_width / --panel_height (in) [3.2 / 3.4]
#
# OUTPUTS
#   comethyl_output/consensus/14b_me_trait_pairs_across_timepoints/<cpg>/<region>/<adj>/
#     <module>_vs_<trait>.pdf/.png, pair_list.tsv (r, p per timepoint), run_parameters.txt
# ================================================================

message("Starting Script 14b")

suppressPackageStartupMessages({
  library(optparse)
  library(dplyr)
  library(ggplot2)
})

option_list <- list(
  make_option("--project_root", type = "character"),
  make_option("--consensus_modules_rds", type = "character"),
  make_option("--regions_file", type = "character"),
  make_option("--adjustment_version", type = "character", default = "unadjusted"),
  make_option("--out_dir", type = "character", default = NULL),
  make_option("--dataset1_label", type = "character"),
  make_option("--dataset1_sample_info", type = "character"),
  make_option("--dataset1_stats_file", type = "character"),
  make_option("--dataset2_label", type = "character", default = NULL),
  make_option("--dataset2_sample_info", type = "character", default = NULL),
  make_option("--dataset2_stats_file", type = "character", default = NULL),
  make_option("--dataset3_label", type = "character", default = NULL),
  make_option("--dataset3_sample_info", type = "character", default = NULL),
  make_option("--dataset3_stats_file", type = "character", default = NULL),
  make_option("--dataset_long_labels", type = "character", default = NULL,
              help = "Panel titles, comma list in dataset order, e.g. '10-16 wk,36-38 wk,3 mo postpartum'."),
  make_option("--sample_id_col", type = "character", default = NULL),
  make_option("--pairs", type = "character", default = NULL),
  make_option("--pair_file", type = "character", default = NULL),
  make_option("--modules", type = "character", default = NULL),
  make_option("--traits", type = "character", default = NULL),
  make_option("--auto", type = "character", default = "FALSE"),
  make_option("--p_thresh", type = "double", default = 0.05),
  make_option("--min_timepoints", type = "integer", default = 2L),
  make_option("--max_pairs", type = "integer", default = 50L),
  make_option("--trait_label_file", type = "character", default = NULL),
  make_option("--shared_y", type = "character", default = "TRUE"),
  make_option("--cat_plot_style", type = "character", default = "violin"),
  make_option("--no_color", type = "character", default = "#4C72B0"),
  make_option("--yes_color", type = "character", default = "#C44E52"),
  make_option("--gene_assoc_file", type = "character", default = NULL),
  make_option("--module_gene_file", type = "character", default = NULL),
  make_option("--kme_file", type = "character", default = NULL),
  make_option("--kme_source", type = "character", default = "primary"),
  make_option("--kme_dir", type = "character", default = NULL),
  make_option("--top_n_genes", type = "integer", default = 2L),
  make_option("--panel_width", type = "double", default = 3.2),
  make_option("--panel_height", type = "double", default = 3.4),
  make_option("--font_family", type = "character", default = "Helvetica"),
  make_option("--font_size", type = "double", default = 11),
  make_option("--helper_file", type = "character", default = NULL)
)
opt <- parse_args(OptionParser(option_list = option_list))
tb <- function(x) toupper(trimws(as.character(x))) %in% c("TRUE", "T", "1", "YES")
opt$auto <- tb(opt$auto); opt$shared_y <- tb(opt$shared_y)
opt$cat_plot_style <- match.arg(tolower(opt$cat_plot_style), c("violin", "dot"))

script_dir <- tryCatch({
  f <- sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE))
  if (length(f) == 1) dirname(normalizePath(f)) else getwd()
}, error = function(e) getwd())
helper_file <- if (!is.null(opt$helper_file)) opt$helper_file else file.path(script_dir, "helper.R")
if (!file.exists(helper_file)) stop("helper.R not found at: ", helper_file, call. = FALSE)
source(helper_file)

# ------------------------------------------------------------
# Inputs
# ------------------------------------------------------------
for (nm in c("project_root", "consensus_modules_rds", "dataset1_label", "dataset1_sample_info", "dataset1_stats_file"))
  stop_if_missing(opt[[nm]], paste0("--", nm))
validate_file_exists(opt$consensus_modules_rds, "--consensus_modules_rds")
ds <- list()
for (k in 1:3) {
  lb <- opt[[paste0("dataset", k, "_label")]]; si <- opt[[paste0("dataset", k, "_sample_info")]]; sf <- opt[[paste0("dataset", k, "_stats_file")]]
  if (is.null(lb) && is.null(si) && is.null(sf)) next
  if (is.null(lb) || is.null(si) || is.null(sf)) stop("Provide --dataset", k, "_label, _sample_info and _stats_file together.", call. = FALSE)
  validate_file_exists(si, paste0("--dataset", k, "_sample_info")); validate_file_exists(sf, paste0("--dataset", k, "_stats_file"))
  ds[[lb]] <- list(sample_info = si, stats_file = sf)
}
tp <- names(ds)
tp_long <- if (!is.null(opt$dataset_long_labels)) trimws(strsplit(opt$dataset_long_labels, ",")[[1]]) else tp
if (length(tp_long) != length(tp)) stop("--dataset_long_labels needs one entry per dataset.", call. = FALSE)
names(tp_long) <- tp
if (opt$min_timepoints > length(tp)) opt$min_timepoints <- length(tp)

cpg_label <- region_label <- NULL
if (!is.null(opt$regions_file)) { rd <- dirname(opt$regions_file); region_label <- basename(rd); cpg_label <- basename(dirname(rd)) }
out_dir <- opt$out_dir
if (is.null(out_dir)) {
  if (is.null(cpg_label)) stop("Give --regions_file or --out_dir.", call. = FALSE)
  out_dir <- file.path(opt$project_root, "comethyl_output", "consensus", "14b_me_trait_pairs_across_timepoints",
                       cpg_label, region_label, opt$adjustment_version)
}
safe_dir_create(out_dir); message("Output: ", out_dir)

key <- function(x) gsub("[[:space:][:punct:]]+", "", tolower(trimws(x)))

# consensus MEs per dataset
cm <- readRDS(opt$consensus_modules_rds)
if (is.null(cm$multiMEs)) stop("consensus_modules_rds has no $multiMEs.", call. = FALSE)
MEs <- lapply(tp, function(d) {
  x <- cm$multiMEs[[d]]
  if (is.null(x)) stop("No MEs for dataset '", d, "' in consensus_modules_rds (has: ", paste(names(cm$multiMEs), collapse = ", "), ").", call. = FALSE)
  as.data.frame(if (!is.null(x$data)) x$data else x, check.names = FALSE)
}); names(MEs) <- tp

# sample info + stats per dataset
SI <- lapply(tp, function(d) readSampleInfo(ds[[d]]$sample_info, sample_id_col = opt$sample_id_col, verbose = FALSE)); names(SI) <- tp
ST <- lapply(tp, function(d) {
  s <- utils::read.delim(ds[[d]]$stats_file, stringsAsFactors = FALSE, check.names = FALSE)
  pc <- intersect(c("p", "pvalue", "p.value", "p_value"), colnames(s))[1]
  rc <- intersect(c("bicor", "cor", "correlation"), colnames(s))[1]
  if (is.na(pc) || is.na(rc)) stop(d, " stats file needs p and bicor/cor columns.", call. = FALSE)
  data.frame(module = sub("^ME", "", as.character(s$module)), trait = as.character(s$trait), tkey = key(s$trait),
             r = as.numeric(s[[rc]]), p = as.numeric(s[[pc]]), stringsAsFactors = FALSE)
}); names(ST) <- tp

trait_lab <- NULL
if (!is.null(opt$trait_label_file)) {
  tl <- utils::read.delim(opt$trait_label_file, stringsAsFactors = FALSE, check.names = FALSE, quote = "", comment.char = "#")
  if (all(c("trait", "label") %in% colnames(tl))) trait_lab <- stats::setNames(tl$label, key(tl$trait))
}
lab_of <- function(t) { l <- if (!is.null(trait_lab)) trait_lab[key(t)] else NA; if (is.na(l)) t else unname(l) }

# ------------------------------------------------------------
# Pairs
# ------------------------------------------------------------
split_csv <- function(x) if (is.null(x)) character(0) else trimws(strsplit(x, ",")[[1]])
pairs <- data.frame(module = character(0), trait = character(0), stringsAsFactors = FALSE)
if (!is.null(opt$pairs)) {
  pp <- split_csv(opt$pairs); bad <- pp[!grepl(":", pp, fixed = TRUE)]
  if (length(bad)) stop("--pairs entries must look like module:trait (bad: ", paste(bad, collapse = ", "), ")", call. = FALSE)
  pairs <- rbind(pairs, data.frame(module = sub(":.*$", "", pp), trait = sub("^[^:]*:", "", pp)))
}
if (!is.null(opt$pair_file)) {
  pf <- utils::read.delim(opt$pair_file, stringsAsFactors = FALSE, check.names = FALSE)
  if (!all(c("module", "trait") %in% colnames(pf))) stop("--pair_file needs columns module and trait.", call. = FALSE)
  pairs <- rbind(pairs, pf[, c("module", "trait")])
}
mods_req <- split_csv(opt$modules); trs_req <- split_csv(opt$traits)
if (length(mods_req) && length(trs_req) && !opt$auto) pairs <- rbind(pairs, expand.grid(module = mods_req, trait = trs_req, stringsAsFactors = FALSE))
if (opt$auto) {
  allst <- do.call(rbind, lapply(tp, function(d) cbind(dataset = d, ST[[d]])))
  allst <- allst[tolower(allst$module) != "grey" & !is.na(allst$p) & allst$p < opt$p_thresh, ]
  if (length(mods_req)) allst <- allst[allst$module %in% mods_req, ]
  if (length(trs_req)) allst <- allst[allst$tkey %in% key(trs_req), ]
  ap <- allst %>% dplyr::group_by(module, tkey) %>%
    dplyr::summarise(trait = trait[1], n = dplyr::n_distinct(dataset), best = min(p), .groups = "drop") %>%
    dplyr::filter(n >= opt$min_timepoints) %>% dplyr::arrange(best)
  message("Auto: ", nrow(ap), " pair(s) with p < ", opt$p_thresh, " in >= ", opt$min_timepoints, " timepoints")
  pairs <- rbind(pairs, ap[, c("module", "trait")])
}
pairs$module <- sub("^ME", "", pairs$module)
pairs <- pairs[!duplicated(paste(pairs$module, key(pairs$trait))), , drop = FALSE]
if (!nrow(pairs)) stop("No pairs: give --pairs, --pair_file, --modules + --traits, or --auto TRUE.", call. = FALSE)
if (nrow(pairs) > opt$max_pairs) { message("Note: plotting the first ", opt$max_pairs, " of ", nrow(pairs), " pairs (--max_pairs)."); pairs <- pairs[seq_len(opt$max_pairs), ] }

# optional gene labels for titles
gl <- NULL
if (!is.null(opt$gene_assoc_file) || !is.null(opt$module_gene_file)) {
  kme_use <- opt$kme_file
  ks <- match.arg(tolower(opt$kme_source), c("primary", "mean", "none"))
  if (is.null(kme_use) && ks != "none" && !is.null(cpg_label)) {
    kf <- find_kme_files(opt$project_root, cpg_label, region_label, opt$adjustment_version, if (ks == "primary") tp[1] else tp, kme_dir = opt$kme_dir)
    if (length(kf)) kme_use <- paste(kf, collapse = ",")
  }
  gl <- module_gene_labels(unique(pairs$module), gene_assoc_file = opt$gene_assoc_file, module_gene_file = opt$module_gene_file,
                           kme_file = kme_use, top_n = opt$top_n_genes)
}

# ------------------------------------------------------------
# Plot each pair
# ------------------------------------------------------------
fs <- opt$font_size
fmt_p <- function(p) if (is.na(p)) "NA" else if (p < 0.001) formatC(p, format = "e", digits = 1) else formatC(p, format = "f", digits = 3)
safe <- function(x) gsub("_+", "_", gsub("[/\\?<>\\:*|\"'` .]+", "_", x))
summary_rows <- list()

for (i in seq_len(nrow(pairs))) {
  mod <- pairs$module[i]; tr <- pairs$trait[i]
  long <- list(); stat_lab <- c(); is_cat <- NA
  for (d in tp) {
    me_col <- intersect(c(paste0("ME", mod), mod), colnames(MEs[[d]]))[1]
    tcol <- colnames(SI[[d]])[match(key(tr), key(colnames(SI[[d]])))]
    s <- ST[[d]][ST[[d]]$module == mod & ST[[d]]$tkey == key(tr), , drop = FALSE]
    r <- if (nrow(s)) s$r[1] else NA; p <- if (nrow(s)) s$p[1] else NA
    summary_rows[[length(summary_rows) + 1]] <- data.frame(module = mod, trait = tr, dataset = d, bicor = r, p = p,
                                                           n = NA_integer_, status = NA_character_)
    if (is.na(me_col) || is.na(tcol)) {
      summary_rows[[length(summary_rows)]]$status <- if (is.na(me_col)) "module not in MEs" else "trait not in sample_info"
      stat_lab[d] <- paste0(tp_long[d], "\n(not available)"); next
    }
    common <- intersect(rownames(MEs[[d]]), rownames(SI[[d]]))
    x <- SI[[d]][common, tcol]; yv <- MEs[[d]][common, me_col]
    ok <- !is.na(x) & !is.na(yv)
    summary_rows[[length(summary_rows)]]$n <- sum(ok); summary_rows[[length(summary_rows)]]$status <- "plotted"
    if (sum(ok) < 3) { stat_lab[d] <- paste0(tp_long[d], "\n(n < 3)"); next }
    cat_here <- is.character(x) || is.factor(x) || is_binary_like(x[ok])
    is_cat <- if (is.na(is_cat)) cat_here else (is_cat || cat_here)
    long[[d]] <- data.frame(dataset = d, x = x[ok], ME = yv[ok], stringsAsFactors = FALSE)
    stat_lab[d] <- paste0(tp_long[d], "\nr = ", ifelse(is.na(r), "NA", sprintf("%.2f", r)), ", p = ", fmt_p(p), ", n = ", sum(ok))
  }
  if (!length(long)) { message("Skipped ", mod, " vs ", tr, ": no data at any timepoint."); next }
  dat <- do.call(rbind, long)
  dat$panel <- factor(stat_lab[dat$dataset], levels = stat_lab[tp])
  sig <- stats::setNames(vapply(tp, function(d) { s <- ST[[d]][ST[[d]]$module == mod & ST[[d]]$tkey == key(tr), ]
    nrow(s) > 0 && !is.na(s$p[1]) && s$p[1] < opt$p_thresh }, TRUE), stat_lab[tp])

  gname <- if (!is.null(gl)) gl$genes[match(mod, gl$module)] else NA
  ttl <- paste0(if (!is.na(gname)) paste0(gname, " (", mod, ")") else mod, " vs ", lab_of(tr))
  ylab <- paste0(mod, " module eigengene")

  if (isTRUE(is_cat)) {
    dat$grp <- relabel_binary_to_yesno(if (is.factor(dat$x)) dat$x else as.character(dat$x))
    if (!is.factor(dat$grp)) dat$grp <- factor(dat$grp)
    pal <- if (all(levels(dat$grp) %in% c("No", "Yes"))) c(No = opt$no_color, Yes = opt$yes_color) else NULL
    p <- ggplot(dat, aes(grp, ME, fill = grp, colour = grp))
    if (opt$cat_plot_style == "violin") p <- p + geom_violin(trim = FALSE, alpha = 0.45, linewidth = 0.5) +
      geom_boxplot(width = 0.16, outlier.shape = NA, fill = "white", colour = "black", linewidth = 0.45)
    p <- p + geom_jitter(width = 0.1, height = 0, size = 1.4, alpha = 0.85, colour = "black") +
      labs(x = lab_of(tr), y = ylab) + theme(legend.position = "none")
    if (!is.null(pal)) p <- p + scale_fill_manual(values = pal, drop = FALSE) + scale_colour_manual(values = pal, drop = FALSE)
  } else {
    dat$x <- suppressWarnings(as.numeric(dat$x))
    p <- ggplot(dat, aes(x, ME)) +
      geom_smooth(method = "lm", formula = y ~ x, se = TRUE, colour = "grey25", fill = "grey80", linewidth = 0.6) +
      geom_point(size = 1.6, alpha = 0.85, colour = "#2F4F6F") +
      labs(x = lab_of(tr), y = ylab)
  }
  p <- p + facet_wrap(~panel, nrow = 1, scales = if (opt$shared_y) "fixed" else "free_y") +
    labs(title = ttl) +
    theme_bw(base_size = fs, base_family = opt$font_family) +
    theme(panel.grid.minor = element_blank(), strip.background = element_rect(fill = "grey95", colour = NA),
          strip.text = element_text(size = fs - 1, lineheight = 1.05), plot.title = element_text(face = "bold", size = fs + 1),
          legend.position = "none")
  # bold strips where p < threshold (built inside each device so text is measured by that device)
  build <- function() {
    g <- ggplot2::ggplotGrob(p)
    strips <- grep("^strip-t", g$layout$name)
    for (k in seq_along(strips)) {
      lab <- levels(dat$panel)[k]
      if (!is.na(lab) && isTRUE(sig[lab])) {
        tx <- tryCatch(g$grobs[[strips[k]]]$grobs[[1]]$children[[2]]$children[[1]], error = function(e) NULL)
        if (!is.null(tx)) g$grobs[[strips[k]]]$grobs[[1]]$children[[2]]$children[[1]]$gp$font <- 2L
      }
    }
    g
  }
  w <- 0.9 + opt$panel_width * length(tp); h <- opt$panel_height + 0.9
  stem <- file.path(out_dir, paste0(safe(mod), "_vs_", safe(tr)))
  if (capabilities("cairo")) grDevices::cairo_pdf(paste0(stem, ".pdf"), width = w, height = h, family = opt$font_family) else
    grDevices::pdf(paste0(stem, ".pdf"), width = w, height = h)
  grid::grid.newpage(); grid::grid.draw(build()); invisible(grDevices::dev.off())
  grDevices::png(paste0(stem, ".png"), width = w, height = h, units = "in", res = 300, bg = "white",
                 type = if (capabilities("cairo")) "cairo" else getOption("bitmapType"))
  grid::grid.newpage(); grid::grid.draw(build()); invisible(grDevices::dev.off())
  message("Saved: ", basename(stem), ".pdf/.png")
}

sm <- do.call(rbind, summary_rows)
utils::write.table(sm, file.path(out_dir, "pair_list.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write_lines_safe(paste0(names(opt), "\t", vapply(opt, function(x) if (is.null(x)) "" else paste(x, collapse = ","), "")),
                 file.path(out_dir, "run_parameters.txt"))
message("Script 14b complete: ", length(unique(paste(sm$module, sm$trait))), " pair(s) -> ", out_dir)