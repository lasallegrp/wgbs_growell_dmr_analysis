#!/usr/bin/env Rscript

# ================================================================
# SCRIPT 12D: Outcome-group module dot plot (+ trait correlation panel)
#
# Pipeline: comethyl WGBS consensus analysis
#
# PURPOSE
#   Figure-style summary of which consensus modules are associated with
#   which GROUPS of outcomes, and whether those associations are shared
#   between groups or specific to one (e.g. modules associated with PPWR
#   only, with PPWR and GWG, or with body size only).
#
#   Panel b (main): modules on x, grouped into blocks by the combination
#   of outcome groups they are associated with ("PPWR only", "PPWR and
#   GWG", "GWG only", ...); traits on y. A dot is drawn where the module-
#   trait association is nominally significant (p < --p_thresh) in the
#   primary dataset; fill = bicor, size = p. Dots with a thick outline
#   replicate (p < --p_thresh, any direction by default; see
#   --replication_direction) in >= --replication_min datasets.
#   Module labels show top genes over the module colour, from
#   script 15 / 12a gene assignments (or a manual file).
#
#   Panel a (optional, --sample_info): lower-triangle correlation matrix
#   between the traits (Spearman by default), to show how related the
#   outcome groups are.
#
# REQUIRED INPUTS
#   --project_root
#   --regions_file       Filtered_Regions.txt (cpg/region labels, as in 12c)
#   --dataset1_label / --dataset1_stats_file   (+ dataset2/3 optional)
#   --trait_group_file   TSV/CSV, one row per trait:
#       trait   column name in the stats file (matched ignoring case,
#               spaces and punctuation, as resolveTraits())
#       group   outcome group, e.g. PPWR, GWG, Body size, Pregnancy
#       label   (optional) display label, e.g. "PPWR 6 mo"
#       show    (optional) TRUE/FALSE: draw this trait as a row. Traits
#               with show = FALSE still count when classifying modules
#               (e.g. classify with BMI at every timepoint, show one).
#     Group order and row order follow the file.
#
# KEY OPTIONS
#   --primary_dataset       dataset whose p-values decide dots/blocks
#                           [default = dataset1, e.g. Baseline]
#   --focus_group           group the figure is about [default = first group]
#   --primary_trait         trait drawn in bold and used to order modules
#   --hide_nonfocus_blocks  TRUE = drop blocks without the focus group
#                           (e.g. "GWG only") [default = FALSE]
#   --block_order           explicit comma list of block labels
#   --max_modules_per_block keep the N strongest modules per block [8]
#   --replication_min       datasets for the thick-outline marker [2]
#   --replication_direction any (default) = p < thresh in that many datasets,
#                           either sign; same = sign must also agree
#   --require_replication   TRUE = only replicated associations count [FALSE]
#   --gene_assoc_file       Region_Gene_Associations.tsv (15/12a) or
#                           Annotated_Regions.tsv, for module gene labels
#   --kme_file              script 11 module membership -> hub-gene labels
#   --module_gene_file      manual labels: columns module, genes
#   --top_n_genes           genes per label [2]
#   --highlight_modules     comma list; boxed and coloured
#   --sample_info / --sample_id_col / --cor_method   panel a
#
# OUTPUTS
#   comethyl_output/consensus/12d_outcome_group_dotplot/<cpg>/<region>/<adj>/<primary_dataset>/
#     Panel_b_module_dotplot.pdf/.png, Panel_a_trait_correlation.pdf/.png,
#     module_classification.tsv, plotted_associations.tsv,
#     module_gene_labels.tsv, trait_correlations.tsv, run_parameters.txt
# ================================================================

message("Starting Script 12d")

suppressPackageStartupMessages({
  library(optparse)
  library(dplyr)
  library(ggplot2)
  library(scales)
})

option_list <- list(
  make_option("--project_root", type = "character"),
  make_option("--regions_file", type = "character"),
  make_option("--adjustment_version", type = "character", default = "unadjusted"),
  make_option("--out_dir", type = "character", default = NULL,
              help = "Override the output folder [default: derived like 12c]."),
  make_option("--dataset1_label", type = "character"),
  make_option("--dataset1_stats_file", type = "character"),
  make_option("--dataset2_label", type = "character", default = NULL),
  make_option("--dataset2_stats_file", type = "character", default = NULL),
  make_option("--dataset3_label", type = "character", default = NULL),
  make_option("--dataset3_stats_file", type = "character", default = NULL),
  make_option("--primary_dataset", type = "character", default = NULL),
  make_option("--trait_group_file", type = "character"),
  make_option("--focus_group", type = "character", default = NULL),
  make_option("--primary_trait", type = "character", default = NULL),
  make_option("--hide_nonfocus_blocks", type = "character", default = "FALSE"),
  make_option("--block_order", type = "character", default = NULL),
  make_option("--max_modules_per_block", type = "integer", default = 8L),
  make_option("--p_thresh", type = "double", default = 0.05),
  make_option("--replication_min", type = "integer", default = 2L),
  make_option("--require_replication", type = "character", default = "FALSE"),
  make_option("--replication_direction", type = "character", default = "any",
              help = paste("any = p < --p_thresh in >= --replication_min datasets, whatever the sign;",
                           "same = also the same sign in those datasets [default %default].")),
  make_option("--drop_grey_modules", type = "character", default = "TRUE"),
  make_option("--gene_assoc_file", type = "character", default = NULL),
  make_option("--kme_file", type = "character", default = NULL,
              help = "Explicit kME file(s), comma-separated. Overrides --kme_source."),
  make_option("--kme_source", type = "character", default = "primary",
              help = "Auto-find script 11 kME: primary = --primary_dataset's file, mean = mean over all datasets, none = region counts [default %default]."),
  make_option("--kme_dir", type = "character", default = NULL,
              help = "Root of script 11 output [default: <project_root>/comethyl_output/consensus/11_consensus_membership]."),
  make_option("--module_gene_file", type = "character", default = NULL),
  make_option("--top_n_genes", type = "integer", default = 2L),
  make_option("--show_module_names", type = "character", default = "TRUE",
              help = "Show the module colour name under the gene label [default %default]."),
  make_option("--highlight_modules", type = "character", default = NULL),
  make_option("--highlight_color", type = "character", default = "#1F3F99"),
  make_option("--sample_info", type = "character", default = NULL),
  make_option("--sample_id_col", type = "character", default = NULL),
  make_option("--cor_method", type = "character", default = "spearman"),
  make_option("--color_limit", type = "double", default = NA),
  make_option("--font_family", type = "character", default = "Helvetica"),
  make_option("--font_size", type = "double", default = 8),
  make_option("--plot_width", type = "double", default = NA),
  make_option("--plot_height", type = "double", default = NA),
  make_option("--preview_dpi", type = "double", default = 300),
  make_option("--helper_file", type = "character", default = NULL)
)
opt <- parse_args(OptionParser(option_list = option_list))
opt$replication_direction <- match.arg(tolower(opt$replication_direction), c("any", "same"))
tb <- function(x, nm) {
  v <- toupper(trimws(as.character(x)))
  if (v %in% c("TRUE", "T", "1", "YES")) return(TRUE)
  if (v %in% c("FALSE", "F", "0", "NO")) return(FALSE)
  stop("--", nm, " must be TRUE or FALSE", call. = FALSE)
}
for (nm in c("hide_nonfocus_blocks", "require_replication", "drop_grey_modules", "show_module_names")) opt[[nm]] <- tb(opt[[nm]], nm)

script_dir <- tryCatch({
  f <- sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE))
  if (length(f) == 1) dirname(normalizePath(f)) else getwd()
}, error = function(e) getwd())
helper_file <- if (!is.null(opt$helper_file)) opt$helper_file else file.path(script_dir, "helper.R")
if (!file.exists(helper_file)) stop("helper.R not found at: ", helper_file, call. = FALSE)
source(helper_file)
has_ggtext <- requireNamespace("ggtext", quietly = TRUE)
if (!has_ggtext) message("Note: ggtext not installed; labels fall back to plain text.")

# ------------------------------------------------------------
# Validate inputs
# ------------------------------------------------------------
for (nm in c("project_root", "dataset1_label", "dataset1_stats_file", "trait_group_file")) stop_if_missing(opt[[nm]], paste0("--", nm))
if (is.null(opt$out_dir)) stop_if_missing(opt$regions_file, "--regions_file (or --out_dir)")
datasets <- list()
for (k in 1:3) {
  lb <- opt[[paste0("dataset", k, "_label")]]; sf <- opt[[paste0("dataset", k, "_stats_file")]]
  if (is.null(lb) && is.null(sf)) next
  if (is.null(lb) || is.null(sf)) stop("Provide both --dataset", k, "_label and --dataset", k, "_stats_file", call. = FALSE)
  validate_file_exists(sf, paste0("--dataset", k, "_stats_file"))
  datasets[[lb]] <- sf
}
ds_labels <- names(datasets)
primary <- if (is.null(opt$primary_dataset)) ds_labels[1] else opt$primary_dataset
if (!primary %in% ds_labels) stop("--primary_dataset '", primary, "' is not one of: ", paste(ds_labels, collapse = ", "), call. = FALSE)
if (opt$replication_min > length(ds_labels)) opt$replication_min <- length(ds_labels)

out_dir <- opt$out_dir
if (is.null(out_dir)) {
  rd <- dirname(opt$regions_file)
  out_dir <- file.path(opt$project_root, "comethyl_output", "consensus", "12d_outcome_group_dotplot",
                       basename(dirname(rd)), basename(rd), opt$adjustment_version, primary)
}
safe_dir_create(out_dir)
message("Output: ", out_dir)

# ------------------------------------------------------------
# Trait groups
# ------------------------------------------------------------
tg <- if (tolower(tools::file_ext(opt$trait_group_file)) == "csv") {
  utils::read.csv(opt$trait_group_file, stringsAsFactors = FALSE, check.names = FALSE)
} else utils::read.delim(opt$trait_group_file, stringsAsFactors = FALSE, check.names = FALSE, quote = "")
if (!all(c("trait", "group") %in% colnames(tg))) stop("--trait_group_file needs columns trait and group", call. = FALSE)
tg$trait <- trimws(tg$trait); tg$group <- trimws(tg$group)
tg <- tg[nzchar(tg$trait) & nzchar(tg$group), , drop = FALSE]
if (!"label" %in% colnames(tg)) tg$label <- tg$trait
tg$label[is.na(tg$label) | !nzchar(trimws(tg$label))] <- tg$trait[is.na(tg$label) | !nzchar(trimws(tg$label))]
tg$show <- if ("show" %in% colnames(tg)) toupper(trimws(as.character(tg$show))) %in% c("TRUE", "T", "1", "YES", "") else TRUE
group_order <- unique(tg$group)
focus <- if (is.null(opt$focus_group)) group_order[1] else opt$focus_group
if (!focus %in% group_order) stop("--focus_group '", focus, "' not in trait_group_file groups: ", paste(group_order, collapse = ", "), call. = FALSE)

# ------------------------------------------------------------
# Stats (trait names resolved like resolveTraits())
# ------------------------------------------------------------
read_stats <- function(f, ds) {
  d <- utils::read.delim(f, stringsAsFactors = FALSE, check.names = FALSE)
  pc <- intersect(c("p", "pvalue", "p.value", "p_value"), colnames(d))[1]
  rc <- intersect(c("bicor", "cor", "correlation"), colnames(d))[1]
  if (is.na(pc) || is.na(rc) || !all(c("module", "trait") %in% colnames(d))) {
    stop(ds, " stats file needs module, trait, p and bicor/cor columns", call. = FALSE)
  }
  data.frame(dataset = ds, module = sub("^ME", "", as.character(d$module)), trait = as.character(d$trait),
             r = as.numeric(d[[rc]]), p = as.numeric(d[[pc]]), stringsAsFactors = FALSE)
}
key <- function(x) gsub("[[:space:][:punct:]]+", "", tolower(trimws(x)))
# Match trait names per dataset (resolveTraits rule: ignore case, spaces,
# punctuation) and relabel every dataset to the trait_group_file spelling,
# so "ppwr_6m" in one file and "ppwr.6m" in another become the same trait.
st <- do.call(rbind, lapply(ds_labels, function(ds) {
  d <- read_stats(datasets[[ds]], ds)
  rt <- resolveTraits(tg$trait, unique(d$trait), label = paste0(ds, " traits"), verbose = TRUE)
  if (length(rt$missing) > 0) message("WARNING: ", ds, ": not in stats file: ", paste(rt$missing, collapse = ", "))
  d$trait <- tg$trait[match(key(d$trait), key(tg$trait))]
  d[!is.na(d$trait), , drop = FALSE]
}))
if (opt$drop_grey_modules) st <- st[tolower(st$module) != "grey", , drop = FALSE]
tg$stats_trait <- tg$trait
tg$stats_trait[!tg$trait %in% st$trait] <- NA
tg <- tg[!is.na(tg$stats_trait), , drop = FALSE]
if (nrow(tg) == 0) stop("No trait in --trait_group_file was found in the stats files.", call. = FALSE)
st <- merge(st, tg[, c("stats_trait", "group", "label", "show")], by.x = "trait", by.y = "stats_trait")
primary_trait <- NULL
if (!is.null(opt$primary_trait)) {
  primary_trait <- tg$stats_trait[match(key(opt$primary_trait), key(tg$trait))]
  if (is.na(primary_trait)) stop("--primary_trait '", opt$primary_trait, "' not in trait_group_file", call. = FALSE)
}

# ------------------------------------------------------------
# Significance + replication per module x trait
# ------------------------------------------------------------
st$sig <- !is.na(st$p) & st$p < opt$p_thresh
rep_tab <- st %>% dplyr::group_by(module, trait) %>%
  dplyr::summarise(n_sig_pos = sum(sig & r > 0), n_sig_neg = sum(sig & r < 0),
                   datasets_sig = paste(dataset[sig], collapse = ";"), .groups = "drop") %>%
  dplyr::mutate(n_rep = pmax(n_sig_pos, n_sig_neg))
pr <- merge(st[st$dataset == primary, ], rep_tab, by = c("module", "trait"))
pr$n_sig_any <- pr$n_sig_pos + pr$n_sig_neg
pr$direction_consistent <- pr$n_sig_pos == 0 | pr$n_sig_neg == 0
pr$replicated <- if (opt$replication_direction == "any") {
  pr$sig & pr$n_sig_any >= opt$replication_min
} else {
  pr$sig & ((pr$r > 0 & pr$n_sig_pos >= opt$replication_min) | (pr$r < 0 & pr$n_sig_neg >= opt$replication_min))
}
pr$counts <- if (opt$require_replication) pr$replicated else pr$sig

# ------------------------------------------------------------
# Classify modules by the set of groups they are associated with
# ------------------------------------------------------------
mod_groups <- pr[pr$counts, ] %>% dplyr::group_by(module) %>%
  dplyr::summarise(groups = paste(group_order[group_order %in% unique(group)], collapse = "|"),
                   n_groups = length(unique(group)), .groups = "drop")
if (nrow(mod_groups) == 0) stop("No module has a nominally significant association with any listed trait in ", primary, call. = FALSE)
block_label <- function(g) {
  gs <- strsplit(g, "|", fixed = TRUE)[[1]]
  if (length(gs) == 1) paste(gs, "only") else if (length(gs) == 2) paste(gs, collapse = " and ") else
    paste(gs, collapse = " + ")
}
mod_groups$block <- vapply(mod_groups$groups, block_label, "")
mod_groups$has_focus <- vapply(strsplit(mod_groups$groups, "|", fixed = TRUE), function(g) focus %in% g, TRUE)

# strength for ordering: primary trait p, else best p within focus group, else best p overall
strength <- function(m) {
  d <- pr[pr$module == m & pr$counts, ]
  if (!is.null(primary_trait) && any(d$trait == primary_trait)) return(min(d$p[d$trait == primary_trait]))
  if (any(d$group == focus)) return(min(d$p[d$group == focus]))
  min(d$p)
}
mod_groups$order_p <- vapply(mod_groups$module, strength, 0)
if (opt$hide_nonfocus_blocks) mod_groups <- mod_groups[mod_groups$has_focus, , drop = FALSE]
if (nrow(mod_groups) == 0) stop("No module is associated with the focus group '", focus, "'.", call. = FALSE)

blocks <- unique(mod_groups$block)
if (!is.null(opt$block_order)) {
  bo <- trimws(strsplit(opt$block_order, ",")[[1]])
  unk <- setdiff(bo, blocks); if (length(unk) > 0) message("Note: --block_order labels not present: ", paste(unk, collapse = ", "))
  blocks <- c(intersect(bo, blocks), setdiff(blocks, bo))
} else {
  bi <- mod_groups[!duplicated(mod_groups$block), c("block", "groups", "has_focus", "n_groups")]
  first_g <- vapply(strsplit(bi$groups, "|", fixed = TRUE), function(g) min(match(g, group_order)), 0)
  # focus-only first, then focus shared with others, then blocks without the focus group
  bi$rank <- ifelse(bi$has_focus & bi$n_groups == 1, 0, ifelse(bi$has_focus, 1, 2))
  blocks <- bi$block[order(bi$rank, bi$n_groups, first_g)]
}
mod_groups$block <- factor(mod_groups$block, levels = blocks)
mod_groups <- mod_groups[order(mod_groups$block, mod_groups$order_p), , drop = FALSE]
mod_groups$n_in_block <- ave(seq_len(nrow(mod_groups)), mod_groups$block, FUN = length)
dropped <- mod_groups %>% dplyr::group_by(block) %>% dplyr::slice(-seq_len(opt$max_modules_per_block)) %>% dplyr::ungroup()
kept <- mod_groups %>% dplyr::group_by(block) %>% dplyr::slice_head(n = opt$max_modules_per_block) %>% dplyr::ungroup()
if (nrow(dropped) > 0) message("Note: ", nrow(dropped), " module(s) beyond --max_modules_per_block not shown (listed in module_classification.tsv)")
modules <- kept$module

# ------------------------------------------------------------
# Gene labels
# ------------------------------------------------------------
kme_use <- opt$kme_file
kme_source <- match.arg(tolower(opt$kme_source), c("primary", "mean", "none"))
if (is.null(kme_use) && kme_source != "none") {
  if (is.null(opt$regions_file)) {
    message("No --regions_file: cannot derive cpg/region labels to find kME files; ranking genes by region count.")
  } else {
    rd <- dirname(opt$regions_file)
    kf <- find_kme_files(opt$project_root, basename(dirname(rd)), basename(rd), opt$adjustment_version,
                         if (kme_source == "primary") primary else ds_labels, kme_dir = opt$kme_dir)
    if (length(kf)) kme_use <- paste(kf, collapse = ",") else message("No kME files found; module genes ranked by region count.")
  }
}
gl <- module_gene_labels(mod_groups$module, gene_assoc_file = opt$gene_assoc_file, module_gene_file = opt$module_gene_file,
                         kme_file = kme_use, top_n = opt$top_n_genes)
mod_out <- merge(mod_groups, gl, by = "module", all.x = TRUE)
mod_out$shown <- mod_out$module %in% modules
mod_out <- mod_out[order(mod_out$block, mod_out$order_p), c("module", "genes", "block", "groups", "order_p", "shown", "source")]
utils::write.table(mod_out, file.path(out_dir, "module_classification.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
utils::write.table(gl, file.path(out_dir, "module_gene_labels.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

highlight <- if (is.null(opt$highlight_modules)) character(0) else trimws(strsplit(opt$highlight_modules, ",")[[1]])
# HTML-escape and encode non-ASCII (en dash, rho, ...) as entities so
# ggtext renders them the same in any locale.
esc <- function(x) vapply(enc2utf8(as.character(x)), function(s) {
  s <- gsub("<", "&lt;", gsub("&", "&amp;", s))
  cp <- utf8ToInt(s); if (all(cp < 128)) return(s)
  paste(ifelse(cp < 128, vapply(cp, function(i) intToUtf8(i), ""), paste0("&#", cp, ";")), collapse = "")
}, "", USE.NAMES = FALSE)
xlab <- vapply(modules, function(m) {
  g <- gl$genes[match(m, gl$module)]
  if (has_ggtext) {
    gtxt <- if (g == "intergenic") paste0("<i>", g, "</i>") else paste0("<b><i>", esc(g), "</i></b>")
    if (m %in% highlight) gtxt <- paste0("<span style='color:", opt$highlight_color, "'>", gtxt, "</span>")
    if (opt$show_module_names) paste0(gtxt, "<br><span style='color:", if (m %in% highlight) opt$highlight_color else "#7F7F7F",
                                     "; font-size:", opt$font_size - 1.5, "pt'>", esc(m), "</span>") else gtxt
  } else if (opt$show_module_names) paste0(g, "\n", m) else g
}, "")

# ------------------------------------------------------------
# Panel b
# ------------------------------------------------------------
rows <- tg[tg$show, , drop = FALSE]
rows$y_label <- if (has_ggtext) esc(rows$label) else rows$label
if (!is.null(primary_trait) && has_ggtext) rows$y_label[rows$stats_trait == primary_trait] <- paste0("<b>", rows$y_label[rows$stats_trait == primary_trait], "</b>")
pdat <- pr[pr$module %in% modules & pr$trait %in% rows$stats_trait & pr$sig, , drop = FALSE]
# x positions: every block spans at least MIN_BLOCK_SPAN columns, so a
# block holding one or two modules still has room for its title
MIN_BLOCK_SPAN <- 2.4
kb <- as.character(kept$block); xpos <- numeric(0); blk_rows <- list(); cursor <- 0
for (b in unique(kb)) {
  nm <- sum(kb == b); span <- max(nm, MIN_BLOCK_SPAN)
  xpos <- c(xpos, cursor + (span - nm) / 2 + seq_len(nm) - 0.5)
  blk_rows[[length(blk_rows) + 1]] <- data.frame(block = b, n = nm, start = cursor, end = cursor + span)
  cursor <- cursor + span
}
blkdf <- do.call(rbind, blk_rows); names(xpos) <- modules
pdat$x <- unname(xpos[pdat$module])
pdat$y <- nrow(rows) + 1 - match(pdat$trait, rows$stats_trait)
pdat$neglogp <- -log10(pmax(pdat$p, 1e-300))
pdat$outline <- ifelse(pdat$replicated, "replicated", "single")
utils::write.table(pdat[order(pdat$x, -pdat$y), c("module", "trait", "label", "group", "r", "p", "replicated", "n_sig_any", "direction_consistent", "datasets_sig")],
                   file.path(out_dir, "plotted_associations.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

nb <- table(factor(kept$block, levels = blocks)); nb <- nb[nb > 0]
# wrap each block title to the width of its block (~5.5 characters per module column)
wrap_to <- pmax(9, round(5.2 * (blkdf$end - blkdf$start)))
blk <- data.frame(label = mapply(function(l, w) paste(strwrap(l, width = w), collapse = "\n"), blkdf$block, wrap_to),
                  x = (blkdf$start + blkdf$end) / 2, stringsAsFactors = FALSE)
blk_lines <- max(lengths(strsplit(blk$label, "\n", fixed = TRUE)))
seps <- utils::head(blkdf$end, -1)
ymax <- nrow(rows)
rlim <- if (is.na(opt$color_limit)) { v <- max(abs(pdat$r), na.rm = TRUE); ceiling(v * 10) / 10 } else opt$color_limit
p_brk <- c(0.05, 0.01, 0.001, 1e-4); p_brk <- p_brk[p_brk < opt$p_thresh * 1.0001 & p_brk >= min(pdat$p) / 3]
if (length(p_brk) < 2) p_brk <- c(opt$p_thresh, signif(min(pdat$p), 1))
grp_rows <- rows$group; grp_sep <- which(grp_rows[-1] != grp_rows[-length(grp_rows)])
hsep <- ymax + 1 - (grp_sep + 0.5)
pt_mm <- 25.4 / 72.27

gp <- ggplot(pdat, aes(x = x, y = y)) +
  geom_hline(yintercept = seq_len(ymax), colour = "grey92", linewidth = 0.3) +
  geom_vline(xintercept = xpos, colour = "grey92", linewidth = 0.3) +
  { if (length(hsep)) geom_hline(yintercept = hsep, colour = "grey70", linewidth = 0.3, linetype = "dotted") } +
  { if (length(seps)) geom_vline(xintercept = seps, colour = "grey45", linewidth = 0.4, linetype = "dashed") } +
  geom_point(aes(fill = r, size = neglogp, stroke = outline), shape = 21, colour = "black") +
  scale_discrete_manual("stroke", values = c(replicated = 1.1, single = 0.3), guide = "none") +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0,
                       limits = c(-rlim, rlim), oob = scales::squish, name = "bicor (r)") +
  scale_size_continuous(range = c(1.8, 6), breaks = -log10(p_brk), labels = format(p_brk, scientific = FALSE, drop0trailing = TRUE),
                        limits = c(-log10(opt$p_thresh), max(c(pdat$neglogp, -log10(min(p_brk))))),
                        name = "Nominal p", guide = guide_legend(override.aes = list(fill = "grey75"))) +
  annotate("text", x = blk$x, y = ymax + 0.65, label = blk$label, fontface = "bold", vjust = 0, lineheight = 0.9,
           size = (opt$font_size + 0.5) * pt_mm, family = opt$font_family) +
  scale_x_continuous(breaks = unname(xpos), labels = unname(xlab), limits = c(0, max(blkdf$end)), expand = expansion(add = 0.15)) +
  scale_y_continuous(breaks = seq_len(ymax), labels = rev(rows$y_label), expand = expansion(add = c(0.6, 0.2))) +
  coord_cartesian(ylim = c(0.5, ymax + 0.75 + 0.42 * blk_lines), clip = "off") +
  labs(x = NULL, y = NULL,
       caption = paste0("Dots: p < ", opt$p_thresh, " in ", primary, ". Thick outline: p < ", opt$p_thresh,
                        " in \u2265", opt$replication_min, " of ", length(ds_labels), " timepoints",
                        if (opt$replication_direction == "same") ", same direction." else " (any direction).")) +
  theme_minimal(base_family = opt$font_family, base_size = opt$font_size) +
  theme(panel.grid = element_blank(), legend.title = element_text(size = opt$font_size),
        legend.text = element_text(size = opt$font_size - 1),
        plot.caption = element_text(size = opt$font_size - 1.5, colour = "grey30", hjust = 0),
        plot.margin = margin(14, 6, 4, 6))
gp <- gp + if (has_ggtext) theme(
  axis.text.x = ggtext::element_markdown(angle = 60, hjust = 1, vjust = 1, size = opt$font_size, lineheight = 1.05),
  axis.text.y = ggtext::element_markdown(size = opt$font_size + 0.5)) else theme(
  axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1, size = opt$font_size),
  axis.text.y = element_text(size = opt$font_size + 0.5))
for (m in intersect(highlight, modules)) {
  i <- xpos[[m]]
  gp <- gp + annotate("rect", xmin = i - 0.42, xmax = i + 0.42, ymin = 0.45, ymax = ymax + 0.5,
                      fill = NA, colour = opt$highlight_color, linewidth = 0.9)
}

w <- if (is.na(opt$plot_width)) max(4.5, 2.4 + 0.36 * max(blkdf$end)) else opt$plot_width
h <- if (is.na(opt$plot_height)) max(3.2, 2.3 + 0.33 * ymax + 0.14 * blk_lines) else opt$plot_height
save_fig <- function(p, stem, w, h) {
  if (capabilities("cairo")) grDevices::cairo_pdf(paste0(stem, ".pdf"), width = w, height = h, family = opt$font_family)
  else grDevices::pdf(paste0(stem, ".pdf"), width = w, height = h)
  print(p); grDevices::dev.off()
  ggsave(paste0(stem, ".png"), p, width = w, height = h, dpi = opt$preview_dpi, bg = "white")
  message("Saved: ", stem, ".pdf/.png (", round(w, 2), " x ", round(h, 2), " in)")
}
save_fig(gp, file.path(out_dir, "Panel_b_module_dotplot"), w, h)

# ------------------------------------------------------------
# Panel a: trait-trait correlation (optional)
# ------------------------------------------------------------
if (!is.null(opt$sample_info)) {
  si <- readSampleInfo(opt$sample_info, sample_id_col = opt$sample_id_col, verbose = FALSE)
  ct <- resolveTraits(rows$trait, colnames(si), label = "sample_info traits", verbose = TRUE)
  rr <- rows[key(rows$trait) %in% key(ct$found), , drop = FALSE]
  rr$col <- colnames(si)[match(key(rr$trait), key(colnames(si)))]
  X <- sapply(rr$col, function(cc) suppressWarnings(as.numeric(si[[cc]])))
  k <- ncol(X); res <- list()
  for (i in seq_len(k)) for (j in seq_len(k)) if (j < i) {
    ok <- stats::complete.cases(X[, i], X[, j])
    tt <- suppressWarnings(stats::cor.test(X[ok, i], X[ok, j], method = opt$cor_method, exact = FALSE))
    res[[length(res) + 1]] <- data.frame(row = rr$label[i], col = rr$label[j], rho = unname(tt$estimate),
                                         p = tt$p.value, n = sum(ok), stringsAsFactors = FALSE)
  }
  cd <- do.call(rbind, res)
  utils::write.table(cd, file.path(out_dir, "trait_correlations.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  cd$row <- factor(cd$row, levels = rev(rr$label)); cd$col <- factor(cd$col, levels = rr$label)
  cd$star <- ifelse(cd$p < 0.001, "***", ifelse(cd$p < 0.01, "**", ifelse(cd$p < 0.05, "*", "")))
  cd$txt <- paste0(sprintf("%.2f", cd$rho), ifelse(nzchar(cd$star), paste0("\n", cd$star), ""))
  ga <- ggplot(cd, aes(col, row, fill = rho)) +
    geom_tile(colour = "white", linewidth = 0.6) +
    geom_text(aes(label = txt, colour = abs(rho) > 0.6), size = (opt$font_size - 1.5) * pt_mm, lineheight = 0.85,
              family = opt$font_family, show.legend = FALSE) +
    scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = "grey15")) +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0, limits = c(-1, 1),
                         name = if (opt$cor_method == "spearman") "Spearman ρ" else paste(opt$cor_method, "r")) +
    labs(x = NULL, y = NULL, caption = paste0("* p < 0.05, ** p < 0.01, *** p < 0.001; pairwise complete n = ",
                                               min(cd$n), "–", max(cd$n))) +
    theme_minimal(base_family = opt$font_family, base_size = opt$font_size) +
    theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1),
          plot.caption = element_text(size = opt$font_size - 1.5, colour = "grey30", hjust = 0))
  save_fig(ga, file.path(out_dir, "Panel_a_trait_correlation"), max(4, 1.8 + 0.42 * k), max(3.2, 1.6 + 0.36 * k))
}

write_lines_safe(c(paste0(names(opt), "\t", vapply(opt, function(x) if (is.null(x)) "" else paste(x, collapse = ","), "")),
                   paste0("primary_dataset_used\t", primary), paste0("focus_group_used\t", focus),
                   paste0("kme_files_used\t", if (is.null(kme_use)) "" else kme_use),
                   paste0("blocks\t", paste(names(nb), nb, sep = "=", collapse = "; ")),
                   paste0("n_modules_shown\t", length(modules))), file.path(out_dir, "run_parameters.txt"))
message("Script 12d complete: ", length(modules), " modules in ", length(nb), " blocks -> ", out_dir)