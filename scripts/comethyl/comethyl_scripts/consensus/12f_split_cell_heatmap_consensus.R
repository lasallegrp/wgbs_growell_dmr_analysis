#!/usr/bin/env Rscript

# ================================================================
# SCRIPT 12F: Split-cell heatmap of module-trait associations
#
# Pipeline: comethyl WGBS consensus analysis
#
# PURPOSE
#   One readable figure instead of three timepoint panels. Each
#   module x trait cell is split into one bar per timepoint (e.g.
#   B | 36 | PP). A bar is coloured by bicor where p < --p_thresh, grey
#   where tested but not significant, and dashed where the trait was not
#   measured at that visit. A module replicated across pregnancy reads as
#   a solid coloured cell.
#
#   Rows = modules associated with a "select" column at >= --min_timepoints
#   timepoints, grouped into blocks by the first select group they
#   replicate in (column-file order), ordered by p within block, labelled
#   "genes  module". Columns = traits from --column_file, under group headers.
#
# INPUTS
#   --project_root, --regions_file, --adjustment_version
#   --dataset1_label / --dataset1_stats_file (+ dataset2/3)   12a stats files
#   --column_file  TSV: group, label, select (TRUE/FALSE), then ONE COLUMN
#                  PER DATASET LABEL holding that visit's trait name (blank =
#                  not measured). Lets exposures map to the measure taken at
#                  each draw (epds_base / blank / epds_mth3). Trait names are
#                  matched ignoring case/spaces/punctuation.
#   gene labels:   --gene_assoc_file, --kme_source primary|mean|none (auto kME),
#                  --kme_file, --module_gene_file, --top_n_genes
#
# KEY OPTIONS
#   --p_thresh 0.05, --strong_p 0.01 (white dot), --min_timepoints 2,
#   --replication_direction any|same, --max_rows_per_block 15,
#   --highlight_modules, --color_limit 0.6, --title
#
# THREE VERSIONS (like 12c): Exposures, Outcomes, Combined.
#   Columns are split by the column file's "category" column (exposure /
#   outcome); without it, select = TRUE -> outcome, FALSE -> exposure.
#   --panel_rows shared (default): every version shows the SAME modules,
#     chosen by the select columns (outcomes), in the same order -- rows
#     line up across the three figures, as in 12c.
#   --panel_rows own: each version picks its own modules from its own
#     columns (the exposures figure then shows exposure-replicated modules).
#   --panels Combined,Outcomes,Exposures chooses which versions to draw.
#
# OUTPUTS  comethyl_output/consensus/12f_split_cell_heatmap/<cpg>/<region>/<adj>/
#   Split_cell_heatmap_<Combined|Outcomes|Exposures>.pdf/.png,
#   rows_<version>.tsv, cell_values_<version>.tsv, run_parameters.txt
# ================================================================

message("Starting Script 12f")

suppressPackageStartupMessages({
  library(optparse)
  library(dplyr)
  library(ggplot2)
})

option_list <- list(
  make_option("--project_root", type = "character"),
  make_option("--regions_file", type = "character", default = NULL),
  make_option("--adjustment_version", type = "character", default = "unadjusted"),
  make_option("--out_dir", type = "character", default = NULL),
  make_option("--dataset1_label", type = "character"),
  make_option("--dataset1_stats_file", type = "character"),
  make_option("--dataset2_label", type = "character", default = NULL),
  make_option("--dataset2_stats_file", type = "character", default = NULL),
  make_option("--dataset3_label", type = "character", default = NULL),
  make_option("--dataset3_stats_file", type = "character", default = NULL),
  make_option("--dataset_short_labels", type = "character", default = NULL,
              help = "Comma list of short sub-column labels, e.g. 'B,36,PP' [default: first letters]."),
  make_option("--dataset_long_labels", type = "character", default = NULL,
              help = "Comma list for the key, e.g. '10-16 wk,36-38 wk,3 mo postpartum' [default: dataset labels]."),
  make_option("--column_file", type = "character"),
  make_option("--p_thresh", type = "double", default = 0.05),
  make_option("--strong_p", type = "double", default = 0.01),
  make_option("--min_timepoints", type = "integer", default = 2L),
  make_option("--replication_direction", type = "character", default = "any"),
  make_option("--max_rows_per_block", type = "integer", default = 15L),
  make_option("--panels", type = "character", default = "Combined,Outcomes,Exposures",
              help = "Versions to draw, comma list of Combined, Outcomes, Exposures [default %default]."),
  make_option("--panel_rows", type = "character", default = "shared",
              help = "shared = same modules (chosen by select columns) in every version; own = each version chooses from its own columns [default %default]."),
  make_option("--drop_grey_modules", type = "character", default = "TRUE"),
  make_option("--highlight_modules", type = "character", default = NULL),
  make_option("--highlight_color", type = "character", default = "#1F3F99"),
  make_option("--color_limit", type = "double", default = 0.6),
  make_option("--gene_assoc_file", type = "character", default = NULL),
  make_option("--module_gene_file", type = "character", default = NULL),
  make_option("--kme_file", type = "character", default = NULL),
  make_option("--kme_source", type = "character", default = "primary"),
  make_option("--kme_dir", type = "character", default = NULL),
  make_option("--top_n_genes", type = "integer", default = 2L),
  make_option("--title", type = "character", default = NULL),
  make_option("--font_family", type = "character", default = "Helvetica"),
  make_option("--font_size", type = "double", default = 9),
  make_option("--plot_width", type = "double", default = NA),
  make_option("--plot_height", type = "double", default = NA),
  make_option("--helper_file", type = "character", default = NULL)
)
opt <- parse_args(OptionParser(option_list = option_list))
opt$drop_grey_modules <- toupper(opt$drop_grey_modules) %in% c("TRUE", "T", "1", "YES")
opt$replication_direction <- match.arg(tolower(opt$replication_direction), c("any", "same"))
opt$panel_rows <- match.arg(tolower(opt$panel_rows), c("shared", "own"))
panels <- trimws(strsplit(opt$panels, ",")[[1]])
bad <- setdiff(panels, c("Combined", "Outcomes", "Exposures"))
if (length(bad)) stop("--panels accepts Combined, Outcomes, Exposures (got: ", paste(bad, collapse = ", "), ")", call. = FALSE)

script_dir <- tryCatch({
  f <- sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE))
  if (length(f) == 1) dirname(normalizePath(f)) else getwd()
}, error = function(e) getwd())
helper_file <- if (!is.null(opt$helper_file)) opt$helper_file else file.path(script_dir, "helper.R")
if (!file.exists(helper_file)) stop("helper.R not found at: ", helper_file, call. = FALSE)
source(helper_file)
has_ggtext <- requireNamespace("ggtext", quietly = TRUE)

# ------------------------------------------------------------
# Inputs
# ------------------------------------------------------------
for (nm in c("project_root", "dataset1_label", "dataset1_stats_file", "column_file")) stop_if_missing(opt[[nm]], paste0("--", nm))
datasets <- list()
for (k in 1:3) {
  lb <- opt[[paste0("dataset", k, "_label")]]; sf <- opt[[paste0("dataset", k, "_stats_file")]]
  if (is.null(lb) && is.null(sf)) next
  if (is.null(lb) || is.null(sf)) stop("Provide both --dataset", k, "_label and --dataset", k, "_stats_file", call. = FALSE)
  validate_file_exists(sf, paste0("--dataset", k, "_stats_file")); datasets[[lb]] <- sf
}
tp <- names(datasets); ntp <- length(tp)
tp_short <- if (!is.null(opt$dataset_short_labels)) trimws(strsplit(opt$dataset_short_labels, ",")[[1]]) else substr(tp, 1, 2)
tp_long  <- if (!is.null(opt$dataset_long_labels))  trimws(strsplit(opt$dataset_long_labels, ",")[[1]])  else tp
if (length(tp_short) != ntp || length(tp_long) != ntp) stop("--dataset_short/long_labels need one entry per dataset.", call. = FALSE)
if (opt$min_timepoints > ntp) opt$min_timepoints <- ntp

cpg_label <- region_label <- NULL
if (!is.null(opt$regions_file)) { rd <- dirname(opt$regions_file); region_label <- basename(rd); cpg_label <- basename(dirname(rd)) }
out_dir <- opt$out_dir
if (is.null(out_dir)) {
  if (is.null(cpg_label)) stop("Give --regions_file or --out_dir.", call. = FALSE)
  out_dir <- file.path(opt$project_root, "comethyl_output", "consensus", "12f_split_cell_heatmap", cpg_label, region_label, opt$adjustment_version)
}
safe_dir_create(out_dir); message("Output: ", out_dir)

cf <- utils::read.delim(opt$column_file, stringsAsFactors = FALSE, check.names = FALSE, quote = "", comment.char = "#", na.strings = c("", "NA"))
miss <- setdiff(c("group", "label", tp), colnames(cf))
if (length(miss)) stop("--column_file is missing column(s): ", paste(miss, collapse = ", "),
                       " (needs group, label, select and one column per dataset label: ", paste(tp, collapse = ", "), ")", call. = FALSE)
if (!"select" %in% colnames(cf)) cf$select <- TRUE
cf$select <- toupper(trimws(as.character(cf$select))) %in% c("TRUE", "T", "1", "YES")
cf$category <- if ("category" %in% colnames(cf)) tolower(trimws(cf$category)) else ifelse(cf$select, "outcome", "exposure")
if (!all(cf$category %in% c("exposure", "outcome"))) stop("column_file category must be exposure or outcome.", call. = FALSE)
cf$col_id <- seq_len(nrow(cf))
grp_levels <- unique(cf$group)

key <- function(x) gsub("[[:space:][:punct:]]+", "", tolower(trimws(x)))
read_stats <- function(f, ds) {
  d0 <- utils::read.delim(f, stringsAsFactors = FALSE, check.names = FALSE)
  if ("dataset" %in% colnames(d0)) d0 <- d0[d0$dataset == ds, , drop = FALSE]
  pc <- intersect(c("p", "pvalue", "p.value", "p_value"), colnames(d0))[1]
  rc <- intersect(c("bicor", "cor", "correlation", "r_value"), colnames(d0))[1]
  if (is.na(pc) || is.na(rc) || !all(c("module", "trait") %in% colnames(d0))) stop(ds, " stats file needs module, trait, p and bicor/cor.", call. = FALSE)
  data.frame(dataset = ds, module = sub("^ME", "", as.character(d0$module)), tkey = key(d0$trait),
             r = as.numeric(d0[[rc]]), p = as.numeric(d0[[pc]]), stringsAsFactors = FALSE)
}
st <- do.call(rbind, lapply(tp, function(ds) read_stats(datasets[[ds]], ds)))
if (opt$drop_grey_modules) st <- st[tolower(st$module) != "grey", , drop = FALSE]
modules_all <- sort(unique(st$module))

# cell table: module x column x timepoint
cells <- do.call(rbind, lapply(seq_len(nrow(cf)), function(i) do.call(rbind, lapply(tp, function(ds) {
  tr <- cf[[ds]][i]
  data.frame(col_id = i, dataset = ds, trait = ifelse(is.na(tr), NA_character_, trimws(tr)), stringsAsFactors = FALSE)
}))))
cells <- merge(expand.grid(module = modules_all, col_id = cf$col_id, stringsAsFactors = FALSE), cells, by = "col_id")
cells$tkey <- key(cells$trait)
cells <- merge(cells, st, by = c("module", "dataset", "tkey"), all.x = TRUE)
cells$measured <- !is.na(cells$trait)
present <- unique(paste(st$dataset, st$tkey))                      # trait exists in that timepoint's stats
cells$tested <- cells$measured & paste(cells$dataset, cells$tkey) %in% present
cells$sig <- cells$tested & !is.na(cells$p) & cells$p < opt$p_thresh
nf <- unique(cells[cells$measured & !cells$tested, c("dataset", "trait")])
if (nrow(nf)) message("WARNING: in column_file but not in that timepoint's stats file (drawn dashed): ",
                      paste(paste0(nf$trait, " [", nf$dataset, "]"), collapse = ", "))
cells$measured[cells$measured & !cells$tested] <- FALSE

# ------------------------------------------------------------
# Rows: modules replicated (>= min_timepoints) in the given columns
# ------------------------------------------------------------
compute_rows <- function(col_ids, label) {
  rep <- cells[cells$sig & cells$col_id %in% col_ids, , drop = FALSE]
  if (!nrow(rep)) { message(label, ": no significant associations in its row-selecting columns."); return(NULL) }
  rep <- rep %>% dplyr::group_by(module, col_id) %>%
    dplyr::summarise(n = dplyr::n_distinct(dataset), n_pos = dplyr::n_distinct(dataset[r > 0]), n_neg = dplyr::n_distinct(dataset[r < 0]),
                     best_p = min(p), .groups = "drop") %>%
    dplyr::mutate(ok = if (opt$replication_direction == "any") n >= opt$min_timepoints else pmax(n_pos, n_neg) >= opt$min_timepoints) %>%
    dplyr::filter(ok)
  if (!nrow(rep)) { message(label, ": no module at p < ", opt$p_thresh, " in >= ", opt$min_timepoints, " timepoints."); return(NULL) }
  rep$group <- cf$group[rep$col_id]
  rows <- rep %>% dplyr::group_by(module) %>%
    dplyr::summarise(block = group[order(match(group, grp_levels), best_p)][1],
                     order_p = min(best_p[group == block]), n_rep_cols = dplyr::n(), .groups = "drop") %>%
    dplyr::mutate(block = factor(block, levels = grp_levels)) %>%
    dplyr::arrange(block, order_p) %>%
    dplyr::group_by(block) %>% dplyr::mutate(rank = dplyr::row_number()) %>% dplyr::ungroup()
  rows$shown <- rows$rank <= opt$max_rows_per_block
  if (any(!rows$shown)) message(label, ": ", sum(!rows$shown), " module(s) beyond --max_rows_per_block not drawn (shown = FALSE in rows table).")
  rows
}
sel_cols <- cf$col_id[cf$select]
shared_rows <- compute_rows(sel_cols, "Row selection (select columns)")
panel_cols <- list(Combined = cf$col_id, Outcomes = cf$col_id[cf$category == "outcome"], Exposures = cf$col_id[cf$category == "exposure"])
panel_rows <- lapply(setNames(panels, panels), function(pn) {
  if (!length(panel_cols[[pn]])) return(NULL)
  if (opt$panel_rows == "shared") return(shared_rows)
  own_sel <- if (pn == "Exposures") panel_cols[[pn]] else intersect(panel_cols[[pn]], sel_cols)
  compute_rows(own_sel, pn)
})
all_mods <- unique(unlist(lapply(panel_rows, function(r) if (is.null(r)) NULL else r$module[r$shown])))
if (!length(all_mods)) stop("No module meets the row rule in any version.", call. = FALSE)

# gene labels (once for every module drawn)
kme_use <- opt$kme_file; gl <- NULL
if (!is.null(opt$gene_assoc_file) || !is.null(opt$module_gene_file)) {
  ks <- match.arg(tolower(opt$kme_source), c("primary", "mean", "none"))
  if (is.null(kme_use) && ks != "none" && !is.null(cpg_label)) {
    kf <- find_kme_files(opt$project_root, cpg_label, region_label, opt$adjustment_version,
                         if (ks == "primary") tp[1] else tp, kme_dir = opt$kme_dir)
    if (length(kf)) kme_use <- paste(kf, collapse = ",")
  }
  gl <- module_gene_labels(all_mods, gene_assoc_file = opt$gene_assoc_file, module_gene_file = opt$module_gene_file,
                           kme_file = kme_use, top_n = opt$top_n_genes)
}

highlight <- if (is.null(opt$highlight_modules)) character(0) else trimws(strsplit(opt$highlight_modules, ",")[[1]])
esc <- function(x) gsub("<", "&lt;", gsub("&", "&amp;", x))
fs <- opt$font_size; ptm <- 25.4 / 72.27
lim <- opt$color_limit

# ------------------------------------------------------------
# Draw one version
# ------------------------------------------------------------
draw_version <- function(pn, cols_in, rows_all) {
  cfp <- cf[cf$col_id %in% cols_in, , drop = FALSE]
  rows <- rows_all[rows_all$shown, , drop = FALSE]
  rows$genes <- if (is.null(gl)) NA_character_ else gl$genes[match(rows$module, gl$module)]
  n_hidden <- sum(!rows_all$shown)
  gl_p <- unique(cfp$group)
  cw <- ntp + 0.35
  cfp$x0 <- (seq_len(nrow(cfp)) - 1) * cw + (match(cfp$group, gl_p) - 1) * 1.6
  blocks <- levels(droplevels(rows$block))
  y <- numeric(nrow(rows)); hdr <- data.frame(); cur <- 0
  for (b in blocks) {
    cur <- cur - 1.25; hdr <- rbind(hdr, data.frame(block = b, y = cur + 0.15))
    idx <- which(rows$block == b); y[idx] <- cur - seq_along(idx); cur <- min(y[idx]) - 0.2
  }
  rows$y <- y
  g <- cells[cells$module %in% rows$module & cells$col_id %in% cfp$col_id, , drop = FALSE]
  g$x <- cfp$x0[match(g$col_id, cfp$col_id)] + match(g$dataset, tp) - 1
  g$y <- rows$y[match(g$module, rows$module)]
  g$fill <- pmax(pmin(g$r, lim), -lim)
  ylab <- vapply(seq_len(nrow(rows)), function(i) {
    m <- rows$module[i]; gn <- rows$genes[i]; hl <- m %in% highlight
    if (!has_ggtext) return(if (is.na(gn)) m else paste0(gn, "  ", m))
    gtxt <- if (is.na(gn)) "" else paste0(if (hl) paste0("<b style='color:", opt$highlight_color, "'>") else "", "<i>", esc(gn), "</i>",
                                          if (hl) "</b>" else "", "  ")
    paste0(gtxt, "<span style='color:", if (hl) opt$highlight_color else "#8a8a8a", "'>", esc(m), "</span>")
  }, "")
  top <- max(hdr$y) + 0.9
  xmin_all <- min(cfp$x0) - 0.6; xmax_all <- max(cfp$x0) + ntp - 0.4
  grp_hdr <- cfp %>% dplyr::group_by(group) %>% dplyr::summarise(xmin = min(x0) - 0.4, xmax = max(x0) + ntp - 0.6, .groups = "drop") %>%
    dplyr::mutate(x = (xmin + xmax) / 2)
  y_grp <- top + 0.35 + 0.33 * max(nchar(cfp$label))
  key_x <- xmax_all + 2.2; key_y <- min(rows$y) + 2.5
  row_rule <- if (opt$panel_rows == "shared" || pn != "Exposures") paste(unique(cf$group[cf$select]), collapse = "/") else "exposure"
  ttl <- if (!is.null(opt$title)) paste0(opt$title, if (length(panels) > 1) paste0(" — ", pn) else "") else NULL

  p <- ggplot(g) +
    annotate("text", x = xmin_all, y = hdr$y, label = hdr$block, hjust = 0, vjust = 0.5, fontface = "bold", size = (fs - 0.5) * ptm,
             family = opt$font_family, colour = "#333333") +
    geom_tile(data = g[g$tested, ], aes(x, y), fill = "#EBEBEB", width = 0.9, height = 0.84) +
    geom_tile(data = g[g$sig, ], aes(x, y, fill = fill), width = 0.9, height = 0.84) +
    geom_tile(data = g[!g$measured, ], aes(x, y), fill = NA, colour = "#D6D6D6", linewidth = 0.25, linetype = "22", width = 0.82, height = 0.76) +
    geom_point(data = g[g$sig & g$p < opt$strong_p, ], aes(x, y), colour = "white", size = 0.75) +
    annotate("text", x = cfp$x0 + (ntp - 1) / 2, y = top, label = cfp$label, angle = 50, hjust = 0, vjust = 0, size = (fs - 0.6) * ptm,
             family = opt$font_family) +
    annotate("segment", x = grp_hdr$xmin, xend = grp_hdr$xmax, y = y_grp, yend = y_grp, linewidth = 0.4) +
    annotate("text", x = grp_hdr$x, y = y_grp + 0.4, label = grp_hdr$group, fontface = "bold", size = fs * ptm, vjust = 0,
             family = opt$font_family) +
    annotate("tile", x = key_x + seq_len(ntp) - 1, y = key_y, fill = "#EBEBEB", width = 0.9, height = 0.84) +
    annotate("text", x = key_x + seq_len(ntp) - 1, y = key_y - 0.7, label = tp_short, size = (fs - 2.4) * ptm, vjust = 1, colour = "#444444") +
    annotate("text", x = key_x + (ntp - 1) / 2, y = key_y + 0.7, label = "Each cell", size = (fs - 1.2) * ptm, fontface = "bold", vjust = 0) +
    annotate("text", x = key_x - 0.5, y = key_y - 1.6, hjust = 0, vjust = 1, size = (fs - 2.2) * ptm, colour = "#555555", lineheight = 0.95,
             label = paste(paste0(tp_short, " = ", tp_long), collapse = "\n")) +
    scale_fill_gradient2(low = "#2166AC", mid = "#F7F7F7", high = "#B2182B", limits = c(-lim, lim),
                         name = paste0("bicor (r)\np < ", opt$p_thresh)) +
    scale_y_continuous(breaks = rows$y, labels = ylab, expand = expansion(add = c(0.6, 0))) +
    scale_x_continuous(expand = expansion(add = c(0.4, 3.6))) +
    coord_cartesian(clip = "off", ylim = c(min(rows$y) - 0.6, y_grp + 1.2)) +
    labs(x = NULL, y = NULL, title = ttl,
         caption = paste0("Coloured: p < ", opt$p_thresh, " (white dot: p < ", opt$strong_p, "). Grey: tested, not significant. ",
                          "Dashed: not measured at that visit.\nRows: modules with ", ifelse(grepl("^[AEIOUaeiou]", row_rule), "an ", "a "), row_rule, " association at p < ", opt$p_thresh,
                          " in ≥", opt$min_timepoints, " of ", ntp, " timepoints (", opt$replication_direction,
                          " direction), grouped by the first such group.",
                          if (n_hidden) paste0(" ", n_hidden, " further module(s) not shown.") else "")) +
    theme_minimal(base_size = fs, base_family = opt$font_family) +
    theme(panel.grid = element_blank(), axis.text.x = element_blank(),
          axis.text.y = if (has_ggtext) ggtext::element_markdown(size = fs - 0.3) else element_text(size = fs - 0.3),
          plot.caption = element_text(size = fs - 2.2, colour = "#555555", hjust = 0), plot.caption.position = "plot",
          plot.title = element_text(size = fs + 1.5, face = "bold"), plot.title.position = "plot",
          legend.title = element_text(size = fs - 1), legend.key.height = unit(14, "pt"), plot.margin = margin(4, 8, 4, 4))
  for (m in intersect(highlight, rows$module)) {
    yy <- rows$y[rows$module == m]
    p <- p + annotate("rect", xmin = xmin_all + 0.2, xmax = xmax_all, ymin = yy - 0.5, ymax = yy + 0.5,
                      fill = NA, colour = opt$highlight_color, linewidth = 0.6)
  }
  w <- if (is.na(opt$plot_width)) 4.6 + 0.115 * (max(cfp$x0) + ntp) + 1.4 else opt$plot_width
  h <- if (is.na(opt$plot_height)) 1.4 + 0.045 * max(nchar(cfp$label)) + 0.2 * (nrow(rows) + 1.3 * length(blocks)) else opt$plot_height
  stem <- file.path(out_dir, paste0("Split_cell_heatmap_", pn))
  if (capabilities("cairo")) grDevices::cairo_pdf(paste0(stem, ".pdf"), width = w, height = h, family = opt$font_family) else
    grDevices::pdf(paste0(stem, ".pdf"), width = w, height = h)
  print(p); invisible(grDevices::dev.off())
  ggsave(paste0(stem, ".png"), p, width = w, height = h, dpi = 300, bg = "white")
  message("Saved: ", basename(stem), ".pdf/.png (", round(w, 2), " x ", round(h, 2), " in)")

  ra <- rows_all; ra$genes <- if (is.null(gl)) NA_character_ else gl$genes[match(ra$module, gl$module)]
  utils::write.table(ra[, c("module", "genes", "block", "order_p", "n_rep_cols", "shown")], file.path(out_dir, paste0("rows_", pn, ".tsv")),
                     sep = "\t", quote = FALSE, row.names = FALSE)
  cv <- g[, c("module", "dataset", "trait", "r", "p", "measured", "sig")]
  cv$column <- cfp$label[match(g$col_id, cfp$col_id)]; cv$group <- cfp$group[match(g$col_id, cfp$col_id)]
  ord <- order(match(cv$module, rows$module), match(g$col_id, cfp$col_id), match(cv$dataset, tp))
  utils::write.table(cv[ord, c("module", "group", "column", "dataset", "trait", "r", "p", "measured", "sig")],
                     file.path(out_dir, paste0("cell_values_", pn, ".tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
  invisible(nrow(rows))
}

drawn <- c()
for (pn in panels) {
  if (!length(panel_cols[[pn]])) { message(pn, ": no columns of that category in the column file; skipped."); next }
  if (is.null(panel_rows[[pn]]) || !any(panel_rows[[pn]]$shown)) { message(pn, ": no modules to show; skipped."); next }
  n <- draw_version(pn, panel_cols[[pn]], panel_rows[[pn]]); drawn <- c(drawn, paste0(pn, "=", n))
}

write_lines_safe(c(paste0(names(opt), "\t", vapply(opt, function(x) if (is.null(x)) "" else paste(x, collapse = ","), "")),
                   paste0("kme_files_used\t", if (is.null(kme_use)) "" else kme_use),
                   paste0("versions_drawn\t", paste(drawn, collapse = "; "))), file.path(out_dir, "run_parameters.txt"))
message("Script 12f complete: ", paste(drawn, collapse = ", "), " -> ", out_dir)