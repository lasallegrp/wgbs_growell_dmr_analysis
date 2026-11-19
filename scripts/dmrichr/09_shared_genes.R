#!/usr/bin/env Rscript
#=====================================================================================
#  09_shared_genes.R
#  Gene lists for enrichment of genes SHARED between DMRichR comparisons, made from
#  the gene-mode output of 05_dmrichr_overlap.R (<group>_overlap_gene.xlsx, sheet Genes).
#
#  Foreground: genes with a DMR in every comparison of the group (--shared all) or in
#              at least two (--shared any2); optionally direction-concordant only, and
#              optionally split into shared-hyper and shared-hypo lists.
#  Background: genes TESTABLE in every DMRichR set of the group = intersection of the
#              geneSymbol columns of each set's DMRs/background_annotated.xlsx (sets with no
#              DMRichR background, e.g. a GWAS gene list, are left out: --no_background_sets).
#              Same DMRichR annotation as the DMR genes, so the test is fair.
#
#  Writes to --out_dir:
#    <list>_genes.txt          one symbol per line, for 07_enrichment_plot.R --genes
#    background_genes.txt      for --background_genes
#    lists.tsv                 list, n_genes, run (TRUE if >= --min_genes)
#
#  Called by the job that 09_make_shared_jobs.py writes; can also be run alone:
#    Rscript scripts/dmrichr/09_shared_genes.R \
#      --overlap_xlsx results/dmrichr/overlap/DMR_AZ_04/custom/custom_AZ04_career_by_timepoint/gene/custom_AZ04_career_by_timepoint_overlap_gene.xlsx \
#      --manifest results/dmrichr/overlap/DMR_AZ_04/manifest.tsv \
#      --out_dir  .../gene/shared_gene_enrichment
#=====================================================================================
suppressPackageStartupMessages({ library(optparse); library(readxl) })

opt <- parse_args(OptionParser(option_list = list(
  make_option("--overlap_xlsx", type = "character", help = "<group>_overlap_gene.xlsx from 05_dmrichr_overlap.R"),
  make_option("--manifest", type = "character", default = NULL, help = "manifest.tsv (set_id, path[, method]); optional with --dmrichr_root"),
  make_option("--dmrichr_root", type = "character", default = NULL,
              help = "Without --manifest: folder of <outcome>/<timepoint>/ (results/dmrichr_output); sets are matched by name"),
  make_option("--out_dir", type = "character"),
  make_option("--background_file", type = "character", default = "DMRs/background_annotated.xlsx"),
  make_option("--shared", type = "character", default = "all",
              help = "all = in every comparison of the group; any2 = in at least two [default %default]"),
  make_option("--concordant", action = "store_true", default = FALSE,
              help = "Keep only genes changed in the same direction in every comparison they are in"),
  make_option("--split_direction", action = "store_true", default = FALSE,
              help = "Also write shared_hyper and shared_hypo lists (direction-concordant genes only)"),
  make_option("--no_background_sets", type = "character", default = "",
              help = "Comma list of set_ids with no DMRichR background (e.g. a GWAS gene list); left out of the background intersection"),
  make_option("--min_genes", type = "integer", default = 20L,
              help = "Lists with fewer genes are written but marked run = FALSE [default %default]")
)))
for (a in c("overlap_xlsx", "out_dir")) if (is.null(opt[[a]])) stop("--", a, " is required")
stopifnot(opt$shared %in% c("all", "any2"))
dir.create(opt$out_dir, recursive = TRUE, showWarnings = FALSE)

# ---- shared genes -------------------------------------------------------------------
gt <- as.data.frame(read_excel(opt$overlap_xlsx, sheet = "Genes"))
dcols <- grep("^dir_", names(gt), value = TRUE)
sets <- sub("^dir_", "", dcols)
if (length(sets) < 2) stop("Expected >= 2 dir_<set> columns in the Genes sheet; found: ", paste(dcols, collapse = ", "))
need <- if (opt$shared == "all") length(sets) else 2L
gt <- gt[!is.na(gt$geneSymbol) & gt$geneSymbol != "" & gt$n_sets >= need, ]
conc <- gt$direction_concordant %in% c(TRUE, "TRUE")
if (opt$concordant) gt <- gt[conc, ]; conc <- gt$direction_concordant %in% c(TRUE, "TRUE")
# direction of a concordant gene = the (single) direction label across its sets
gene_dir <- apply(gt[, dcols, drop = FALSE], 1, function(v) { v <- unique(na.omit(v)); if (length(v) == 1) v else "mixed" })

lists <- list(shared = gt$geneSymbol)
if (opt$split_direction) {
  lists$shared_hyper <- gt$geneSymbol[conc & gene_dir == "hyper"]
  lists$shared_hypo  <- gt$geneSymbol[conc & gene_dir == "hypo"]
}
lists <- lapply(lists, function(x) sort(unique(x)))

# ---- background: genes testable in every comparison ----------------------------------
# ---- set_id -> DMRichR folder: from --manifest, or matched by name under --dmrichr_root -----
# Name matching (case-insensitive; "-" and "_" are the same): a set is <outcome>_<timepoint>
# (or <timepoint>_<outcome>), or just the timepoint when the group is an outcome
# (within_outcome_across_timepoints/egwg), or just the outcome when the group is a timepoint
# (within_timepoint_across_outcomes/Baseline). Sets matching no folder (e.g. a GWAS list) have
# no DMRichR background.
norm <- function(x) gsub("[^a-z0-9]+", "_", tolower(x))
auto_manifest <- function(sets, root, group) {
  si <- Sys.glob(file.path(root, "*", "*", "sample_info.xlsx"))
  if (!length(si)) stop("No <outcome>/<timepoint>/sample_info.xlsx under --dmrichr_root ", root)
  d <- dirname(si); o <- basename(dirname(d)); t <- basename(d); g <- norm(group)
  path <- vapply(sets, function(s) {
    k <- norm(s)
    hit <- which(k == norm(paste(o, t)) | k == norm(paste(t, o)) |
                 (k == norm(t) & norm(o) == g) | (k == norm(o) & norm(t) == g))
    if (length(hit) > 1) stop("Set ", s, " matches several folders: ", paste(d[hit], collapse = ", "))
    if (length(hit)) d[hit] else NA_character_
  }, "")
  data.frame(set_id = sets, method = ifelse(is.na(path), "external", "dmrichr"), path = path,
             stringsAsFactors = FALSE)
}
if (!is.null(opt$manifest)) {
  man <- read.delim(opt$manifest, stringsAsFactors = FALSE)
  if (!all(c("set_id", "path") %in% names(man))) stop(opt$manifest, " needs set_id and path columns")
} else {
  if (is.null(opt$dmrichr_root)) stop("Give --manifest or --dmrichr_root")
  man <- auto_manifest(sets, opt$dmrichr_root, basename(dirname(dirname(normalizePath(opt$overlap_xlsx)))))
  for (i in seq_len(nrow(man))) message(sprintf("  set %-28s -> %s", man$set_id[i],
                                                if (is.na(man$path[i])) "no DMRichR folder (no background)" else man$path[i]))
  write.table(man, file.path(opt$out_dir, "resolved_sets.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
}
miss <- setdiff(sets, man$set_id)
if (length(miss)) stop("Sets not in ", opt$manifest, ": ", paste(miss, collapse = ", "))
# a set's folder = the DMRichR timepoint folder (e.g. dmrichr_output/ppwr/36-38wk), whether the
# manifest gives that folder, its DMRs/ folder, or a file in either; relative paths are tried
# as given and then relative to the manifest's folder
set_dir <- function(s) {
  p <- man$path[man$set_id == s][1]
  if (is.na(p)) return(NA_character_)
  if (!grepl("^/", p) && !file.exists(p) && !is.null(opt$manifest)) {
    p2 <- file.path(dirname(normalizePath(opt$manifest)), p); if (file.exists(p2)) p <- p2
  }
  if (!dir.exists(p) && (file.exists(p) || grepl("[.][A-Za-z0-9]+$", basename(p)))) p <- dirname(p)   # a file -> its folder
  if (basename(p) == "DMRs") p <- dirname(p)
  p
}
# sets with no DMRichR background (external lists, e.g. GWAS genes) stay out of the background
nobg <- trimws(strsplit(opt$no_background_sets, ",")[[1]]); nobg <- nobg[nzchar(nobg)]
if ("method" %in% names(man)) nobg <- union(nobg, man$set_id[!is.na(man$method) & tolower(man$method) != "dmrichr"])
bg_ids <- setdiff(sets, nobg)
if (!length(bg_ids)) stop("No set in this group has a DMRichR background (all are in --no_background_sets / method != dmrichr)")
if (length(intersect(sets, nobg))) message("  Not used for the background (no DMRichR background): ",
                                           paste(intersect(sets, nobg), collapse = ", "))
bg_sets <- lapply(bg_ids, function(s) {
  f <- file.path(set_dir(s), opt$background_file)
  if (!file.exists(f)) stop("Missing background for ", s, ": ", f)
  g <- as.character(read_excel(f)[["geneSymbol"]])
  message(sprintf("  %s: %s background genes", s, format(length(unique(na.omit(g))), big.mark = ",")))
  unique(na.omit(g[g != ""]))
})
bg <- Reduce(intersect, bg_sets)
not_in_bg <- setdiff(unique(unlist(lists)), bg)   # e.g. symbols renamed by the overlap step
if (length(not_in_bg)) {
  message("  ", length(not_in_bg), " shared genes not in the common background were added to it: ",
          paste(head(not_in_bg, 10), collapse = ", "), if (length(not_in_bg) > 10) ", ...")
  bg <- union(bg, not_in_bg)
}
bg <- sort(bg)
writeLines(bg, file.path(opt$out_dir, "background_genes.txt"))

# ---- write lists ------------------------------------------------------------------
summ <- data.frame(list = names(lists), n_genes = lengths(lists),
                   run = lengths(lists) >= opt$min_genes, stringsAsFactors = FALSE)
for (n in names(lists)) writeLines(lists[[n]], file.path(opt$out_dir, paste0(n, "_genes.txt")))
write.table(summ, file.path(opt$out_dir, "lists.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
writeLines(c(paste("overlap_xlsx:", opt$overlap_xlsx), paste("sets:", paste(sets, collapse = ", ")),
             paste("shared:", opt$shared), paste("concordant:", opt$concordant),
             paste("background sets:", paste(bg_ids, collapse = ", ")),
             paste("background genes (testable in all background sets):", length(bg))),
           file.path(opt$out_dir, "selection.txt"))
message(sprintf("Sets: %s\nBackground (testable in all): %s genes", paste(sets, collapse = " & "),
                format(length(bg), big.mark = ",")))
print(summ, row.names = FALSE)