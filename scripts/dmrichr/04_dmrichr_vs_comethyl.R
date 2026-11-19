#!/usr/bin/env Rscript
#=====================================================================================
#  04_dmrichr_vs_comethyl.R
#  Do DMRichR DMRs converge on comethyl CONSENSUS modules that associate with the same
#  outcome? (GROWELL outcome x timepoint layout)
#
#  DMR sets (discovered automatically, as in 03):
#    <results_dir>/<outcome>/<timepoint>/DMRs/DMRs_annotated.xlsx (+ background_annotated.xlsx)
#  Consensus network (<comethyl_root> = folder holding the numbered comethyl steps):
#    regions + module + genes (shared)   15_module_annotation/<cpg>/<filter>/<variant>/Annotated_Regions.tsv
#    kME, per timepoint                  11_consensus_membership/<tp>/<cpg>/<filter>/<variant>/<tp>_Consensus_Region_Module_Membership.txt
#    ME-trait bicor, per timepoint       12a_me_trait_analysis/<tp>/<cpg>/<filter>/<variant>/ME_Trait_Correlation_Stats_Bicor.tsv
#  Each DMR set is paired with the consensus modules using ITS timepoint's kME and
#  ME-trait table; the outcome is mapped to the comethyl trait with --trait_map.
#
#  Per pair:
#   1. Map DMRs onto network regions (coordinates), within DMRichR-tested territory
#   2. Module enrichment: Fisher per module, region level (and gene level)
#   3. Link to ME-trait bicor for the outcome's trait
#   4. Direction concordance: expected DMR direction = sign(bicor) * sign(kME)
#   5. kME: are DMR-overlapping regions more central within a module?
#   6. All-trait profile of the DMR-enriched modules (which other traits they track)
#  Summaries: per outcome across timepoints (persistence) and per timepoint across outcomes.
#
#  Example:
#   Rscript 04_dmrichr_vs_comethyl.R --project_dir /quobyte/.../2025_DBS_dmrich_George \
#       --comethyl_root /quobyte/.../results/comethyl/consensus
#=====================================================================================
suppressPackageStartupMessages({ library(optparse) })

opt <- parse_args(OptionParser(option_list = list(
  make_option("--project_dir", type = "character", help = "Project root [required]"),
  make_option("--helper", type = "character", default = NA, help = "Path to helper_overlap.R"),
  # ---- DMRichR sets ----
  make_option("--results_dir", type = "character", default = NA,
              help = "Folder holding <outcome>/<timepoint>/ [default <project_dir>/results/dmrichr_output]"),
  make_option("--outcomes", type = "character", default = NA, help = "Comma list [default: all discovered]"),
  make_option("--timepoints", type = "character", default = "Baseline,36-38wk,Postpartum",
              help = "Comma list, in display order; must match comethyl timepoint folders [default %default]"),
  make_option("--subdir", type = "character", default = "", help = "Extra level above DMRs/ [default none]"),
  make_option("--exclude_dirs", type = "character", default = "logs,metadata,DMRichR_overlap,DMRichR_vs_comethyl",
              help = "Folders in results_dir that are not outcomes [default %default]"),
  make_option("--dmr_file", type = "character", default = "DMRs/DMRs_annotated.xlsx"),
  make_option("--background_file", type = "character", default = "DMRs/background_annotated.xlsx",
              help = "Restrict network regions to DMRichR-tested territory. 'none' to skip"),
  make_option("--effect_col", type = "character", default = "auto"),
  make_option("--gene_col", type = "character", default = "geneSymbol"),
  make_option("--gene_sep", type = "character", default = NA),
  # ---- comethyl consensus network ----
  make_option("--comethyl_root", type = "character", default = NA,
              help = "Folder with 11_consensus_membership/, 12a_me_trait_analysis/, 15_module_annotation/ [default: auto under <project_dir>/results]"),
  make_option("--cpg_filter", type = "character", default = "cov3_75pct", help = "[default %default]"),
  make_option("--region_filter", type = "character", default = "covMin14_methSD0p06", help = "[default %default]"),
  make_option("--variant", type = "character", default = "v2_exclude_protected_pcs", help = "[default %default]"),
  make_option("--regions_file", type = "character", default = "Annotated_Regions.tsv"),
  make_option("--membership_file", type = "character", default = "{tp}_Consensus_Region_Module_Membership.txt",
              help = "{tp} is replaced by the timepoint [default %default]"),
  make_option("--me_trait_file", type = "character", default = "ME_Trait_Correlation_Stats_Bicor.tsv"),
  make_option("--trait_map", type = "character", default = "apo=apo,egwg=egwg,ppwr=ppwr_e,experimental_group=group",
              help = "outcome=comethyl_trait pairs [default %default]"),
  make_option("--exclude_traits", type = "character", default = "Percent_Duplicated_Reads,Percent_Aligned_Total,Total_Reads",
              help = "Traits left out of the all-trait profile (ME* traits are always left out)"),
  # ---- analysis ----
  make_option("--out_dir", type = "character", default = NA,
              help = "[default <results_dir>/DMRichR_vs_comethyl/<region_filter>/<variant><out_suffix>]"),
  make_option("--out_suffix", type = "character", default = ""),
  make_option("--overlap_mode", type = "character", default = "both", help = "coord | gene | both"),
  make_option("--max_gap", type = "integer", default = 0L),
  make_option("--min_overlap_bp", type = "integer", default = 1L),
  make_option("--exclude_sex_chr", action = "store_true", default = FALSE),
  make_option("--module_trait_p", type = "double", default = 0.05,
              help = "p threshold calling a module trait-associated [default %default]"),
  make_option("--fdr", type = "double", default = 0.05),
  make_option("--grey_module", type = "character", default = "grey"),
  make_option("--min_dmrs_plot", type = "integer", default = 3L,
              help = "Modules with fewer DMRs are faded in plots and left out of the Spearman test [default %default]"),
  make_option("--module_genes", type = "integer", default = 10L,
              help = "Gene names labeled per module in *_module_DMR_genes.pdf (0 = no plot) [default %default]"),
  make_option("--min_dmrs_concordance", type = "integer", default = 5L,
              help = "Min DMRs in a module to report a direction test"),
  make_option("--profile_modules", type = "integer", default = 6L,
              help = "If no module is DMR-enriched, all-trait profile of this many top modules (by enrichment p) [default %default]")
)))
if (is.null(opt$project_dir)) stop("--project_dir is required")
.a <- grep("^--file=", commandArgs(FALSE), value = TRUE)
self <- if (length(.a)) normalizePath(sub("^--file=", "", .a[1])) else NA_character_
if (is.na(opt$helper)) opt$helper <- file.path(dirname(self), "helper_overlap.R")
source(opt$helper)
suppressPackageStartupMessages(library(openxlsx))
stopifnot(opt$overlap_mode %in% c("coord", "gene", "both"))
do_coord <- opt$overlap_mode %in% c("coord", "both")
do_gene  <- opt$overlap_mode %in% c("gene", "both")
has_repel <- requireNamespace("ggrepel", quietly = TRUE)
split_opt <- function(x) if (is.na(x) || !nzchar(x)) NULL else trimws(strsplit(x, ",")[[1]])
clean_id  <- function(x) gsub("^_|_$", "", gsub("[^A-Za-z0-9]+", "_", x))
read_tsv_fast <- function(p) {
  if (requireNamespace("data.table", quietly = TRUE))
    as.data.frame(data.table::fread(p, sep = "\t", quote = "", data.table = FALSE, check.names = FALSE))
  else read.delim(p, quote = "", check.names = FALSE, stringsAsFactors = FALSE)
}
# genes in a cell may be "A, B" or "A;B"
split_gene_cells <- function(x) {
  x <- ifelse(is.na(x), "", as.character(x))
  lapply(strsplit(x, "[,;|] *"), function(v) { v <- trimws(v); v[v != "" & v != "NA"] })
}

if (is.na(opt$results_dir)) opt$results_dir <- file.path(opt$project_dir, "results", "dmrichr_output")
timepoints <- split_opt(opt$timepoints)

# ----------------------------- comethyl paths ----------------------------------------
regions_rel <- file.path("15_module_annotation", opt$cpg_filter, opt$region_filter, opt$variant, opt$regions_file)
if (is.na(opt$comethyl_root)) {
  cand <- file.path(opt$project_dir, "results",
                    c("comethyl/consensus", "comethyl_output/consensus", "comethyl/comethyl_output/consensus",
                      "comethyl", "comethyl_output"))
  hit <- cand[file.exists(file.path(cand, regions_rel))]
  if (!length(hit)) stop("Could not find ", regions_rel, " under any of:\n  ", paste(cand, collapse = "\n  "),
                         "\nGive --comethyl_root (the folder holding 15_module_annotation/)")
  opt$comethyl_root <- hit[1]
}
cpath <- function(step, tp, file) file.path(opt$comethyl_root, step, tp, opt$cpg_filter, opt$region_filter,
                                             opt$variant, gsub("{tp}", tp, file, fixed = TRUE))
regions_path <- file.path(opt$comethyl_root, regions_rel)
if (!file.exists(regions_path)) stop("Annotated regions not found: ", regions_path)

# ----------------------------- DMR sets -----------------------------------------------
discover_sets <- function(results_dir, dmr_file, subdir, outcomes, timepoints, exclude) {
  if (!dir.exists(results_dir)) stop("results_dir not found: ", results_dir)
  found <- setdiff(list.dirs(results_dir, recursive = FALSE, full.names = FALSE), exclude)
  if (!is.null(outcomes)) {
    miss <- setdiff(outcomes, found)
    if (length(miss)) message("Outcome folders not found, skipped: ", paste(miss, collapse = ", "))
    found <- outcomes[outcomes %in% found]
  }
  rows <- list()
  for (o in found) {
    tp <- list.dirs(file.path(results_dir, o), recursive = FALSE, full.names = FALSE)
    extra <- setdiff(tp, timepoints)
    if (length(extra)) message(o, ": folders not in --timepoints, skipped: ", paste(extra, collapse = ", "))
    for (t in timepoints[timepoints %in% tp]) {
      p <- if (nzchar(subdir)) file.path(results_dir, o, t, subdir) else file.path(results_dir, o, t)
      if (file.exists(file.path(p, dmr_file))) {
        rows[[length(rows) + 1]] <- data.frame(set_id = clean_id(paste(o, t, sep = "_")), outcome = o,
                                               timepoint = t, path = p, stringsAsFactors = FALSE)
      } else message("No ", dmr_file, " in ", p, " -> skipped")
    }
  }
  if (!length(rows)) stop("No DMR sets found under ", results_dir)
  do.call(rbind, rows)
}
dm <- discover_sets(opt$results_dir, opt$dmr_file, opt$subdir, split_opt(opt$outcomes), timepoints,
                    split_opt(opt$exclude_dirs))
tmap <- split_opt(opt$trait_map)
tmap <- setNames(sub("^[^=]*=", "", tmap), sub("=.*$", "", tmap))
dm$trait <- unname(tmap[dm$outcome])
if (any(is.na(dm$trait))) message("No --trait_map entry for: ", paste(unique(dm$outcome[is.na(dm$trait)]), collapse = ", "),
                                  " -> trait-linked statistics will be NA for those sets")
dm$dmr_path <- file.path(dm$path, opt$dmr_file)
dm$bg_path  <- if (tolower(opt$background_file) == "none") NA else file.path(dm$path, opt$background_file)

if (is.na(opt$out_dir))
  opt$out_dir <- file.path(opt$results_dir, "DMRichR_vs_comethyl", opt$region_filter, paste0(opt$variant, opt$out_suffix))
dir.create(opt$out_dir, recursive = TRUE, showWarnings = FALSE)
write.table(dm, file.path(opt$out_dir, "sets_used.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# ----------------------------- load consensus network -----------------------------
message("Comethyl root: ", opt$comethyl_root)
message("Reading consensus regions: ", regions_path)
reg0 <- read_tsv_fast(regions_path)
miss <- setdiff(c("RegionID", "chr", "start", "end", "module"), names(reg0))
if (length(miss)) stop("Annotated regions missing: ", paste(miss, collapse = ", "))
reg0$start <- as.integer(reg0$start); reg0$end <- as.integer(reg0$end)
if (!"gene_symbol" %in% names(reg0)) reg0$gene_symbol <- NA
if (!"hubRegion" %in% names(reg0)) reg0$hubRegion <- FALSE
reg0$hubRegion <- as.logical(reg0$hubRegion)
if (opt$exclude_sex_chr) reg0 <- reg0[!reg0$chr %in% c("chrX", "chrY"), ]
reg0_genes <- if (do_gene) split_gene_cells(reg0$gene_symbol) else NULL
message("  ", nrow(reg0), " regions, ", length(unique(reg0$module)), " modules (incl. ", opt$grey_module, ")")

tps_needed <- unique(dm$timepoint)
ex_traits <- split_opt(opt$exclude_traits)
nets <- list(); input_files <- regions_path
for (tp in tps_needed) {
  mp <- cpath("11_consensus_membership", tp, opt$membership_file)
  tpf <- cpath("12a_me_trait_analysis", tp, opt$me_trait_file)
  reg <- reg0
  if (file.exists(mp)) {
    message("Reading ", tp, " membership: ", mp)
    mem <- read_tsv_fast(mp)
    mcol <- detect_col(mem, c("Module", "module"), "module")
    m <- as.matrix(mem[, setdiff(names(mem), c("RegionID", mcol)), drop = FALSE])
    j <- match(mem[[mcol]], colnames(m))
    own <- rep(NA_real_, nrow(m)); ok <- !is.na(j)
    own[ok] <- m[cbind(which(ok), j[ok])]
    k <- match(reg$RegionID, mem$RegionID)
    reg$membership <- own[k]
    agree <- mean(reg$module[!is.na(k)] == mem[[mcol]][k[!is.na(k)]])
    message(sprintf("  kME for %d of %d regions; module labels agree for %.1f%%",
                    sum(!is.na(reg$membership)), nrow(reg), 100 * agree))
    if (!is.na(agree) && agree < 0.99) warning(tp, ": module labels differ between Annotated_Regions and membership file", call. = FALSE)
    rm(mem, m); input_files <- c(input_files, mp)
  } else {
    message(tp, ": membership file not found (", mp, ") -> kME tests skipped, expected direction = sign(bicor)")
    reg$membership <- NA_real_
  }
  if (!file.exists(tpf)) stop(tp, ": ME-trait table not found: ", tpf)
  mt <- read_tsv_fast(tpf)
  miss <- setdiff(c("module", "trait", "bicor", "p"), names(mt))
  if (length(miss)) stop("ME-trait table missing columns: ", paste(miss, collapse = ", "), " in ", tpf)
  mt$bicor <- as.numeric(mt$bicor); mt$p <- as.numeric(mt$p)
  input_files <- c(input_files, tpf)
  nets[[tp]] <- list(reg = reg, gr = to_gr(reg, "chr", "start", "end"), mt = mt,
                     mt_all = mt[!startsWith(mt$trait, "ME") & !mt$trait %in% ex_traits, ])
  miss_t <- setdiff(na.omit(unique(dm$trait[dm$timepoint == tp])), mt$trait)
  if (length(miss_t)) warning(tp, ": trait(s) not in ME-trait table: ", paste(miss_t, collapse = ", "), call. = FALSE)
}
log_run_metadata(opt$out_dir, opt, c(dm$dmr_path, na.omit(dm$bg_path), input_files), self)
message("\n", nrow(dm), " DMR sets x consensus network (matched timepoint)")

col_ok <- function(x) vapply(x, function(z) !inherits(try(col2rgb(z), silent = TRUE), "try-error"), TRUE)
all_mod <- list(); all_sum <- list(); all_prof <- list()

# ----------------------------- per DMR set ----------------------------------------
for (k in seq_len(nrow(dm))) {
  drow <- dm[k, ]; sid <- drow$set_id; tp <- drow$timepoint; net <- nets[[tp]]
  nid <- paste0("consensus_", clean_id(tp))
  tag <- paste0(sid, "__vs__", nid)
  trait_lab <- ifelse(is.na(drow$trait), "(no trait)", drow$trait)
  odir <- file.path(opt$out_dir, "per_set", drow$outcome, tag); dir.create(odir, recursive = TRUE, showWarnings = FALSE)
  message("\n== ", tag, "  (trait: ", trait_lab, ")")
  d <- read_dmr_set(drow$dmr_path, sid, opt$effect_col, opt$gene_col, opt$gene_sep, opt$exclude_sex_chr)
  dgr <- to_gr(d)
  reg <- net$reg; rgr <- net$gr; rgenes <- reg0_genes

  # restrict network to territory DMRichR could test (else non-overlap is uninformative)
  bg <- read_background(drow$bg_path, opt$exclude_sex_chr)
  in_u <- if (is.null(bg)) rep(TRUE, nrow(reg)) else overlapsAny(rgr, bg)
  reg <- reg[in_u, ]; rgr <- rgr[in_u]; if (do_gene) rgenes <- rgenes[in_u]
  message("  network regions in DMRichR territory: ", nrow(reg), " of ", nrow(net$reg))

  # module-trait for this outcome's trait
  mt <- net$mt[!is.na(drow$trait) & net$mt$trait == drow$trait, c("module", "bicor", "p")]
  names(mt) <- c("module", "ME_trait_bicor", "ME_trait_p")

  # ---- 1. mapping ----
  h <- ov_hits(dgr, rgr, opt$max_gap, opt$min_overlap_bp)
  hits <- data.frame(dmr_id = d$.dmr_id[queryHits(h)], dmr_effect = d$.effect[queryHits(h)],
                     dmr_gene = d$.gene[queryHits(h)], reg[subjectHits(h),
                     intersect(c("RegionID", "chr", "start", "end", "module", "membership", "hubRegion",
                                 "gene_symbol", "genomic_location"), names(reg)), drop = FALSE])
  hits <- merge(hits, mt, by = "module", all.x = TRUE)
  hits$expected_dir <- sign(hits$ME_trait_bicor) * ifelse(is.na(hits$membership), 1, sign(hits$membership))
  hits$observed_dir <- sign(hits$dmr_effect)
  reg$dmr_hit <- seq_len(nrow(reg)) %in% subjectHits(h)

  dmr_tab <- if (nrow(hits)) do.call(rbind, lapply(split(hits, hits$dmr_id), function(x) {
    ng <- x[x$module != opt$grey_module, ]
    tabm <- sort(table(ng$module), decreasing = TRUE)
    data.frame(dmr_id = x$dmr_id[1], n_network_regions = nrow(x),
               modules_hit = paste(sort(unique(x$module)), collapse = ";"),
               primary_module = if (length(tabm)) names(tabm)[1] else opt$grey_module,
               max_membership = suppressWarnings(round(max(x$membership, na.rm = TRUE), 3)),
               overlaps_hub = any(x$hubRegion %in% TRUE))
  })) else data.frame(dmr_id = character(0), primary_module = character(0))
  dmr_out <- cbind(dmr_id = d$.dmr_id, d[, !startsWith(names(d), "."), drop = FALSE])
  dmr_out <- merge(dmr_out, dmr_tab, by = "dmr_id", all.x = TRUE)
  dmr_out$primary_module[is.na(dmr_out$primary_module)] <- "(no network region)"

  # ---- 2-5. per-module table ----
  mods <- sort(unique(reg$module))
  modtab <- do.call(rbind, lapply(mods, function(m) {
    inm <- reg$module == m
    a <- sum(inm & reg$dmr_hit); b <- sum(inm & !reg$dmr_hit)
    c_ <- sum(!inm & reg$dmr_hit); d_ <- sum(!inm & !reg$dmr_hit)
    ft <- fisher.test(matrix(c(a, b, c_, d_), 2))
    hx <- hits[hits$module == m, ]
    dm1 <- hx[!duplicated(hx$dmr_id), ]
    n_hyper <- sum(dm1$observed_dir > 0); n_hypo <- sum(dm1$observed_dir < 0)
    conc <- dm1$observed_dir == dm1$expected_dir
    nc <- sum(conc, na.rm = TRUE); nt <- sum(!is.na(conc))
    bp <- if (nt >= opt$min_dmrs_concordance) binom.test(nc, nt, 0.5)$p.value else NA
    km1 <- na.omit(reg$membership[inm & reg$dmr_hit]); km0 <- na.omit(reg$membership[inm & !reg$dmr_hit])
    wp <- if (m != opt$grey_module && length(km1) >= 3 && length(km0) >= 3)
      suppressWarnings(wilcox.test(km1, km0)$p.value) else NA
    data.frame(module = m, n_regions = sum(inm), n_regions_with_DMR = a,
               pct_regions_with_DMR = round(100 * a / sum(inm), 3),
               n_DMRs = nrow(dm1), n_DMR_hyper = n_hyper, n_DMR_hypo = n_hypo,
               region_OR = unname(ft$estimate), region_p = ft$p.value,
               region_log2OR_adj = log2(((a + .5) * (d_ + .5)) / ((b + .5) * (c_ + .5))),
               n_dir_tested = nt, pct_dir_concordant = if (nt) round(100 * nc / nt, 1) else NA,
               dir_binom_p = bp,
               median_kME_DMR_regions = if (length(km1)) round(median(km1), 3) else NA,
               median_kME_other_regions = if (length(km0)) round(median(km0), 3) else NA, kME_wilcox_p = wp,
               hub_overlapped = any(reg$hubRegion[inm] & reg$dmr_hit[inm]))
  }))
  ng <- modtab$module != opt$grey_module
  modtab$region_FDR <- NA; modtab$region_FDR[ng] <- p.adjust(modtab$region_p[ng], "BH")

  if (do_gene) {
    ug <- unique(unlist(rgenes)); dg <- intersect(dmr_genes(d), ug)
    gstat <- do.call(rbind, lapply(mods, function(m) {
      mg <- unique(unlist(rgenes[reg$module == m]))
      tb <- table(factor(ug %in% mg, c(TRUE, FALSE)), factor(ug %in% dg, c(TRUE, FALSE)))
      ft <- fisher.test(tb)
      data.frame(module = m, n_module_genes = length(mg), n_shared_genes = sum(mg %in% dg),
                 gene_OR = unname(ft$estimate), gene_p = ft$p.value,
                 gene_log2OR_adj = log2(((tb[1, 1] + .5) * (tb[2, 2] + .5)) / ((tb[1, 2] + .5) * (tb[2, 1] + .5))),
                 shared_genes = paste(head(sort(intersect(mg, dg)), 50), collapse = ";"))
    }))
    gstat$gene_FDR <- NA; gstat$gene_FDR[gstat$module != opt$grey_module] <-
      p.adjust(gstat$gene_p[gstat$module != opt$grey_module], "BH")
    modtab <- merge(modtab, gstat, by = "module", all.x = TRUE)
  }
  modtab <- merge(modtab, mt, by = "module", all.x = TRUE)
  modtab$trait_associated <- modtab$ME_trait_p < opt$module_trait_p
  lvl <- if (do_coord) "region" else "gene"
  modtab$enrich_level  <- lvl
  modtab$enrich_log2OR <- modtab[[paste0(lvl, "_log2OR_adj")]]
  modtab$enrich_FDR    <- modtab[[paste0(lvl, "_FDR")]]
  modtab$DMR_enriched  <- modtab$enrich_FDR < opt$fdr & modtab[[paste0(lvl, "_OR")]] > 1
  modtab <- modtab[order(modtab[[paste0(lvl, "_p")]]), ]

  # convergence across modules: does enrichment track trait association?
  mm <- modtab[modtab$module != opt$grey_module & !is.na(modtab$ME_trait_bicor), ]
  mr <- mm[mm$n_DMRs >= opt$min_dmrs_plot, ]
  rho <- if (nrow(mr) >= 5) suppressWarnings(cor.test(mr$enrich_log2OR, abs(mr$ME_trait_bicor),
                                                      method = "spearman")) else NULL
  cross <- table(enriched = factor(mm$DMR_enriched, c(TRUE, FALSE)),
                 trait_assoc = factor(mm$trait_associated, c(TRUE, FALSE)))
  cross_p <- if (all(dim(cross) == 2) && sum(cross)) fisher.test(cross)$p.value else NA

  cov <- data.frame(
    dmr_set = sid, network = nid, outcome = drow$outcome, timepoint = tp, trait = trait_lab,
    region_filter = opt$region_filter, variant = opt$variant,
    network_regions_total = nrow(net$reg), network_regions_in_DMRichR_territory = nrow(reg),
    n_DMRs = nrow(d), n_DMRs_hitting_network = length(unique(queryHits(h))),
    pct_DMRs_hitting_network = round(100 * length(unique(queryHits(h))) / max(1, nrow(d)), 1),
    n_DMRs_hitting_module = sum(!dmr_out$primary_module %in% c(opt$grey_module, "(no network region)")),
    n_modules_enriched = sum(modtab$DMR_enriched, na.rm = TRUE),
    n_modules_trait_associated = sum(mm$trait_associated, na.rm = TRUE),
    n_enriched_and_trait_associated = sum(mm$DMR_enriched & mm$trait_associated, na.rm = TRUE),
    enriched_vs_trait_fisher_p = cross_p,
    spearman_rho_logOR_vs_absBicor = if (is.null(rho)) NA else unname(rho$estimate),
    spearman_p = if (is.null(rho)) NA else rho$p.value,
    spearman_n_modules = nrow(mr),
    pct_dir_concordant_trait_modules = with(modtab[modtab$trait_associated %in% TRUE, ],
      if (sum(n_dir_tested, na.rm = TRUE)) round(100 * sum(pct_dir_concordant * n_dir_tested / 100, na.rm = TRUE) /
                                                   sum(n_dir_tested, na.rm = TRUE), 1) else NA),
    stringsAsFactors = FALSE)
  print(t(cov))

  # ---- 6. all-trait profile of DMR-enriched modules (or top modules if none) ----
  prof_mods <- modtab$module[modtab$DMR_enriched %in% TRUE & modtab$module != opt$grey_module]
  prof_note <- "DMR-enriched modules (FDR < cutoff)"
  if (!length(prof_mods)) {
    cand <- modtab[modtab$module != opt$grey_module & modtab$n_DMRs >= opt$min_dmrs_plot, ]
    prof_mods <- head(cand$module, opt$profile_modules)
    prof_note <- paste0("no module DMR-enriched; top ", length(prof_mods), " by enrichment p (>= ",
                        opt$min_dmrs_plot, " DMRs)")
  }
  prof <- net$mt_all[net$mt_all$module %in% prof_mods, c("module", "trait", "nObs", "bicor", "p")]
  if (nrow(prof)) {
    prof <- merge(prof, modtab[, c("module", "n_DMRs", "enrich_log2OR", "enrich_FDR", "DMR_enriched")], by = "module")
    all_prof[[tag]] <- cbind(dmr_set = sid, outcome = drow$outcome, timepoint = tp, selection = prof_note, prof)
    pr <- prof
    pr$mlab <- paste0(pr$module, " (", pr$n_DMRs, " DMRs", ifelse(pr$DMR_enriched %in% TRUE, ", enriched", ""), ")")
    pr$star <- ifelse(pr$p < 0.001, "***", ifelse(pr$p < 0.01, "**", ifelse(pr$p < 0.05, "*", "")))
    pr$is_outcome <- !is.na(drow$trait) & pr$trait == drow$trait
    tr_order <- unique(c(if (!is.na(drow$trait)) drow$trait, sort(unique(pr$trait))))
    pr$trait <- factor(pr$trait, tr_order[tr_order %in% pr$trait])
    pp <- ggplot(pr, aes(trait, mlab, fill = bicor)) + geom_tile(colour = "white") +
      geom_tile(data = pr[pr$is_outcome, ], colour = "black", linewidth = .8, fill = NA) +
      geom_text(aes(label = star), size = 3.2, vjust = .75) +
      scale_fill_gradient2(low = "#2C7BB6", mid = "white", high = "#D7191C", limits = c(-1, 1) * max(.3, max(abs(pr$bicor), na.rm = TRUE)),
                           name = "ME-trait bicor") +
      labs(x = NULL, y = NULL, title = paste(sid, "- what else do these modules track?"),
           subtitle = paste0(prof_note, "; ", tp, " ME-trait table. Box = this set's outcome (", trait_lab,
                             "); * p<0.05 ** p<0.01 *** p<0.001 (unadjusted)")) +
      theme_minimal(10) + theme(axis.text.x = element_text(angle = 60, hjust = 1),
                                plot.subtitle = element_text(size = 8.5))
    ggsave(file.path(odir, paste0(tag, "_enriched_modules_all_traits.pdf")), pp,
           width = max(8, .28 * length(levels(pr$trait)) + 4), height = 2.5 + .4 * length(unique(pr$mlab)),
           limitsize = FALSE)
  }

  mod_out <- if (do_coord) modtab else modtab[, !grepl("^region_|kME|dir_", names(modtab))]
  write_xlsx_safe(list(Summary = cov, Modules = mod_out, DMR_to_module = dmr_out,
                       DMR_region_hits = hits[order(hits$dmr_id), ],
                       Enriched_x_trait = data.frame(DMR_enriched = rownames(cross), as.data.frame.matrix(cross),
                                                     check.names = FALSE),
                       Module_all_traits = if (nrow(prof)) prof else NULL),
                  file.path(odir, paste0(tag, ".xlsx")))

  # ---- plots ----
  pm <- mm
  if (nrow(pm)) {
    pm$fill <- ifelse(col_ok(pm$module), pm$module, "grey50")
    pm$enr <- pm$DMR_enriched %in% TRUE
    pm$solid <- pm$enr | pm$n_DMRs >= opt$min_dmrs_plot
    pm$lab <- ifelse(pm$enr | (pm$trait_associated %in% TRUE & pm$n_DMRs >= opt$min_dmrs_plot),
                     paste0(pm$module, ifelse(pm$enr, "*", ""), " (", pm$n_DMRs, ")"), "")
    pm$shape <- ifelse(pm$trait_associated %in% TRUE, "trait-associated (ME-trait p < cutoff)", "not trait-associated")
    p <- ggplot(pm, aes(ME_trait_bicor, enrich_log2OR)) +
      geom_hline(yintercept = 0, colour = "grey70") + geom_vline(xintercept = 0, colour = "grey70") +
      geom_point(aes(size = n_DMRs, shape = shape, alpha = solid), fill = pm$fill, colour = "black",
                 stroke = ifelse(pm$enr, 1.2, .3)) +
      scale_shape_manual(values = c("trait-associated (ME-trait p < cutoff)" = 24, "not trait-associated" = 21)) +
      scale_size(range = c(2, 10), name = "DMRs in module") +
      scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = .3), guide = "none") +
      labs(x = paste0("Module eigengene bicor with ", trait_lab, " (", tp, ")"),
           y = paste0("log2 OR (Haldane-corrected), DMR enrichment (", lvl, " level)"), shape = NULL,
           title = paste(sid, "DMRs vs", nid, "modules"),
           subtitle = paste0("* = enrichment FDR < ", opt$fdr, " (thick outline); label = module (n DMRs)",
                             if (!is.null(rho)) sprintf("\nSpearman(log2 OR, |bicor|) = %.2f, p = %.2g, over %d modules with >= %d DMRs",
                                                        rho$estimate, rho$p.value, nrow(mr), opt$min_dmrs_plot)),
           caption = paste0("Faded = fewer than ", opt$min_dmrs_plot, " DMRs and not significant. DMRs cover a tiny ",
                            "fraction of network regions,\nso even 1 DMR in a small module gives a large OR.")) +
      theme_bw(11) + theme(plot.caption = element_text(hjust = 0, size = 8, colour = "grey30"),
                           plot.subtitle = element_text(size = 9))
    p <- p + scale_x_continuous(expand = expansion(mult = .18)) + scale_y_continuous(expand = expansion(mult = .1)) +
      if (has_repel) ggrepel::geom_text_repel(aes(label = lab), size = 3.2, max.overlaps = Inf, seed = 1,
                                              box.padding = 1, point.padding = .8, min.segment.length = 0,
                                              segment.colour = "grey50") else
      geom_text(aes(label = lab), size = 3.2, vjust = -1.6)
    ggsave(file.path(odir, paste0(tag, "_enrichment_vs_trait.pdf")), p, width = 8.5, height = 6.8)

    # ---- genes behind the DMR hits in enriched / trait-associated modules ----
    sel <- pm$module[pm$enr | (pm$trait_associated %in% TRUE & pm$n_DMRs >= opt$min_dmrs_plot)]
    if (opt$module_genes > 0 && do_coord && length(sel) && nrow(hits)) {
      hx <- hits[hits$module %in% sel, ]
      hx$gene <- ifelse(!is.na(hx$dmr_gene) & hx$dmr_gene != "", hx$dmr_gene, hx$gene_symbol)
      hx <- hx[order(hx$dmr_id, hx$module, -ifelse(is.na(hx$membership), -Inf, hx$membership)), ]
      gx <- hx[!duplicated(hx[, c("dmr_id", "module")]), ]
      gx$direction <- ifelse(is.na(gx$expected_dir), "no ME-trait value",
                      ifelse(gx$observed_dir == gx$expected_dir, "matches ME-trait direction", "opposite"))
      gx$rank <- ave(-ifelse(is.na(gx$membership), 0, gx$membership) - abs(gx$dmr_effect) / 1e6, gx$module,
                     FUN = function(v) rank(v, ties.method = "first"))
      gx$lab <- ifelse(gx$rank <= opt$module_genes & !is.na(gx$gene), gx$gene, "")
      gx$facet <- paste0(gx$module, "  (", ave(gx$dmr_id, gx$module, FUN = length), " DMRs; bicor ",
                         round(gx$ME_trait_bicor, 2), ")")
      xvar <- if (all(is.na(gx$membership))) "rank" else "membership"
      p4 <- ggplot(gx, aes(.data[[xvar]], dmr_effect, colour = direction)) +
        geom_hline(yintercept = 0, colour = "grey70") + geom_point(size = 2.4) +
        facet_wrap(~facet, scales = "free_x") +
        scale_colour_manual(values = c("matches ME-trait direction" = "#1A9850", "opposite" = "#D73027",
                                       "no ME-trait value" = "grey50"), name = NULL) +
        labs(x = if (xvar == "membership") paste0("Module membership (kME, ", tp, ") of the overlapping region") else "Rank",
             y = "DMR methylation difference (DMRichR)",
             title = paste(sid, "- genes at DMRs in", nid, "modules"),
             subtitle = paste0("Labeled: top ", opt$module_genes, " per module by kME (most central first)")) +
        theme_bw(11) + theme(legend.position = "bottom")
      p4 <- p4 + if (has_repel) ggrepel::geom_text_repel(aes(label = lab), size = 3, fontface = "italic",
                                                         colour = "black", max.overlaps = Inf,
                                                         min.segment.length = 0, seed = 1) else
        geom_text(aes(label = lab), size = 3, vjust = -0.8, colour = "black")
      nf <- length(unique(gx$facet))
      ggsave(file.path(odir, paste0(tag, "_module_DMR_genes.pdf")), p4,
             width = min(16, 4.5 * min(nf, 3) + 1), height = 4.2 * ceiling(nf / 3) + 1.2, limitsize = FALSE)
      write_xlsx_safe(list(Module_DMR_genes = gx[order(gx$module, gx$rank),
                             intersect(c("module", "rank", "gene", "dmr_id", "dmr_effect", "membership", "hubRegion",
                                         "RegionID", "chr", "start", "end", "ME_trait_bicor", "direction",
                                         "genomic_location"), names(gx))]),
                      file.path(odir, paste0(tag, "_module_DMR_genes.xlsx")))
    }

    dc <- if (do_coord) pm[pm$n_dir_tested >= opt$min_dmrs_concordance, ] else pm[0, ]
    if (nrow(dc)) {
      dc$module <- factor(dc$module, dc$module[order(dc$pct_dir_concordant)])
      p2 <- ggplot(dc, aes(module, pct_dir_concordant, fill = trait_associated)) +
        geom_col() + geom_hline(yintercept = 50, linetype = 2) + coord_flip() +
        geom_text(aes(label = sprintf("%d/%d", round(pct_dir_concordant * n_dir_tested / 100), n_dir_tested)),
                  hjust = -0.1, size = 3) +
        scale_fill_manual(values = c(`TRUE` = "#2C7BB6", `FALSE` = "#BDBDBD"), name = "Trait-associated") +
        scale_y_continuous(limits = c(0, 115), breaks = seq(0, 100, 25)) +
        labs(x = NULL, y = "% DMRs with direction = sign(ME-trait bicor) x sign(kME)",
             title = paste(sid, "- direction concordance by module"),
             subtitle = "Expected DMR direction = sign(bicor) x sign(kME); 50% = chance") +
        theme_bw(11)
      ggsave(file.path(odir, paste0(tag, "_direction_concordance.pdf")), p2, width = 7, height = 5)
    }
    enr <- if (do_coord && any(!is.na(reg$membership))) pm$module[pm$DMR_enriched %in% TRUE] else character(0)
    if (length(enr)) {
      kd <- reg[reg$module %in% enr, c("module", "membership", "dmr_hit")]
      kd$group <- ifelse(kd$dmr_hit, "overlaps DMR", "no DMR")
      p3 <- ggplot(kd, aes(module, membership, fill = group)) + geom_boxplot(outlier.size = .5) +
        scale_fill_manual(values = c("overlaps DMR" = "#D7191C", "no DMR" = "#F0F0F0")) +
        labs(x = NULL, y = paste0("Module membership (kME, ", tp, ")"), fill = NULL,
             title = paste(sid, "- kME of DMR-overlapping regions in enriched modules")) + theme_bw(11)
      ggsave(file.path(odir, paste0(tag, "_kME_enriched_modules.pdf")), p3, width = 7, height = 4.5)
    }
  }
  all_mod[[tag]] <- cbind(dmr_set = sid, network = nid, outcome = drow$outcome, timepoint = tp,
                          trait = trait_lab, modtab)
  all_sum[[tag]] <- cov
}

# ----------------------------- summaries -------------------------------------------
am <- bind_rows_any(all_mod); sum_all <- bind_rows_any(all_sum); ap <- bind_rows_any(all_prof)
sdir <- file.path(opt$out_dir, "summary"); dir.create(sdir, showWarnings = FALSE)

# heatmap: modules (rows) x DMR sets (columns); fill = log2 OR, label = bicor + n DMRs
summary_heatmap <- function(x, xvar, xlevels, file, title) {
  x <- x[x$module != opt$grey_module, ]
  keep <- unique(x$module[x$DMR_enriched %in% TRUE | (x$trait_associated %in% TRUE & x$n_DMRs >= 1)])
  if (!length(keep)) return(invisible(NULL))
  x <- x[x$module %in% keep, ]
  x$l2or <- pmax(-4, pmin(4, x$enrich_log2OR))
  x$star <- ifelse(x$DMR_enriched %in% TRUE, "*", "")
  x$box <- x$trait_associated %in% TRUE
  x$l2or[x$n_DMRs < opt$min_dmrs_plot & !(x$DMR_enriched %in% TRUE)] <- NA
  x$lab <- paste0(ifelse(is.na(x$ME_trait_bicor), "NA", sprintf("%.2f", x$ME_trait_bicor)), x$star, "\n",
                  x$n_DMRs, " DMR", ifelse(x$n_DMRs == 1, "", "s"))
  x$xx <- factor(x[[xvar]], xlevels[xlevels %in% x[[xvar]]])
  p <- ggplot(x, aes(xx, module, fill = l2or)) + geom_tile(colour = "white") +
    geom_tile(data = x[x$box, ], colour = "black", linewidth = .7, fill = NA) +
    geom_text(aes(label = lab, colour = !is.na(l2or) & abs(l2or) > 2.5), size = 2.8, lineheight = .85,
              show.legend = FALSE) +
    scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = "black")) +
    scale_fill_gradient2(low = "#2C7BB6", mid = "white", high = "#D7191C", name = "log2 OR",
                         na.value = "grey92") +
    labs(x = NULL, y = "Consensus module", title = title,
         subtitle = paste0("Label = ME bicor with the outcome trait (that timepoint), then DMRs in module; * enrichment FDR < ",
                           opt$fdr, "\nBox = ME-trait p < ", opt$module_trait_p, "; grey = fewer than ",
                           opt$min_dmrs_plot, " DMRs (OR not estimable). Rows: modules enriched or trait-associated in >= 1 column")) +
    theme_minimal(11) + theme(axis.text.x = element_text(angle = 30, hjust = 1),
                              plot.subtitle = element_text(size = 8.5))
  ggsave(file, p, width = max(6.5, 1.3 * length(unique(x$xx)) + 3.5), height = 2.2 + .45 * length(keep),
         limitsize = FALSE)
}

persist <- NULL
if (nrow(am)) {
  for (o in unique(am$outcome))
    summary_heatmap(am[am$outcome == o, ], "timepoint", timepoints,
                    file.path(sdir, paste0("heatmap_", o, "_across_timepoints.pdf")),
                    paste0(o, " DMRs in consensus modules, across timepoints (trait: ", tmap[o], ")"))
  for (tp in unique(am$timepoint))
    summary_heatmap(am[am$timepoint == tp, ], "outcome", sort(unique(am$outcome)),
                    file.path(sdir, paste0("heatmap_", clean_id(tp), "_across_outcomes.pdf")),
                    paste0(tp, ": DMRs of each outcome in consensus modules"))

  # persistence: per outcome x module, enrichment / trait association at each timepoint
  keep_cols <- c("n_DMRs", "enrich_log2OR", "enrich_FDR", "DMR_enriched", "ME_trait_bicor", "ME_trait_p",
                 "trait_associated", "pct_dir_concordant")
  persist <- bind_rows_any(lapply(unique(am$outcome), function(o) {
    x <- am[am$outcome == o & am$module != opt$grey_module, c("module", "timepoint", keep_cols)]
    if (!nrow(x)) return(NULL)
    tps <- timepoints[timepoints %in% x$timepoint]
    w <- Reduce(function(a, b) merge(a, b, by = "module", all = TRUE),
                lapply(tps, function(t) {
                  y <- x[x$timepoint == t, c("module", keep_cols)]
                  names(y)[-1] <- paste0(names(y)[-1], "_", clean_id(t)); y
                }))
    en <- sapply(tps, function(t) w[[paste0("DMR_enriched_", clean_id(t))]] %in% TRUE)
    ta <- sapply(tps, function(t) w[[paste0("trait_associated_", clean_id(t))]] %in% TRUE)
    if (is.null(dim(en))) { en <- matrix(en, nrow = nrow(w)); ta <- matrix(ta, nrow = nrow(w)) }
    w$n_timepoints_enriched <- rowSums(en)
    w$enriched_at <- apply(en, 1, function(v) paste(tps[v], collapse = ";"))
    w$trait_associated_at <- apply(ta, 1, function(v) paste(tps[v], collapse = ";"))
    w$enriched_and_trait_at <- apply(en & ta, 1, function(v) paste(tps[v], collapse = ";"))
    w <- w[w$n_timepoints_enriched > 0 | w$trait_associated_at != "", ]
    if (!nrow(w)) return(NULL)
    cbind(outcome = o, trait = unname(tmap[o]), w[order(-w$n_timepoints_enriched, w$module), ])
  }))
}

write_xlsx_safe(list(Summary = sum_all, Modules_all = am, Persistence_by_outcome = persist,
                     Enriched_module_all_traits = ap),
                file.path(opt$out_dir, "dmrichr_vs_comethyl_summary.xlsx"))
if (!is.null(persist) && nrow(persist)) {
  message("\nModules DMR-enriched at >= 1 timepoint:")
  y <- persist[persist$n_timepoints_enriched > 0, ]
  for (r in seq_len(nrow(y))) message(sprintf("  %-20s %-15s enriched: %-30s trait-assoc: %s", y$outcome[r],
                                              y$module[r], y$enriched_at[r], y$trait_associated_at[r]))
}
message("\nDone. Outputs in ", normalizePath(opt$out_dir))