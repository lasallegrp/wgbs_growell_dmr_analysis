#!/usr/bin/env Rscript

# ============================================================
# 15_Annotate_modules_consensus.R
#
# Consensus-pipeline counterpart of 12a_annotate_modules.R. The
# annotation, gene assignment, enrichment, plots and output files are
# IDENTICAL to 12a (the code below is copied verbatim from 12a apart from
# input handling and output paths), so consensus and standard results
# can be compared directly and fed to the same downstream code.
#
# Differences from 12a:
#   - input is one Consensus_Region_Assignments.tsv (script 09 shared/)
#     instead of module RDS objects -- consensus module assignments are
#     shared across datasets, so annotation runs once per adjustment
#     version
#   - output goes under comethyl_output/consensus/15_module_annotation/
#
# What it does (same as 12a):
#   - annotatr genic / FANTOM5 enhancer / CpG context for ALL regions
#     (grey included); every overlapping label + one primary label
#   - gene assignment = comethyl::annotateModule()'s GREAT rule (DEFAULT):
#     basalPlusExt, 5 kb up / 1 kb down, extension up to 1 Mb, up to two
#     genes per region ("A | B"), offline via rGREAT with GREAT v4's hg38
#     gene table. --gene_assignment nearest = single nearest gene body.
#   - per-module Fisher's tests of genic/CpG categories, BH-FDR per family
#   - Region_Gene_Associations.tsv, gene lists, background genes,
#     distribution + enrichment tables, plots, run logs
#
# Build the annotatr cache ONCE on a login node (internet needed unless
# local files are given), exactly as for 12a -- the same cache works for
# both pipelines:
#   pixi run Rscript 15_Annotate_modules_consensus.R --build_cache_only TRUE \
#     --genome hg38 --annotatr_cache /path/reference/annotatr_hg38.rds
#
# Typical run:
#   Rscript 15_Annotate_modules_consensus.R \
#     --project_root /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_DBS_dmrich_George \
#     --consensus_regions_tsv .../09_consensus_modules/cov3_75pct/covMin4_methSD0p08/v1_all_pcs/shared/Consensus_Region_Assignments.tsv \
#     --annotatr_cache /path/reference/annotatr_hg38.rds \
#     --plot_modules yellow,brown
#
# Path derivation: cpg_label / region_label / adjustment_version are read
# from the TSV path. Both layouts work:
#   .../<cpg_label>/<region_label>/<adjustment_version>/shared/Consensus_Region_Assignments.tsv
#   .../shared/<cpg_label>/<region_label>/<adjustment_version>/Consensus_Region_Assignments.tsv
# --adjustment_version overrides the derived label for the output folder.
#
# Old arguments still accepted (for existing SLURM scripts):
#   --annotation_mode great|auto -> --gene_assignment great (offline GREAT;
#                    no live web call any more); offline -> nearest
#
# Output (same files as 12a):
#   <project_root>/comethyl_output/consensus/15_module_annotation/<cpg_label>/<region_label>/<adjustment_version>/
#     Annotated_Regions.tsv/.xlsx, Region_Gene_Associations.tsv,
#     Module_Gene_List.tsv, Module_Gene_Summary.tsv/.xlsx,
#     Background_Genes.txt/.tsv, Module_Annotation_Distribution.tsv/.xlsx,
#     Module_Genic_Enrichment.tsv/.xlsx, Module_CpG_Enrichment.tsv/.xlsx,
#     plots/, run_parameters.txt, run_log.txt, session_info.txt
# ============================================================

message("Starting Script 15")

suppressPackageStartupMessages({
  library(optparse)
  library(dplyr)
  library(openxlsx)
  library(GenomicRanges)
  library(GenomeInfoDb)
  library(IRanges)
  library(AnnotationDbi)
  library(stringr)
  library(ggplot2)
})

# ============================================================
# 1) CLI
# ============================================================
option_list <- list(
  make_option("--project_root", type = "character", default = NULL,
              help = "Project root directory. Outputs are written under <project_root>/comethyl_output/consensus/15_module_annotation/. Required unless --build_cache_only TRUE."),
  make_option("--consensus_regions_tsv", type = "character", default = NULL,
              help = "Consensus_Region_Assignments.tsv from script 09 shared/ (columns RegionID, chr, start, end, module). Required unless --build_cache_only TRUE."),
  make_option("--adjustment_version", type = "character", default = NULL,
              help = "Label for the output folder (e.g. v1_all_pcs) [default: derived from the --consensus_regions_tsv path]."),
  make_option("--annotation_mode", type = "character", default = NULL,
              help = "Deprecated (old script 15): great/auto -> --gene_assignment great; offline -> --gene_assignment nearest."),
  make_option("--sample_info", type = "character", default = NULL,
              help = "Optional sample info file (xlsx/csv/tsv). Not used by annotation; only its dimensions and checksum are logged."),
  make_option("--sample_id_col", type = "character", default = NULL,
              help = "Optional sample ID column name in --sample_info (logged only)."),
  make_option("--genome", type = "character", default = "hg38",
              help = "Genome assembly with annotatr support [default %default]."),
  make_option("--helper_file", type = "character", default = NULL,
              help = "Path to helper.R [default: helper.R next to this script, else working directory]."),
  make_option("--annotatr_cache", type = "character", default = NULL,
              help = "RDS of built annotatr annotations. Loaded if it exists; otherwise annotations are built and saved here. Strongly recommended on HPC."),
  make_option("--rebuild_annotatr_cache", type = "character", default = "FALSE",
              help = "TRUE to rebuild --annotatr_cache even if it exists [default %default]."),
  make_option("--build_cache_only", type = "character", default = "FALSE",
              help = "TRUE to build/save --annotatr_cache and exit without annotating modules. Run on a node with internet unless local files are given [default %default]."),
  make_option("--cpg_island_file", type = "character", default = NULL,
              help = "Local UCSC cpgIslandExt.txt(.gz) for the genome. If given, CpG features are built from it with annotatr's algorithm instead of downloading."),
  make_option("--enhancer_file", type = "character", default = NULL,
              help = "Local FANTOM5 enhancer BED (e.g. human_permissive_enhancers_phase_1_and_2.bed.gz, hg19). If given, used instead of downloading."),
  make_option("--enhancer_chain_file", type = "character", default = NULL,
              help = "liftOver chain applied to --enhancer_file (e.g. hg19ToHg38.over.chain.gz). Omit if the BED is already in --genome."),
  make_option("--include_enhancers", type = "character", default = "TRUE",
              help = "TRUE to annotate and test FANTOM5 enhancers [default %default]."),
  make_option("--gene_assignment", type = "character", default = "great",
              help = "How genes are assigned to regions: 'great' = comethyl's GREAT rule (up to two genes per region, offline via rGREAT); 'nearest' = single nearest gene body [default %default]."),
  make_option("--tss_source", type = "character", default = "GREAT:hg38",
              help = "Gene/TSS set for --gene_assignment great: 'GREAT:hg38' = GREAT v4's own gene table shipped with rGREAT (matches comethyl); 'TxDb:hg38' = TxDb.Hsapiens.UCSC.hg38.knownGene protein-coding genes [default %default]."),
  make_option("--great_mode", type = "character", default = "basalPlusExt",
              help = "GREAT association rule: basalPlusExt (comethyl default), twoClosest or oneClosest [default %default]."),
  make_option("--basal_upstream", type = "double", default = 5000,
              help = "basalPlusExt basal domain upstream of TSS, bp [default %default]."),
  make_option("--basal_downstream", type = "double", default = 1000,
              help = "basalPlusExt basal domain downstream of TSS, bp [default %default]."),
  make_option("--extension", type = "double", default = 1e6,
              help = "Maximum domain extension, bp [default %default]."),
  make_option("--enrichment_background", type = "character", default = "all",
              help = "Background for Fisher's tests: 'all' = every region passed into network construction, grey included; 'assigned' = non-grey regions only (module vs rest of network) [default %default]."),
  make_option("--adjust_width", type = "character", default = "FALSE",
              help = "TRUE to add width-adjusted ORs from logistic regression (overlap ~ in_module + log(width)) as a sensitivity analysis [default %default]."),
  make_option("--plot_modules", type = "character", default = NULL,
              help = "Comma-separated modules of interest, e.g. 'yellow,brown,red,lightsteelblue1'. Produces *_selected distribution plots and forest plots."),
  make_option("--color_overrides", type = "character", default = NULL,
              help = "Module bar colour overrides, e.g. 'yellow=darkgoldenrod3,white=grey80'."),
  make_option("--log2OR_cap", type = "double", default = 3,
              help = "Cap for the heatmap log2(OR) colour scale; values beyond are squished [default %default]."),
  make_option("--fdr_cut", type = "double", default = 0.05,
              help = "FDR threshold for significance marks in plots [default %default].")
)
opt <- parse_args(OptionParser(option_list = option_list,
                               description = "Annotate consensus comethyl modules with annotatr genic/enhancer/CpG context and test enrichment."))

to_bool <- function(x, nm) {
  v <- toupper(trimws(as.character(x)))
  if (v %in% c("TRUE", "T", "1", "YES")) return(TRUE)
  if (v %in% c("FALSE", "F", "0", "NO")) return(FALSE)
  stop("--", nm, " must be TRUE or FALSE (got '", x, "')", call. = FALSE)
}
opt$rebuild_annotatr_cache <- to_bool(opt$rebuild_annotatr_cache, "rebuild_annotatr_cache")
opt$build_cache_only       <- to_bool(opt$build_cache_only, "build_cache_only")
opt$include_enhancers      <- to_bool(opt$include_enhancers, "include_enhancers")
opt$adjust_width           <- to_bool(opt$adjust_width, "adjust_width")

# ============================================================
# 2) Helper
# ============================================================
script_dir <- tryCatch({
  f <- sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE))
  if (length(f) == 1) dirname(normalizePath(f)) else getwd()
}, error = function(e) getwd())
helper_file <- if (!is.null(opt$helper_file)) opt$helper_file else file.path(script_dir, "helper.R")
if (!file.exists(helper_file)) {
  stop("helper.R not found at: ", helper_file, "\nPass --helper_file /path/to/helper.R.", call. = FALSE)
}
source(helper_file)
if (exists("setup_publication_font", mode = "function")) setup_publication_font()

# ============================================================
# 3) Validation
# ============================================================
opt$enrichment_background <- match.arg(opt$enrichment_background, c("all", "assigned"))
if (!is.null(opt$annotation_mode)) {
  am <- match.arg(opt$annotation_mode, c("auto", "great", "offline"))
  opt$gene_assignment <- if (am == "offline") "nearest" else "great"
  message("Note: --annotation_mode ", am, " is deprecated; using --gene_assignment ", opt$gene_assignment,
          if (am != "offline") " (GREAT rule, run offline via rGREAT)" else "")
}
opt$gene_assignment <- match.arg(opt$gene_assignment, c("great", "nearest"))
opt$great_mode <- match.arg(opt$great_mode, c("basalPlusExt", "twoClosest", "oneClosest"))
if (opt$gene_assignment == "great") invisible(great_parse_tss_source(opt$tss_source, opt$genome))
for (nm in c("cpg_island_file", "enhancer_file", "enhancer_chain_file")) {
  if (!is.null(opt[[nm]])) validate_file_exists(opt[[nm]], paste0("--", nm))
}
if (!is.null(opt$enhancer_chain_file) && is.null(opt$enhancer_file)) {
  stop("--enhancer_chain_file given without --enhancer_file.", call. = FALSE)
}
if (opt$build_cache_only && is.null(opt$annotatr_cache)) {
  stop("--build_cache_only TRUE requires --annotatr_cache.", call. = FALSE)
}
if (!is.finite(opt$log2OR_cap) || opt$log2OR_cap <= 0) stop("--log2OR_cap must be > 0.", call. = FALSE)
if (!is.finite(opt$fdr_cut) || opt$fdr_cut <= 0 || opt$fdr_cut >= 1) stop("--fdr_cut must be in (0, 1).", call. = FALSE)

plot_modules <- if (!is.null(opt$plot_modules)) trimws(strsplit(opt$plot_modules, ",")[[1]]) else NULL
color_overrides <- annot_parse_color_overrides(opt$color_overrides)

if (!opt$build_cache_only) {
  stop_if_missing(opt$project_root, "--project_root")
  if (!dir.exists(opt$project_root)) stop("project_root not found: ", opt$project_root, call. = FALSE)
  stop_if_missing(opt$consensus_regions_tsv, "--consensus_regions_tsv")
  validate_file_exists(opt$consensus_regions_tsv, "--consensus_regions_tsv")
  if (!is.null(opt$sample_info)) validate_file_exists(opt$sample_info, "--sample_info")
  AnnotationHub::setAnnotationHubOption("CACHE", value = file.path(opt$project_root, ".cache"))
}

# ============================================================
# 4) Annotations (built once, shared across variants)
# ============================================================
annotations <- build_annotatr_annotations(
  genome = opt$genome,
  include_enhancers = opt$include_enhancers,
  cache_file = opt$annotatr_cache,
  rebuild = opt$rebuild_annotatr_cache,
  cpg_island_file = opt$cpg_island_file,
  enhancer_file = opt$enhancer_file,
  enhancer_chain_file = opt$enhancer_chain_file,
  verbose = TRUE
)
annot_prov <- attr(annotations, "provenance")
if (is.null(annot_prov)) annot_prov <- list(note = "no provenance attribute on cached annotations")

if (opt$build_cache_only) {
  message("Annotation cache ready: ", opt$annotatr_cache)
  message(paste(names(annot_prov), unlist(annot_prov), sep = ": ", collapse = "\n"))
  message("Script 12a: cache build finished")
  quit(save = "no", status = 0)
}

# GREAT regulatory domains (built once, shared across variants)
extended_tss <- NULL
if (opt$gene_assignment == "great") {
  extended_tss <- great_build_extended_tss(opt$tss_source, genome = opt$genome, mode = opt$great_mode,
                                           basal_upstream = opt$basal_upstream,
                                           basal_downstream = opt$basal_downstream,
                                           extension = opt$extension, verbose = TRUE)
}

sample_info_df <- NULL
if (!is.null(opt$sample_info)) {
  ext <- tolower(tools::file_ext(opt$sample_info))
  sample_info_df <- switch(ext,
    xlsx = openxlsx::read.xlsx(opt$sample_info), xls = openxlsx::read.xlsx(opt$sample_info),
    csv = utils::read.csv(opt$sample_info, check.names = FALSE),
    utils::read.delim(opt$sample_info, check.names = FALSE))
}

# ============================================================
# 5) Consensus input + output paths
# ============================================================
# cpg_label / region_label / adjustment_version from the TSV path; a
# trailing "shared" folder is skipped, so both script-09 layouts work.
derive_consensus_dirs <- function(regions_tsv, project_root, adjustment_version = NULL,
                                  step_name = "15_module_annotation") {
  d <- dirname(normalizePath(regions_tsv))
  if (basename(d) == "shared") d <- dirname(d)
  info <- list(variant_name = basename(d), region_label = basename(dirname(d)),
               cpg_label = basename(dirname(dirname(d))))
  if (!is.null(adjustment_version) && adjustment_version != info$variant_name) {
    message("Note: --adjustment_version '", adjustment_version, "' differs from the folder name in the TSV path ('",
            info$variant_name, "'); using '", adjustment_version, "' for the output folder.")
    info$variant_name <- adjustment_version
  }
  info$pipeline_root <- file.path(project_root, "comethyl_output", "consensus")
  info$step_dir <- file.path(info$pipeline_root, step_name)
  info$out_dir <- file.path(info$step_dir, info$cpg_label, info$region_label, info$variant_name)
  safe_dir_create(info$out_dir)
  info
}

read_consensus_regions <- function(path) {
  df <- utils::read.delim(path, stringsAsFactors = FALSE, check.names = FALSE)
  validate_regions_df(df, "--consensus_regions_tsv")
  df$chr <- as.character(df$chr); df$module <- as.character(df$module)
  df$start <- as.integer(df$start); df$end <- as.integer(df$end)
  # same five columns 12a takes from Modules.rds $regions
  df[, c("RegionID", "chr", "start", "end", "module")]
}

# ============================================================
# 5b) Annotation runner (verbatim from 12a run_annotation_for_variant)
# ============================================================
write_tsv_xlsx <- function(df, stem, xlsx = TRUE) {
  utils::write.table(df, paste0(stem, ".tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "")
  if (xlsx) openxlsx::write.xlsx(df, paste0(stem, ".xlsx"), rowNames = FALSE)
}

run_annotation_consensus <- function(regions, di, variant_label, source_file) {
  out_dir <- di$out_dir
  plot_dir <- file.path(out_dir, "plots"); safe_dir_create(plot_dir)
  log_file <- file.path(out_dir, "run_log.txt")
  L <- function(...) append_log(log_file, ...)

  L("Starting annotation for variant: ", variant_label, " (annotatr engine)")
  L("consensus_regions_tsv: ", source_file)
  L("cpg_label: ", di$cpg_label, "; region_label: ", di$region_label, "; adjustment_version: ", di$variant_name)
  L("out_dir: ", out_dir)
  L("genome: ", opt$genome, "; enrichment_background: ", opt$enrichment_background,
    "; include_enhancers: ", opt$include_enhancers, "; adjust_width: ", opt$adjust_width)
  for (nm in names(annot_prov)) L("annotations ", nm, ": ", annot_prov[[nm]])
  if (!is.null(sample_info_df)) L("sample_info rows: ", nrow(sample_info_df), "; cols: ", ncol(sample_info_df))

  L("Loaded regions: ", nrow(regions), "; modules (incl. grey): ", length(unique(regions$module)))

  # --- annotatr genic / enhancer / CpG context ---
  annotated <- annotate_regions_annotatr(regions, annotations, genome = opt$genome, verbose = TRUE)

  # --- gene assignment ---
  if (opt$gene_assignment == "great") {
    assoc <- great_region_gene_associations(regions, extended_tss)
    ginfo <- tryCatch(gene_info_from_orgdb(unique(assoc$gene_entrezID), opt$genome),
                      error = function(e) { L("WARNING: gene descriptions not added: ", conditionMessage(e)); NULL })
    collapsed <- great_collapse_associations(assoc, regions, ginfo)
    annotated <- cbind(annotated, collapsed[, setdiff(colnames(collapsed), "RegionID")])
    gp <- attr(extended_tss, "great_params")
    L("Gene assignment: GREAT ", gp$mode, " (", gp$tss_source, "; basal ", gp$basal_upstream, " up / ",
      gp$basal_downstream, " down; extension ", gp$extension, "; rGREAT ", gp$rGREAT_version, "; ",
      gp$n_genes, " genes)")
    tab <- table(factor(pmin(collapsed$n_genes, 3), levels = 0:3, labels = c("0", "1", "2", "3+")))
    L("Genes per region: ", paste(names(tab), tab, sep = "=", collapse = "; "))
  } else {
    gr <- GenomicRanges::GRanges(regions$chr, IRanges::IRanges(regions$start, regions$end))
    ng <- offline_nearest_gene(gr, genome = opt$genome, verbose = TRUE)
    ng <- ng[!duplicated(ng[, c("chr", "start", "end")]), ]
    annotated <- dplyr::left_join(annotated, ng[, c("chr", "start", "end", "gene_symbol", "gene_entrezID")],
                                  by = c("chr", "start", "end"))
    assoc <- annotated[!is.na(annotated$gene_symbol) & annotated$gene_symbol != "",
                       c("RegionID", "chr", "start", "end", "gene_entrezID", "gene_symbol", "module")]
    assoc$distance_to_TSS <- NA_integer_
    L("Gene assignment: single nearest gene (offline_nearest_gene)")
  }
  stopifnot(nrow(annotated) == nrow(regions))
  L("Regions with gene_context: ", sum(!is.na(annotated$gene_context)),
    "; with CpG_context: ", sum(!is.na(annotated$CpG_context)),
    "; with >=1 gene: ", length(unique(assoc$RegionID)))

  write_tsv_xlsx(annotated, file.path(out_dir, "Annotated_Regions"))
  utils::write.table(assoc[, c("RegionID", "chr", "start", "end", "module", "gene_symbol",
                               "gene_entrezID", "distance_to_TSS")],
                     file.path(out_dir, "Region_Gene_Associations.tsv"),
                     sep = "\t", quote = FALSE, row.names = FALSE, na = "")

  # --- gene lists / background genes (one row per module-gene; "A | B" split) ---
  gl <- assoc %>% dplyr::filter(!is.na(gene_symbol), gene_symbol != "") %>%
    dplyr::distinct(module, gene_symbol) %>% dplyr::arrange(module, gene_symbol)
  utils::write.table(gl, file.path(out_dir, "Module_Gene_List.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  gs <- annotated %>% dplyr::count(module, name = "n_regions") %>%
    dplyr::left_join(gl %>% dplyr::count(module, name = "n_unique_genes"), by = "module") %>%
    dplyr::mutate(n_unique_genes = dplyr::coalesce(n_unique_genes, 0L)) %>%
    dplyr::arrange(dplyr::desc(n_regions), module)
  write_tsv_xlsx(gs, file.path(out_dir, "Module_Gene_Summary"))
  bg_genes <- sort(unique(gl$gene_symbol))
  write_lines_safe(bg_genes, file.path(out_dir, "Background_Genes.txt"))
  utils::write.table(data.frame(gene_symbol = bg_genes), file.path(out_dir, "Background_Genes.tsv"),
                     sep = "\t", quote = FALSE, row.names = FALSE)

  # --- distribution + enrichment ---
  dist <- summarizeAnnotationDistribution(annotated)
  write_tsv_xlsx(dist, file.path(out_dir, "Module_Annotation_Distribution"))

  enr <- list()
  for (fam in c("genic", "cpg")) {
    enr[[fam]] <- moduleAnnotatrEnrichment(annotated, family = fam,
                                           background = opt$enrichment_background,
                                           adjust_width = opt$adjust_width, verbose = TRUE)
    stem <- file.path(out_dir, if (fam == "genic") "Module_Genic_Enrichment" else "Module_CpG_Enrichment")
    write_tsv_xlsx(enr[[fam]], stem)
    L("Wrote ", basename(stem), " (", nrow(enr[[fam]]), " tests; ",
      sum(enr[[fam]]$fdr < opt$fdr_cut), " at FDR < ", opt$fdr_cut, ")")
  }

  # --- plots ---
  fam_name <- c(genic = "Genic", cpg = "CpG")
  for (fam in c("genic", "cpg")) {
    plotAnnotatrEnrichmentHeatmap(enr[[fam]], p_cut = opt$fdr_cut, log2OR_cap = opt$log2OR_cap,
                                  title = paste0(fam_name[[fam]], " enrichment by module (", variant_label, ")"),
                                  file_stem = file.path(plot_dir, paste0(fam_name[[fam]], "_Enrichment_Heatmap")))
    if (opt$adjust_width) {
      plotAnnotatrEnrichmentHeatmap(enr[[fam]], p_cut = opt$fdr_cut, log2OR_cap = opt$log2OR_cap,
                                    use_width_adjusted = TRUE,
                                    file_stem = file.path(plot_dir, paste0(fam_name[[fam]], "_Enrichment_Heatmap_width_adjusted")))
    }
    for (mode in c("stacked", "overlap")) {
      plotAnnotationDistribution(dist, family = fam, mode = mode, color_overrides = color_overrides,
                                 file_stem = file.path(plot_dir, paste0(fam_name[[fam]], "_Distribution_", mode, "_all")))
      if (!is.null(plot_modules)) {
        plotAnnotationDistribution(dist, family = fam, mode = mode, groups = plot_modules,
                                   color_overrides = color_overrides,
                                   file_stem = file.path(plot_dir, paste0(fam_name[[fam]], "_Distribution_", mode, "_selected")))
      }
    }
    if (!is.null(plot_modules)) {
      plotAnnotatrEnrichmentForest(enr[[fam]], groups = plot_modules, p_cut = opt$fdr_cut,
                                   color_overrides = color_overrides,
                                   file_stem = file.path(plot_dir, paste0(fam_name[[fam]], "_Enrichment_Forest_selected")))
    }
  }
  L("Wrote plots to: ", plot_dir)

  # --- provenance ---
  params <- c(
    paste0("timestamp\t", timestamp_now()),
    paste0("annotation_engine\tannotatr"),
    paste0("variant_label\t", variant_label),
    paste0("pipeline\tconsensus"),
    paste0("consensus_regions_tsv\t", source_file),
    paste0("cpg_label\t", di$cpg_label),
    paste0("region_label\t", di$region_label),
    paste0("adjustment_version\t", di$variant_name),
    paste0(names(opt), "\t", vapply(opt, function(x) if (is.null(x)) "" else paste(x, collapse = ","), "")),
    paste0("annotations_", names(annot_prov), "\t", unlist(annot_prov)),
    paste0("n_regions\t", nrow(regions)),
    paste0("n_modules_incl_grey\t", length(unique(regions$module))),
    paste0("n_background_genes\t", length(bg_genes)),
    if (!is.null(extended_tss)) paste0("great_", names(attr(extended_tss, "great_params")), "\t",
                                       unlist(attr(extended_tss, "great_params"))) else character(0)
  )
  write_lines_safe(params, file.path(out_dir, "run_parameters.txt"))
  write_session_info(file.path(out_dir, "session_info.txt"),
                     extra_files = list(consensus_regions_tsv = source_file, sample_info = opt$sample_info,
                                        annotatr_cache = opt$annotatr_cache,
                                        cpg_island_file = opt$cpg_island_file,
                                        enhancer_file = opt$enhancer_file))
  L("Finished annotation for variant: ", variant_label)
  invisible(out_dir)
}

# ============================================================
# 6) Run
# ============================================================
di <- derive_consensus_dirs(opt$consensus_regions_tsv, opt$project_root, opt$adjustment_version)
message("Derived cpg_label: ", di$cpg_label, "; region_label: ", di$region_label,
        "; adjustment_version: ", di$variant_name)
message("Output directory: ", di$out_dir)
regions <- read_consensus_regions(opt$consensus_regions_tsv)
run_annotation_consensus(regions, di, di$variant_name, opt$consensus_regions_tsv)

message("Script 15 complete: consensus module annotation finished")