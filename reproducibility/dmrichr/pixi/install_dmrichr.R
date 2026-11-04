options(repos = c(CRAN = "https://cloud.r-project.org"))
Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS = "true")

if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}

# Install missing non-conda packages if needed
cran_pkgs <- c("hablar")
for (pkg in cran_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
}

required_pkgs <- c(
  "GenomeInfoDbData",
  "GO.db",
  "org.Hs.eg.db",
  "TxDb.Hsapiens.UCSC.hg19.knownGene",
  "TxDb.Hsapiens.UCSC.hg38.knownGene",
  "ChIPseeker",
  "lsmeans",
  "R2HTML",
  "hablar",
  "wesanderson",
  "PCAtools"
)

missing <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]

if (length(missing) > 0) {
  stop(
    "These required packages are missing from the Pixi environment: ",
    paste(missing, collapse = ", ")
  )
}

remotes::install_github(
  "ben-laufer/DMRichR",
  upgrade = "never",
  dependencies = FALSE
)

if (requireNamespace("DMRichR", quietly = TRUE)) {
  message("DMRichR installation complete")
} else {
  stop("DMRichR installation failed")
}