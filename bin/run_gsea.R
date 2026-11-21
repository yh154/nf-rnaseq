#!/usr/bin/env Rscript
# fgsea per contrast and GMT collection, ranked by the unshrunken DESeq2 Wald statistic.
VERSION <- "0.1.0"
HELP <- "Usage: run_gsea.R --results_dir DIR --contrasts FILE --gmt FILE[,FILE...] --outdir DIR
                  [--gmt_id_type gene_name|gene_id] [--min_set_size 5] [--max_set_size 500]

Pre-ranked gene-set enrichment (fgsea, fgseaMultilevel) for every contrast with
run_deseq2 = true and every GMT file.

Ranking metric: the DESeq2 Wald `stat` column of Deseq2_res_<treatment>_vs_<control>.csv
(unshrunken; genes with NA stat are dropped). The sign follows treatment/control, so a
positive NES means enriched in <treatment>. Gene ids are taken from the gene_name or gene_id
column (--gmt_id_type); when several genes share a name the one with the largest |stat| is kept.

Arguments
  --results_dir DIR    directory with the DESeq2 result CSV files
  --contrasts FILE     validated contrasts (id, variable, treatment, control, run_deseq2)
  --gmt LIST           comma-separated GMT files
  --outdir DIR         output directory (created)
  --gmt_id_type TYPE   gene_name (default) or gene_id
  --min_set_size N     smallest gene set tested, default 5
  --max_set_size N     largest gene set tested, default 500
  -h, --help           show this help and exit 0
  --version            print the version and exit 0

Output: GSEA_<treatment>_vs_<control>_<gmt>.csv per contrast x GMT with pathway, pval, padj,
log2err, ES, NES, size, leadingEdge (genes joined by ';'), sorted by padj then pval.

Exit codes: 0 ok, 1 failure, 2 usage error.

Example
  run_gsea.R --results_dir diffexp --contrasts contrasts.valid.csv --gmt h.all.gmt \\
      --gmt_id_type gene_name --outdir GSEA
"
usage_error <- function(msg) {
  cat("ERROR: ", msg, "\n\n", file = stderr())
  cat(strsplit(HELP, "\n\n", fixed = TRUE)[[1]][1], "\n", file = stderr())
  quit(status = 2)
}
args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 1 && args[1] %in% c("-h", "--help")) { cat(HELP); quit(status = 0) }
if (length(args) >= 1 && args[1] == "--version") { cat("run_gsea.R", VERSION, "\n"); quit(status = 0) }
opt <- list(gmt_id_type = "gene_name", min_set_size = "5", max_set_size = "500")
valid <- c("results_dir", "contrasts", "gmt", "outdir", "gmt_id_type", "min_set_size", "max_set_size")
i <- 1
while (i <= length(args)) {
  key <- sub("^--", "", args[i])
  if (!startsWith(args[i], "--") || !(key %in% valid) || i == length(args))
    usage_error(paste("bad argument:", args[i]))
  opt[[key]] <- args[i + 1]
  i <- i + 2
}
for (r in c("results_dir", "contrasts", "gmt", "outdir"))
  if (is.null(opt[[r]])) usage_error(paste0("--", r, " is required"))
if (!(opt$gmt_id_type %in% c("gene_name", "gene_id"))) usage_error("--gmt_id_type must be gene_name or gene_id")
min_size <- suppressWarnings(as.integer(opt$min_set_size)); max_size <- suppressWarnings(as.integer(opt$max_set_size))
if (is.na(min_size) || min_size < 1 || is.na(max_size) || max_size < min_size)
  usage_error("--min_set_size / --max_set_size must be integers with 1 <= min <= max")

read_gmt <- function(path) {
  pw <- list()
  for (l in readLines(path, warn = FALSE)) {
    if (!nzchar(trimws(l))) next
    f <- strsplit(l, "\t", fixed = TRUE)[[1]]
    pw[[f[1]]] <- unique(trimws(f[-(1:2)][nzchar(trimws(f[-(1:2)]))]))
  }
  pw
}

main <- function() {
  suppressPackageStartupMessages(library(fgsea))
  dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)
  contrasts <- read.csv(opt$contrasts, stringsAsFactors = FALSE, colClasses = "character")
  gmts <- trimws(strsplit(opt$gmt, ",", fixed = TRUE)[[1]])
  cols <- c("pathway", "pval", "padj", "log2err", "ES", "NES", "size", "leadingEdge")
  for (k in seq_len(nrow(contrasts))) {
    ct <- contrasts[k, ]
    if (!identical(ct$run_deseq2, "true")) next
    r <- read.csv(file.path(opt$results_dir, sprintf("Deseq2_res_%s_vs_%s.csv", ct$treatment, ct$control)),
                  check.names = FALSE, stringsAsFactors = FALSE)
    r <- r[!is.na(r$stat), , drop = FALSE]
    r <- r[order(-abs(r$stat)), , drop = FALSE]
    r <- r[!duplicated(r[[opt$gmt_id_type]]), , drop = FALSE]
    stats <- setNames(r$stat, r[[opt$gmt_id_type]])
    for (g in gmts) {
      pw <- read_gmt(g)
      set.seed(1)
      res <- suppressWarnings(fgseaMultilevel(pw, stats, minSize = min_size, maxSize = max_size,
                                              BPPARAM = BiocParallel::SerialParam()))
      res <- as.data.frame(res)
      if (nrow(res) > 0) {
        res$leadingEdge <- vapply(res$leadingEdge, paste, character(1), collapse = ";")
        res <- res[order(res$padj, res$pval), cols, drop = FALSE]
      } else {
        res <- setNames(data.frame(matrix(ncol = length(cols), nrow = 0)), cols)
      }
      fn <- sprintf("GSEA_%s_vs_%s_%s.csv", ct$treatment, ct$control, sub("\\.gmt$", "", basename(g)))
      write.csv(res, file.path(opt$outdir, fn), row.names = FALSE, na = "NA")
      message(sprintf("%s: %d of %d gene sets tested, %d with padj < 0.05", fn, nrow(res), length(pw),
                      sum(res$padj < 0.05, na.rm = TRUE)))
    }
  }
  writeLines(paste0("fgsea: ", as.character(packageVersion("fgsea"))), file.path(opt$outdir, "versions_fgsea.yml"))
}
tryCatch(main(), error = function(e) { cat("ERROR: ", conditionMessage(e), "\n", file = stderr()); quit(status = 1) })
