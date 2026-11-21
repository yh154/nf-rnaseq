#!/usr/bin/env Rscript
# GSVA scores per GMT collection and a per-gene-set test per contrast.
VERSION <- "0.1.0"
HELP <- "Usage: run_gsva.R --expr FILE --input_type vst|log2tpm --samplesheet FILE --contrasts FILE
                  --gmt FILE[,FILE...] --outdir DIR [--gmt_id_type gene_name|gene_id]
                  [--min_set_size 5] [--max_set_size 500]

Gene-set variation analysis (GSVA >= 2.0, gsvaParam, Gaussian kernel) of every GMT collection,
then a two-group test of the scores for every contrast with run_deseq2 = true.

Arguments
  --expr FILE          expression table: gene_id, gene_biotype, gene_name, <samples>
                       (vst_counts.txt, or TPM.txt with --input_type log2tpm)
  --input_type TYPE    vst: used as is; log2tpm: log2(value + 1) is applied first
  --samplesheet FILE   validated samplesheet (sample + the contrast variables)
  --contrasts FILE     validated contrasts (id, variable, treatment, control, run_deseq2)
  --gmt LIST           comma-separated GMT files
  --outdir DIR         output directory (created)
  --gmt_id_type TYPE   gene_name (default) or gene_id; duplicated names keep the row with the
                       highest mean expression
  --min_set_size N     smallest gene set scored, default 5
  --max_set_size N     largest gene set scored, default 500
  -h, --help / --version

A Welch t-test of the scores of <treatment> vs <control> samples is reported (base R; pvalue, padj;
the table is sorted by it). It needs at least 2 samples per group (guaranteed: single samples are
skipped); a gene set whose scores are constant gives NA. The design covariates (e.g. batch) are NOT
adjusted for here.

Output
  GSVA_scores_<gmt>.csv                          pathway + one column per sample
  GSVA_test_<treatment>_vs_<control>_<gmt>.csv   pathway, mean_treatment, mean_control,
                                                 mean_diff (treatment - control), pvalue, padj (BH)
                                                 [Welch t-test]

Exit codes: 0 ok, 1 failure, 2 usage error.

Example
  run_gsva.R --expr vst_counts.txt --input_type vst --samplesheet samplesheet.valid.csv \\
      --contrasts contrasts.valid.csv --gmt h.all.gmt --outdir GSVA
"
usage_error <- function(msg) {
  cat("ERROR: ", msg, "\n\n", file = stderr())
  cat(strsplit(HELP, "\n\n", fixed = TRUE)[[1]][1], "\n", file = stderr())
  quit(status = 2)
}
args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 1 && args[1] %in% c("-h", "--help")) { cat(HELP); quit(status = 0) }
if (length(args) >= 1 && args[1] == "--version") { cat("run_gsva.R", VERSION, "\n"); quit(status = 0) }
opt <- list(gmt_id_type = "gene_name", min_set_size = "5", max_set_size = "500", input_type = "vst")
valid <- c("expr", "input_type", "samplesheet", "contrasts", "gmt", "outdir", "gmt_id_type", "min_set_size", "max_set_size")
i <- 1
while (i <= length(args)) {
  key <- sub("^--", "", args[i])
  if (!startsWith(args[i], "--") || !(key %in% valid) || i == length(args))
    usage_error(paste("bad argument:", args[i]))
  opt[[key]] <- args[i + 1]
  i <- i + 2
}
for (r in c("expr", "samplesheet", "contrasts", "gmt", "outdir"))
  if (is.null(opt[[r]])) usage_error(paste0("--", r, " is required"))
if (!(opt$gmt_id_type %in% c("gene_name", "gene_id"))) usage_error("--gmt_id_type must be gene_name or gene_id")
if (!(opt$input_type %in% c("vst", "log2tpm"))) usage_error("--input_type must be vst or log2tpm")
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
  suppressPackageStartupMessages(library(GSVA))
  dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)
  tab <- read.table(opt$expr, header = TRUE, sep = "\t", check.names = FALSE, quote = "",
                    comment.char = "", stringsAsFactors = FALSE)
  samples <- setdiff(colnames(tab), c("gene_id", "gene_biotype", "gene_name"))
  m <- as.matrix(tab[, samples, drop = FALSE])
  if (opt$input_type == "log2tpm") m <- log2(m + 1)
  ids <- tab[[opt$gmt_id_type]]
  o <- order(-rowMeans(m)); m <- m[o, , drop = FALSE]; ids <- ids[o]
  keep <- !duplicated(ids)
  m <- m[keep, , drop = FALSE]; rownames(m) <- ids[keep]

  ss <- read.csv(opt$samplesheet, check.names = FALSE, stringsAsFactors = FALSE)
  meta <- ss[!duplicated(ss$sample), , drop = FALSE]
  rownames(meta) <- meta$sample
  contrasts <- read.csv(opt$contrasts, stringsAsFactors = FALSE, colClasses = "character")

  for (g in trimws(strsplit(opt$gmt, ",", fixed = TRUE)[[1]])) {
    name <- sub("\\.gmt$", "", basename(g))
    pw <- read_gmt(g)
    param <- gsvaParam(m, pw, minSize = min_size, maxSize = max_size, kcdf = "Gaussian")
    sc <- suppressMessages(gsva(param, verbose = FALSE))
    write.csv(data.frame(pathway = rownames(sc), sc, check.names = FALSE),
              file.path(opt$outdir, sprintf("GSVA_scores_%s.csv", name)), row.names = FALSE)
    for (k in seq_len(nrow(contrasts))) {
      ct <- contrasts[k, ]
      if (!identical(ct$run_deseq2, "true")) next
      grp <- as.character(meta[colnames(sc), ct$variable])
      tr <- sc[, grp == ct$treatment, drop = FALSE]; cr <- sc[, grp == ct$control, drop = FALSE]
      pv <- vapply(seq_len(nrow(sc)), function(i)
        tryCatch(t.test(tr[i, ], cr[i, ])$p.value, error = function(e) NA_real_), numeric(1))
      pv[is.nan(pv)] <- NA_real_
      out <- data.frame(pathway = rownames(sc), mean_treatment = rowMeans(tr), mean_control = rowMeans(cr),
                        mean_diff = rowMeans(tr) - rowMeans(cr), pvalue = pv, padj = p.adjust(pv, "BH"))
      out <- out[order(out$pvalue, -abs(out$mean_diff)), ]
      write.csv(out, file.path(opt$outdir, sprintf("GSVA_test_%s_vs_%s_%s.csv", ct$treatment, ct$control, name)),
                row.names = FALSE, na = "NA")
    }
    message(sprintf("%s: %d of %d gene sets scored", name, nrow(sc), length(pw)))
  }
  writeLines(paste0("GSVA: ", as.character(packageVersion("GSVA"))), file.path(opt$outdir, "versions_gsva.yml"))
}
tryCatch(main(), error = function(e) { cat("ERROR: ", conditionMessage(e), "\n", file = stderr()); quit(status = 1) })
