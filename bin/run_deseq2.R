#!/usr/bin/env Rscript
# One DESeq2 fit on round(RSEM expected_count) and one results() per contrast.
VERSION <- "0.1.0"
HELP <- "Usage: run_deseq2.R --counts FILE --samplesheet FILE --contrasts FILE --outdir DIR
                    [--design STR] [--padj 0.05] [--min_count 10]

One DESeq2 fit on all samples, then one results() per contrast.

Method
  counts = round(expected_count) -> DESeqDataSetFromMatrix(design). Genes with a total count
  below --min_count are left out of the fit (they stay in every output table with NA
  statistics). One DESeq() fit with the full design; every contrast is extracted from it
  with results(contrast = c(variable, treatment, control), independentFiltering = TRUE,
  alpha = --padj). log2FoldChange = log2(treatment / control), unshrunken. All design
  variables are treated as factors.

Arguments
  --counts FILE       RSEM_expected_counts.txt (gene_id, gene_biotype, gene_name, samples)
  --samplesheet FILE  validated samplesheet (sample, condition, design variables)
  --contrasts FILE    validated contrasts: id, variable, treatment, control, run_deseq2
  --outdir DIR        output directory (created)
  --design STR        additive formula, default '~ condition'
  --padj FLOAT        alpha of results() and of independent filtering, default 0.05
  --min_count INT     pre-filter on the total count over all samples, default 10
  -h, --help          show this help and exit 0
  --version           print the version and exit 0

Output (in --outdir)
  dds.rds                           the fitted DESeqDataSet
  Deseq2_res_<treatment>_vs_<control>.csv   one per contrast with run_deseq2 = true;
      gene_id, gene_name, gene_biotype, normalized counts of the two groups,
      baseMean, log2FoldChange, lfcSE, stat, pvalue, padj
  DESeq2_normalized_counts.txt      size-factor normalized counts, all genes, all samples
  vst_counts.txt                    variance-stabilised counts (blind), fitted genes only (the pipeline does not publish it)
  versions_deseq2.yml               R and DESeq2 versions

Exit codes: 0 ok, 1 failure, 2 usage error.

Example
  run_deseq2.R --counts RSEM_expected_counts.txt --samplesheet samplesheet.valid.csv \\
      --contrasts contrasts.valid.csv --design '~ batch + condition' --outdir diffexp
"

usage_error <- function(msg) {
  cat("ERROR: ", msg, "\n\n", file = stderr())
  cat(strsplit(HELP, "\n\n", fixed = TRUE)[[1]][1], "\n", file = stderr())
  quit(status = 2)
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 1 && args[1] %in% c("-h", "--help")) { cat(HELP); quit(status = 0) }
if (length(args) >= 1 && args[1] == "--version") { cat("run_deseq2.R", VERSION, "\n"); quit(status = 0) }

opt <- list(design = "~ condition", padj = "0.05", min_count = "10")
valid <- c("counts", "samplesheet", "contrasts", "outdir", "design", "padj", "min_count")
i <- 1
while (i <= length(args)) {
  key <- sub("^--", "", args[i])
  if (!startsWith(args[i], "--") || !(key %in% valid) || i == length(args))
    usage_error(paste("bad argument:", args[i]))
  opt[[key]] <- args[i + 1]
  i <- i + 2
}
for (r in c("counts", "samplesheet", "contrasts", "outdir"))
  if (is.null(opt[[r]])) usage_error(paste0("--", r, " is required"))
padj <- suppressWarnings(as.numeric(opt$padj))
min_count <- suppressWarnings(as.numeric(opt$min_count))
if (is.na(padj) || padj <= 0 || padj > 1) usage_error("--padj must be in (0, 1]")
if (is.na(min_count) || min_count < 0) usage_error("--min_count must be >= 0")

main <- function() {
  suppressPackageStartupMessages(library(DESeq2))
  dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)

  cnt <- read.table(opt$counts, header = TRUE, sep = "\t", check.names = FALSE, quote = "",
                    comment.char = "", stringsAsFactors = FALSE)
  anno <- cnt[, c("gene_id", "gene_biotype", "gene_name")]
  samples <- setdiff(colnames(cnt), colnames(anno))
  counts <- as.matrix(cnt[, samples, drop = FALSE])
  rownames(counts) <- cnt$gene_id
  counts <- round(counts)
  storage.mode(counts) <- "integer"

  ss <- read.csv(opt$samplesheet, check.names = FALSE, stringsAsFactors = FALSE)
  meta <- ss[!duplicated(ss$sample), , drop = FALSE]
  meta <- meta[match(samples, meta$sample), , drop = FALSE]
  if (anyNA(meta$sample)) stop("samples of the counts table are missing from the samplesheet")
  rownames(meta) <- meta$sample
  design <- as.formula(opt$design)
  for (v in all.vars(design)) {
    if (!(v %in% colnames(meta))) stop("design variable not in samplesheet: ", v)
    meta[[v]] <- factor(meta[[v]])
  }
  contrasts <- read.csv(opt$contrasts, stringsAsFactors = FALSE, colClasses = "character")

  dds <- DESeqDataSetFromMatrix(countData = counts, colData = meta, design = design)
  keep <- rowSums(counts(dds)) >= min_count
  message(sprintf("genes: %d total, %d with total count >= %g are fitted", nrow(dds), sum(keep), min_count))
  if (sum(keep) < 2) stop("fewer than 2 genes pass --min_count")
  dds <- DESeq(dds[keep, ], quiet = TRUE)
  saveRDS(dds, file.path(opt$outdir, "dds.rds"))

  norm_all <- sweep(counts, 2, sizeFactors(dds)[colnames(counts)], "/")
  write_table <- function(m, genes, file) {
    out <- data.frame(anno[match(genes, anno$gene_id), ], m[genes, , drop = FALSE], check.names = FALSE)
    write.table(out, file, sep = "\t", quote = FALSE, row.names = FALSE)
  }
  write_table(norm_all, anno$gene_id, file.path(opt$outdir, "DESeq2_normalized_counts.txt"))

  vsd <- if (nrow(dds) >= 1000) vst(dds, blind = TRUE) else varianceStabilizingTransformation(dds, blind = TRUE)
  write_table(assay(vsd), rownames(dds), file.path(opt$outdir, "vst_counts.txt"))

  base_all <- rowMeans(norm_all)
  for (k in seq_len(nrow(contrasts))) {
    ct <- contrasts[k, ]
    if (!identical(ct$run_deseq2, "true")) {
      message("skipping ", ct$id, " (run_deseq2 = false)")
      next
    }
    res <- results(dds, contrast = c(ct$variable, ct$treatment, ct$control),
                   independentFiltering = TRUE, alpha = padj)
    rdf <- as.data.frame(res)[match(anno$gene_id, rownames(res)), , drop = FALSE]
    rdf$baseMean[is.na(rdf$baseMean)] <- base_all[is.na(rdf$baseMean)]
    sel <- samples[as.character(meta[[ct$variable]]) %in% c(ct$treatment, ct$control)]
    out <- data.frame(anno[, c("gene_id", "gene_name", "gene_biotype")],
                      norm_all[, sel, drop = FALSE],
                      rdf[, c("baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")],
                      check.names = FALSE, row.names = NULL)
    fn <- sprintf("Deseq2_res_%s_vs_%s.csv", ct$treatment, ct$control)
    write.csv(out, file.path(opt$outdir, fn), row.names = FALSE, na = "NA")
    sig <- !is.na(out$padj) & out$padj < padj
    message(sprintf("%s: %d genes tested, padj < %g: %d up, %d down (log2FC = %s/%s)",
                    fn, sum(!is.na(out$padj)), padj, sum(sig & out$log2FoldChange > 0),
                    sum(sig & out$log2FoldChange < 0), ct$treatment, ct$control))
  }
  writeLines(c(paste0("R: ", R.version.string), paste0("DESeq2: ", as.character(packageVersion("DESeq2")))),
             file.path(opt$outdir, "versions_deseq2.yml"))
}

tryCatch(main(), error = function(e) {
  cat("ERROR: ", conditionMessage(e), "\n", file = stderr())
  quit(status = 1)
})
