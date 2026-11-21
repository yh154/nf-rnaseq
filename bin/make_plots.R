#!/usr/bin/env Rscript
# All plots of the pipeline in ONE multi-page PDF.
VERSION <- "0.1.0"
HELP <- "Usage: make_plots.R --dds FILE --contrasts FILE --results_dir DIR --out FILE
                    [--padj 0.05] [--lfc 1] [--pca_ntop 500]
                    [--gsea_dir DIR] [--gsva_dir DIR] [--gsva_top 10]

Draw every plot of the pipeline into one multi-page PDF (one plot per page, same page size),
batched by plot type, contrasts in the order of the contrasts file:
  1. PCA of ALL samples, coloured by condition (first page)
  2. PCA per contrast (the two contrasted groups only)
  3. MA plot per contrast
  4. Volcano plot per contrast
  5. padj histogram (bar plot) per contrast
  6. GSEA bubble plots   (only with --gsea_dir; per contrast x GMT collection)
  7. GSVA boxplots       (only with --gsva_dir; per contrast x GMT collection)
Page count = 1 + 4 x C + (GSEA files) + (GSVA test files), C = contrasts with a results CSV.

Direction: '<A>_vs_<B>' always means fold change A/B; Up = higher in A (treatment).
log2 fold changes are the unshrunken DESeq2 values.

Arguments
  --dds FILE          dds.rds written by run_deseq2.R
  --contrasts FILE    validated contrasts (id, variable, treatment, control, run_deseq2)
  --results_dir DIR   directory with Deseq2_res_<treatment>_vs_<control>.csv
  --out FILE          output PDF
  --padj FLOAT        adjusted p-value cut (MA colour, volcano, default 0.05)
  --lfc FLOAT         |log2FC| cut of the volcano plot only (default 1)
  --pca_ntop INT      number of most variable genes for PCA (default 500)
  --gsea_dir DIR      run_gsea.R output (GSEA_<treatment>_vs_<control>_<gmt>.csv): bubble plots,
                      top 10 positive and top 10 negative NES; positive NES = enriched in treatment
  --gsva_dir DIR      run_gsva.R output (GSVA_scores_<gmt>.csv, GSVA_test_*.csv): boxplots of the
                      --gsva_top gene sets with the smallest p-values
  --gsva_top INT      gene sets per GSVA page (default 10)
  -h, --help          show this help and exit 0
  --version           print the version and exit 0

Plots use VST (blind = TRUE) computed here from dds.rds; the output is deterministic.
Exit codes: 0 ok, 1 failure, 2 usage error.

Example
  make_plots.R --dds dds.rds --contrasts contrasts.valid.csv --results_dir diffexp --out plots.pdf
"

usage_error <- function(msg) {
  cat("ERROR: ", msg, "\n\n", file = stderr())
  cat(strsplit(HELP, "\n\n", fixed = TRUE)[[1]][1], "\n", file = stderr())
  quit(status = 2)
}
args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 1 && args[1] %in% c("-h", "--help")) { cat(HELP); quit(status = 0) }
if (length(args) >= 1 && args[1] == "--version") { cat("make_plots.R", VERSION, "\n"); quit(status = 0) }

opt <- list(padj = "0.05", lfc = "1", pca_ntop = "500", gsva_top = "10")
valid <- c("dds", "contrasts", "results_dir", "out", "padj", "lfc", "pca_ntop", "gsea_dir", "gsva_dir", "gsva_top")
i <- 1
while (i <= length(args)) {
  key <- sub("^--", "", args[i])
  if (!startsWith(args[i], "--") || !(key %in% valid) || i == length(args))
    usage_error(paste("bad argument:", args[i]))
  opt[[key]] <- args[i + 1]
  i <- i + 2
}
for (r in c("dds", "contrasts", "results_dir", "out"))
  if (is.null(opt[[r]])) usage_error(paste0("--", r, " is required"))
padj_cut <- suppressWarnings(as.numeric(opt$padj))
lfc_cut <- suppressWarnings(as.numeric(opt$lfc))
ntop <- suppressWarnings(as.integer(opt$pca_ntop))
if (is.na(padj_cut) || padj_cut <= 0 || padj_cut > 1) usage_error("--padj must be in (0, 1]")
if (is.na(lfc_cut) || lfc_cut < 0) usage_error("--lfc must be >= 0")
if (is.na(ntop) || ntop < 2) usage_error("--pca_ntop must be an integer >= 2")
gsva_top <- suppressWarnings(as.integer(opt$gsva_top))
if (is.na(gsva_top) || gsva_top < 1) usage_error("--gsva_top must be an integer >= 1")

main <- function() {
  suppressPackageStartupMessages({ library(DESeq2); library(ggplot2) })
  set.seed(1)
  dds <- readRDS(opt$dds)
  contrasts <- read.csv(opt$contrasts, stringsAsFactors = FALSE, colClasses = "character")
  contrasts <- contrasts[contrasts$run_deseq2 == "true", , drop = FALSE]
  pdf_w <- 9; pdf_h <- 7

  do_vst <- function(d) {
    f <- function(type) if (nrow(d) >= 1000) vst(d, blind = TRUE, fitType = type) else
      varianceStabilizingTransformation(d, blind = TRUE, fitType = type)
    tryCatch(suppressMessages(f("parametric")), error = function(e) suppressMessages(f("mean")))
  }

  pca_plot <- function(d, group_col, title) {
    v <- assay(do_vst(d))
    rv <- apply(v, 1, var)
    top <- order(rv, decreasing = TRUE)[seq_len(min(ntop, nrow(v)))]
    p <- prcomp(t(v[top, , drop = FALSE]))
    pct <- round(100 * p$sdev^2 / sum(p$sdev^2))
    df <- data.frame(PC1 = p$x[, 1], PC2 = if (ncol(p$x) > 1) p$x[, 2] else 0,
                     group = factor(colData(d)[[group_col]]), name = colnames(d))
    ggplot(df, aes(PC1, PC2, color = group)) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
      geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
      geom_point(size = 2.5) +
      geom_text(aes(label = name), vjust = -0.9, size = 3, show.legend = FALSE) +
      scale_x_continuous(expand = expansion(mult = 0.15)) +
      scale_y_continuous(expand = expansion(mult = 0.15)) +
      labs(title = title, x = sprintf("PC1: %d%% variance", pct[1]),
           y = sprintf("PC2: %d%% variance", if (length(pct) > 1) pct[2] else 0), color = "group") +
      theme_bw() + theme(legend.position = "right", panel.grid.minor = element_blank())
  }

  read_res <- function(ct) {
    f <- file.path(opt$results_dir, sprintf("Deseq2_res_%s_vs_%s.csv", ct$treatment, ct$control))
    if (!file.exists(f)) stop("missing results file: ", f)
    read.csv(f, check.names = FALSE, stringsAsFactors = FALSE)
  }
  cname <- function(ct) sprintf("%s_vs_%s", ct$treatment, ct$control)

  ma_plot <- function(ct) {
    r <- read_res(ct)
    r <- r[!is.na(r$log2FoldChange) & r$baseMean > 0, ]
    L <- max(2, ceiling(quantile(abs(r$log2FoldChange), 0.995, names = FALSE)))
    r$y <- pmax(pmin(r$log2FoldChange, L), -L)
    r$clipped <- abs(r$log2FoldChange) > L
    r$sig <- !is.na(r$padj) & r$padj < padj_cut
    ggplot(r, aes(baseMean, y, color = sig, shape = clipped)) +
      geom_point(size = 1.2, alpha = 0.8) +
      geom_hline(yintercept = 0, color = "black") +
      scale_x_log10() + ylim(-L, L) +
      scale_color_manual(values = c("FALSE" = "grey60", "TRUE" = "#1f4fd8"), guide = "none") +
      scale_shape_manual(values = c("FALSE" = 16, "TRUE" = 17), guide = "none") +
      labs(title = cname(ct), x = "mean of normalized counts", y = "log fold change") +
      theme_bw() + theme(panel.grid.minor = element_blank())
  }

  volcano_plot <- function(ct) {
    r <- read_res(ct)
    r <- r[!is.na(r$padj) & !is.na(r$log2FoldChange), ]
    r$p <- -log10(pmax(r$padj, 1e-300))
    sig <- r$padj < padj_cut & abs(r$log2FoldChange) > lfc_cut
    nup <- sum(sig & r$log2FoldChange > 0); ndn <- sum(sig & r$log2FoldChange < 0)
    lv <- c(sprintf("Up (%d)", nup), sprintf("Down (%d)", ndn), "NS")
    r$cls <- factor(ifelse(sig & r$log2FoldChange > 0, lv[1], ifelse(sig & r$log2FoldChange < 0, lv[2], lv[3])), levels = lv)
    lab <- r[sig, , drop = FALSE]
    lab <- head(lab[order(lab$padj, -abs(lab$log2FoldChange)), ], 20)
    cols <- setNames(c("#d62728", "#1f77b4", "grey70"), lv)
    ggplot(r, aes(log2FoldChange, p, color = cls)) +
      geom_point(size = 1.1, alpha = 0.8) +
      geom_vline(xintercept = c(-lfc_cut, lfc_cut), linetype = "dashed") +
      geom_hline(yintercept = -log10(padj_cut), linetype = "dashed") +
      geom_text(data = lab, aes(label = gene_name), vjust = -0.7, size = 2.6, check_overlap = TRUE, show.legend = FALSE) +
      scale_color_manual(values = cols, name = NULL) +
      labs(title = sprintf("%s vs %s", ct$treatment, ct$control),
           subtitle = sprintf("padj < %g, log2FC > %g: %d up, %d down", padj_cut, lfc_cut, nup, ndn),
           x = sprintf("log2 fold change (%s/%s)", ct$treatment, ct$control), y = "-log10 adjusted p") +
      theme_bw() + theme(legend.position = "top", panel.grid.minor = element_blank())
  }

  padj_bar_plot <- function(ct) {
    r <- read_res(ct)
    n_na <- sum(is.na(r$padj))
    p <- r$padj[!is.na(r$padj)]
    bin <- pmin(floor(p / 0.05), 19)
    cnt <- tabulate(bin + 1, nbins = 20)
    df <- data.frame(mid = (0:19) * 0.05 + 0.025, n = cnt)
    ggplot(df, aes(mid, n)) +
      geom_col(width = 0.05, fill = "steelblue", color = "white") +
      geom_text(aes(label = n), vjust = -0.4, size = 3) +
      scale_x_continuous(breaks = seq(0, 1, 0.1), limits = c(0, 1)) +
      scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
      labs(title = cname(ct), subtitle = sprintf("%d genes with padj = NA excluded", n_na),
           x = "Adjusted-pvalue", y = "Gene count") +
      theme_bw() + theme(panel.grid.minor = element_blank())
  }

  wrap <- function(x, w = 36) unname(vapply(x, function(z) paste(strwrap(gsub("_", " ", z), w), collapse = "\n"), ""))
  message_page <- function(title, msg) {
    print(ggplot() + annotate("text", x = 0, y = 0, label = msg) + labs(title = title) + theme_void())
  }

  gsea_plot <- function(ct, f) {
    coll <- sub("\\.csv$", "", sub(sprintf("^GSEA_%s_vs_%s_", ct$treatment, ct$control), "", basename(f)))
    ttl <- sprintf("%s vs %s", ct$treatment, ct$control)
    r <- read.csv(f, check.names = FALSE, stringsAsFactors = FALSE)
    r <- r[!is.na(r$NES), , drop = FALSE]
    if (nrow(r) == 0) return(message_page(ttl, sprintf("GSEA %s: no gene set passed the size limits", coll)))
    pos <- r[r$NES > 0, , drop = FALSE]; pos <- head(pos[order(-pos$NES), , drop = FALSE], 10)
    neg <- r[r$NES < 0, , drop = FALSE]; neg <- head(neg[order(neg$NES), , drop = FALSE], 10)
    d <- rbind(pos, neg)
    d$label <- wrap(d$pathway)
    d$label <- factor(d$label, levels = d$label[order(d$NES)])
    d$FDR <- ifelse(is.na(d$padj), 1, d$padj)
    d$absNES <- abs(d$NES)
    ggplot(d, aes(NES, label, size = absNES)) +
      geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
      geom_point(data = d[d$NES > 0, , drop = FALSE], aes(fill = FDR), shape = 21, stroke = 0) +
      geom_point(data = d[d$NES < 0, , drop = FALSE], aes(color = FDR), shape = 16) +
      scale_fill_gradient(low = "#b2182b", high = "#fddbc7", limits = c(0, 1), name = "FDR (NES > 0)") +
      scale_color_gradient(low = "#2166ac", high = "#d1e5f0", limits = c(0, 1), name = "FDR (NES < 0)") +
      scale_y_discrete(limits = levels(d$label)) +   # explicit: layers on subsets would reorder the axis
      scale_size(range = c(3, 9), name = "abs(NES)") +
      labs(title = ttl, subtitle = sprintf("GSEA %s: top and bottom 10 by NES", coll), x = "NES", y = NULL) +
      theme_bw() + theme(panel.grid.minor = element_blank(), axis.text.y = element_text(size = 8))
  }

  gsva_plot <- function(ct, f) {
    coll <- sub("\\.csv$", "", sub(sprintf("^GSVA_test_%s_vs_%s_", ct$treatment, ct$control), "", basename(f)))
    ttl <- sprintf("%s vs %s", ct$treatment, ct$control)
    tst <- read.csv(f, check.names = FALSE, stringsAsFactors = FALSE)
    sc <- read.csv(file.path(opt$gsva_dir, sprintf("GSVA_scores_%s.csv", coll)), check.names = FALSE, stringsAsFactors = FALSE)
    tst <- tst[!is.na(tst$pvalue), , drop = FALSE]
    if (nrow(tst) == 0) return(message_page(ttl, sprintf("GSVA %s: no gene set could be tested", coll)))
    top <- head(tst, gsva_top)
    grp <- as.character(colData(dds)[[ct$variable]]); names(grp) <- colnames(dds)
    smp <- names(grp)[grp %in% c(ct$treatment, ct$control)]
    long <- do.call(rbind, lapply(seq_len(nrow(top)), function(i) {
      row <- sc[sc$pathway == top$pathway[i], , drop = FALSE]
      data.frame(pathway = wrap(top$pathway[i], 28), sample = smp, score = as.numeric(row[1, smp]),
                 group = factor(grp[smp], levels = c(ct$control, ct$treatment)), stringsAsFactors = FALSE)
    }))
    long$pathway <- factor(long$pathway, levels = unique(long$pathway))
    ann <- data.frame(pathway = factor(wrap(top$pathway, 28), levels = levels(long$pathway)),
                      label = sprintf("p = %.2g", top$pvalue))
    ggplot(long, aes(group, score, color = group)) +
      geom_boxplot(outlier.shape = NA, fill = NA) +
      geom_jitter(width = 0.15, height = 0, size = 1.6) +
      geom_text(data = ann, aes(x = 1.5, y = Inf, label = label), inherit.aes = FALSE, vjust = 1.4, size = 2.8) +
      facet_wrap(~pathway, scales = "free_y") +
      labs(title = ttl, subtitle = sprintf("GSVA %s: top %d gene sets by p-value", coll, nrow(top)),
           x = ct$variable, y = "GSVA score", color = "group") +
      theme_bw() + theme(legend.position = "right", strip.text = element_text(size = 7),
                         axis.text.x = element_text(size = 7))
  }

  pdf(opt$out, width = pdf_w, height = pdf_h, onefile = TRUE)
  on.exit(dev.off(), add = TRUE)
  pages <- 0
  print(pca_plot(dds, "condition", "All samples")); pages <- pages + 1
  groups <- function(ct) {
    sel <- as.character(colData(dds)[[ct$variable]]) %in% c(ct$treatment, ct$control)
    dds[, sel]
  }
  for (k in seq_len(nrow(contrasts))) {
    ct <- contrasts[k, ]
    print(pca_plot(groups(ct), ct$variable, cname(ct))); pages <- pages + 1
  }
  for (k in seq_len(nrow(contrasts))) { print(ma_plot(contrasts[k, ])); pages <- pages + 1 }
  for (k in seq_len(nrow(contrasts))) { print(volcano_plot(contrasts[k, ])); pages <- pages + 1 }
  for (k in seq_len(nrow(contrasts))) { print(padj_bar_plot(contrasts[k, ])); pages <- pages + 1 }
  if (!is.null(opt$gsea_dir)) {
    for (k in seq_len(nrow(contrasts))) {
      ct <- contrasts[k, ]
      for (f in sort(Sys.glob(file.path(opt$gsea_dir, sprintf("GSEA_%s_vs_%s_*.csv", ct$treatment, ct$control))))) {
        print(gsea_plot(ct, f)); pages <- pages + 1
      }
    }
  }
  if (!is.null(opt$gsva_dir)) {
    for (k in seq_len(nrow(contrasts))) {
      ct <- contrasts[k, ]
      for (f in sort(Sys.glob(file.path(opt$gsva_dir, sprintf("GSVA_test_%s_vs_%s_*.csv", ct$treatment, ct$control))))) {
        print(gsva_plot(ct, f)); pages <- pages + 1
      }
    }
  }
  message(sprintf("wrote %s: %d page(s)", opt$out, pages))
}

tryCatch(main(), error = function(e) {
  cat("ERROR: ", conditionMessage(e), "\n", file = stderr())
  quit(status = 1)
})
