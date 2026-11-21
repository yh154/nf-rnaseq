#!/usr/bin/env Rscript
# Usage: check_dds.R DDS_RDS DESIGN_FORMULA N_SAMPLES   (run inside the deseq2 container)
suppressPackageStartupMessages(library(DESeq2))
a <- commandArgs(trailingOnly = TRUE)
if (length(a) != 3) { cat("Usage: check_dds.R DDS_RDS DESIGN_FORMULA N_SAMPLES\n"); quit(status = 2) }
dds <- readRDS(a[1])
stopifnot(is(dds, "DESeqDataSet"))
stopifnot(identical(deparse(design(dds)), deparse(as.formula(a[2]))))
stopifnot(ncol(dds) == as.integer(a[3]))
stopifnot(!is.null(sizeFactors(dds)), !is.null(dispersions(dds)))
stopifnot(all(counts(dds) == round(counts(dds))))
cat("ok  dds.rds loads, design =", deparse(design(dds)), ",", ncol(dds), "samples,", nrow(dds), "genes fitted\n")
