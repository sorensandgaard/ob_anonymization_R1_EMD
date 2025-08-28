#!/usr/bin/env Rscript
library("tidyverse")
library("Seurat")
library("Matrix")
library("EMDomics")
library("jsonlite")
library("EMDomics")

rm(list=ls())

args = commandArgs(trailingOnly=TRUE)
case_pos <- args[1]
ctrl_pos <- args[2]
outdir <- args[3]


# setwd("/home/projects/dtu_00062/people/sorsan/ob_anonymization_dataloss")
# case_pos <- "out/data/D1/default/init_align/P1/param_0/second_align/M1/param_0/D1_case.rds"
# ctrl_pos <- "out/data/D1/default/init_align/P1/param_0/second_align/M1/param_0/D1_ctrl.rds"
# outdir <- "out/data/D1/default/init_align/P1/param_0/second_align/M1/param_0/metrics/m1/default/D1.somefile.txt"

case_obj <- readRDS(case_pos)
ctrl_obj <- readRDS(ctrl_pos)

# Normalize
print("Normalizing")
case_obj <- NormalizeData(case_obj)
ctrl_obj <- NormalizeData(ctrl_obj)

# Convert Seurat objects to expression matrices
print("Extracting expression matrices")
expr_orig <- case_obj[["RNA"]]$data
expr_anon <- ctrl_obj[["RNA"]]$data

# Check genelists are identical and save genelist and genecount
print("Check whether genelists are identical")
genelist_orig <- rownames(expr_orig)
genelist_anon <- rownames(expr_anon)

try(if(!identical(genelist_orig,genelist_anon)) stop("Genelists are not identical"))
genelist <- genelist_orig
rm(genelist_orig,genelist_anon)

# Keep only genes that are non-zero in at least one matrix
print("Remove genes that have zero expression in both matrices")
non_zero_anon <- rowSums(expr_anon > 0) > 0
non_zero_orig <- rowSums(expr_orig > 0) > 0
common_genes <- non_zero_anon & non_zero_orig

expr_anon <- as.matrix(expr_anon[common_genes, ])
expr_orig <- as.matrix(expr_orig[common_genes, ])

# Combine the expression matrices
print("Combine matrices to one")
combined_matrix <- cbind(expr_anon, expr_orig)
colnames(combined_matrix) <- paste("sample", 1:ncol(combined_matrix), sep="")

# Create labels vector
print("Creating label vectors")
labels <- c(rep("Group1", ncol(expr_anon)), rep("Group2", ncol(expr_orig)))
names(labels) <- colnames(combined_matrix)

# Calculate EMD
print("Calculating EMD")
emd_results <- calculate_emd(combined_matrix, labels,nperm = 2,parallel = F)

emd_results <- data.frame(emd_results$emd)
poor_genes <- filter(emd_results,emd > 1.0)

print("Creating JSON object")
json_obj <- list(
  module = "R1_EMD",
  timestamp = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ"),
  metrics = list(
    earth_movers_distance = list(
      gene_wise = list(
        mean = mean(emd_results$emd),
        median = median(emd_results$emd),
        sd = sd(emd_results$emd),
        min = min(emd_results$emd),
        max = max(emd_results$emd),
        q25 = quantile(emd_results$emd, 0.25),
        q75 = quantile(emd_results$emd, 0.75)
      ),
      poorly_matched_genes = list(
        count = sum(emd_results$emd > 1.0),
        threshold = 1.0,
        gene_ids = rownames(poor_genes),
        emd_values = poor_genes$emd
      ),
      all_results = list(
        gene_ids = rownames(emd_results),
        emd_values = emd_results$emd
      )
    )
  )
)

print("Writing JSON")
write_json(json_obj, "/home/projects/dtu_00062/people/sorsan/ob_anonymization_dataloss/test.json", pretty = TRUE, auto_unbox = TRUE)
write_json(json_obj, outdir, pretty = TRUE, auto_unbox = TRUE)
