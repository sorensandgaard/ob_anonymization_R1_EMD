#!/usr/bin/env Rscript
library("Seurat")
library("tidyverse")

args = commandArgs(trailingOnly=TRUE)
case_pos <- args[1]
ctrl_pos <- args[2]
outdir <- args[3]

case_obj <- readRDS(case_pos)
ctrl_obj <- readRDS(ctrl_pos)

# Do something - write earth movers distance
