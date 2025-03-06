#!/usr/bin/env Rscript
library("tidyverse")
library("Seurat")
library("cyCombine")
library("parallel")
library("Matrix")


args = commandArgs(trailingOnly=TRUE)
case_pos <- args[1]
ctrl_pos <- args[2]
outdir <- args[3]

# case_pos = out/data/D1/default/init_align/P1/param_0/second_align/M1/param_0/D1_case.rds
# ctrl_pos = out/data/D1/default/init_align/P1/param_0/second_align/M1/param_0/D1_ctrl.rds
# outdir = out/data/D1/default/init_align/P1/param_0/second_align/M1/param_0/metrics/m1/default/D1.somefile.txt

case_obj <- readRDS(case_pos)
ctrl_obj <- readRDS(ctrl_pos)

# Do something - write earth movers distance

# Convert Seurat objects to expression matrices
expr_orig <- case_obj[["RNA"]]$counts
expr_anon <- case_obj[["RNA"]]$counts


# Check genelists are identical and save genelist and genecount
genelist_orig <- rownames(expr_orig)
genelist_anon <- rownames(expr_anon)

try(if(!identical(genelist_orig,genelist_anon)) stop("Genelists are not identical"))
genelist <- genelist_orig
rm(genelist_orig,genelist_anon)

genecount <- length(genelist)

# Get cellcounts for each matrix
cellcount_orig <- dim(expr_orig)[2]
cellcount_anon <- dim(expr_anon)[2]

# Create clusters for parallel processing
cl <- makeCluster(detectCores())
print("Number of clusters")
print(length(cl))
clusterExport(cl, c("expr_orig","expr_anon","genelist","genecount","cellcount_orig","cellcount_anon","compute_emd"))

# Ensure proper initialization of clusters
tmp <- parLapply(cl = cl, X = 1:length(cl),function(i){
  try(expr_orig[1,1])
  TRUE
})
rm(tmp)

# Run parallel script for EMD calculation
emds <- parLapply(cl = cl,X = 1:genecount, function(i){
  # Gene expression for i'th gene
  tmp_1 <- expr_orig[i,]
  tmp_2 <- expr_anon[i,]
  
  # Save expression as batches
  tmp_df <- data.frame(
    exp = c(tmp_1,tmp_2),
    batch = rep(c("orig","anon"),times=c(cellcount_orig,cellcount_anon)),
    label = 1
  )
  
  # Compute EMD
  tmp <- compute_emd(tmp_df,markers = "exp", cell_col = 'label')
  
  # Save EMD to dataframe
  tmp_emd <- data.frame(tmp)[1,2]
  data.frame(gene = genelist[i],
             emd = tmp_emd,
             read_count_orig = sum(tmp_1),
             read_count_anon = sum(tmp_2))
})

# Stop parallelization
stopCluster(cl = cl)
rm(cl)

# Create dataframe of EMDs
emds <- Reduce("rbind", emds)
rownames(emds) <- emds$gene
emds <- emds %>% 
  mutate(difference = read_count_orig - read_count_anon)

### Write EMD dataframe to file ###
write.table(emds, file = outdir, sep="\t")
