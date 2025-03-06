library(tidyverse)
library(Seurat)
library(cyCombine)
library(parallel)
library(Matrix)

# Load and preprocess data
set.seed(1)
expr_orig <- ReadMtx(mtx="/home/projects/dtu_00062/people/sorsan/test3/orig_processed_hg38-2024-A/outs/raw_feature_bc_matrix/matrix.mtx.gz",
                     cells = "/home/projects/dtu_00062/people/sorsan/test3/orig_processed_hg38-2024-A/outs/raw_feature_bc_matrix/barcodes.tsv.gz",
                     features = "/home/projects/dtu_00062/people/sorsan/test3/orig_processed_hg38-2024-A/outs/raw_feature_bc_matrix/features.tsv.gz")
expr_anon <- ReadMtx(mtx="/home/projects/dtu_00062/people/sorsan/test3/anon_processed_hg38-2024-A/outs/raw_feature_bc_matrix/matrix.mtx.gz",
                     cells = "/home/projects/dtu_00062/people/sorsan/test3/anon_processed_hg38-2024-A/outs/raw_feature_bc_matrix/barcodes.tsv.gz",
                     features = "/home/projects/dtu_00062/people/sorsan/test3/anon_processed_hg38-2024-A/outs/raw_feature_bc_matrix/features.tsv.gz")

# Load metadata
# Check genelists are identical and save the
genelist_orig <- rownames(expr_orig)
genelist_anon <- rownames(expr_anon)

try(if(!identical(genelist_orig,genelist_anon)) stop("Genelists are not identical"))
genelist <- genelist_orig
rm(genelist_orig,genelist_anon)

genecount <- dim(expr_orig)[1]

cellcount_orig <- dim(expr_orig)[2]
cellcount_anon <- dim(expr_anon)[2]

# Initialize clusters
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

# Run parallel script
emds <- parLapply(cl = cl,X = 1:genecount, function(i){
  tmp_1 <- expr_orig[i,]
  tmp_2 <- expr_anon[i,]
  tmp_df <- data.frame(
    exp = c(tmp_1,tmp_2),
    batch = rep(c("orig","anon"),times=c(cellcount_orig,cellcount_anon)),
    label = 1
  )
  
  tmp <- compute_emd(tmp_df,markers = "exp", cell_col = 'label')
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

# Write EMD dataframe to file
write.table(emds, "/home/projects/dtu_00062/people/sorsan/test3/EMDs_2020_v_2024.txt", sep="\t")





