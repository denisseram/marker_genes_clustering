# Apply SC3 #
library(SingleCellExperiment)
library(SC3)
library(scater)
library(igraph)
suppressPackageStartupMessages(library(SummarizedExperiment))
suppressPackageStartupMessages(library(rngtools))


pre_clean <- function(sce) {
  ave.counts <- rowMeans(counts(sce))
  keep <- ave.counts >= 1
  sce_clean <- sce[keep,]
  sce_clean <- addPerCellQC(sce_clean)
  sce_clean <- scater::logNormCounts(sce_clean)
  return(sce_clean)
}

main <- function(filename_dataset) {
  dataset = paste0(strsplit(filename_dataset, "_")[[1]][1], "_", strsplit(filename_dataset, "_")[[1]][2])
  folder = paste0("./", dataset, "/sc3")
  if (!file.exists(folder)) {
    dir.create(folder, recursive = TRUE)
  }
  sce<-readRDS(paste0("./data/sce/", filename_dataset, ".RDS"))

  
  start<-Sys.time()
  sce1<-pre_clean(sce)
  # define feature names in feature_symbol column
  rowData(sce1)$feature_symbol <- rownames(sce1)
  # remove features with duplicated names
  sce1 <- sce1[!duplicated(rowData(sce1)$feature_symbol), ]
  
  sce1 <- runPCA(sce1)
  pca<-plotPCA(sce1, colour_by = "cell_type1")
  ggsave(filename = paste0(folder, "/pca.png"), plot = pca)
  
  
  sce1 <- sc3_estimate_k(sce1)
  k_est <- sce1@metadata$sc3$k_estimation
  sce1 <- sc3(sce1, ks = k_est, biology = TRUE)
  
  sce1 <- runPCA(sce1)
  pca1<-plotPCA(
    sce1, 
    colour_by = paste0("sc3_",k_est,"_clusters"), 
  )
  ggsave(filename = paste0(folder, paste0("/pca_clust_",filename_dataset,".png")), plot = pca1)
  
  end <- Sys.time()
  eval(parse(text=paste0("colData(sce1)$SC3 <- colData(sce1)$sc3_", k_est, "_clusters")))
  final <-data.frame(row.names(colData(sce1)), colData(sce1)$SC3)
  colnames(final) <- c("cell", "sc3")
  
  time <- end - start
  write.csv(final, file=paste0(folder, "/sc3.csv"), quote=FALSE, row.names = FALSE, col.names=TRUE)
  cat(paste0(time, " ", dataset, " sc3"), file="time.txt", append=TRUE, sep = "\n")
  
 
  
}