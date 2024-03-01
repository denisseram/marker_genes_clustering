# CIDR #
library(cidr)
library(SingleCellExperiment)


main <- function(filename_dataset){

  dataset = paste0(strsplit(filename_dataset, "_")[[1]][1], "_", strsplit(filename_dataset, "_")[[1]][2])
  folder = paste0("./", dataset, "/cidr")
  if (!file.exists(folder)) {
    dir.create(folder, recursive = TRUE)
  }
  sce<-readRDS(paste0("./data/sce/", filename_dataset, ".RDS"))
  start<- Sys.time()
  
  # Make an object
  sce_obj <- scDataConstructor(assay(sce), tagType = "raw")
  
  # Identify the dropouts candidates
  sce_obj <- determineDropoutCandidates(sce_obj, min1 = 3, min2 = 8, N = 2000, alpha = 0.1, fast = TRUE, zerosOnly = FALSE,
                                        bw_adjust = 1)
  
  # wThreshold function is applied to set weighted gene expression thresholds, helping to improve dropout detection
  sce_obj <- wThreshold(sce_obj)
  # Calculate the dissimilarity between cells
  sce_obj <- scDissim(sce_obj)
  # A PCA is made
  sce_obj <- scPCA(sce_obj)
  # Calculate the optimal number of PC on the next analysis
  sce_obj <- nPC(sce_obj)
  # Number of clusters
  nCluster(sce_obj)
  
  sce_obj <- scCluster(sce_obj)
  plot.new()
  #legend("center", legend=types, col=scols, pch=15)
  plot(sce_obj@PC[,c(1,2)], col=sce_obj@clusters,pch=sce_obj@clusters, 
       main="CIDR", xlab="PC1", ylab="PC2")
  dev.off()
  gg <- ggplot(data.frame(PC1 = sce_obj@PC[, 1], PC2 = sce_obj@PC[, 2], 
                          clusters = as.factor(sce_obj@clusters)), 
               aes(x = PC1, y = PC2, color = clusters)) +
    geom_point() +
    labs(title = "CIDR", x = "PC1", y = "PC2")
  
  ggsave(paste0(folder, paste0("/pca_clust_cidr_",filename_dataset,".png")), gg, width = 6, height = 4, units = "in")
  

  end<-Sys.time()
  time<-end-start
  
  final<-data.frame(row.names(colData(sce)), sce_obj@clusters)
  colnames(final) <- c("cell", "cidr")
  
  write.csv(final, file=paste0(folder, "/cidr.csv"), quote=FALSE, row.names = FALSE, col.names=TRUE)
  cat(paste0(time, " ", dataset, " cidr"), file="time.txt", append=TRUE, sep = "\n")
  
  
}