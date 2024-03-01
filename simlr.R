library(SIMLR)
library(SingleCellExperiment)


main <- function(filename_dataset){
  dataset = paste0(strsplit(filename_dataset, "_")[[1]][1], "_", strsplit(filename_dataset, "_")[[1]][2])
  folder = paste0("./", dataset, "/simlr")
  if (!file.exists(folder)) {
    dir.create(folder, recursive = TRUE)
  }
  sce<-readRDS(paste0("./data/sce/", filename_dataset, ".RDS"))
  start<-Sys.time()
  
  number_clusters<-SIMLR_Estimate_Number_of_Clusters(X = as.matrix(logcounts(sce)), NUMC = 2:10, cores.ratio = 1)
  number_cluster_opt <- (min(which.min(number_clusters$K1), which.min(number_clusters$K2))+1)
  cluster_output <- SIMLR_Large_Scale(X = as.matrix(counts(sce)), c = number_cluster_opt, k = 10, kk = 100)
  

  png(paste0(folder, paste0("/pca_clust_simlr_",filename_dataset,".png")), width = 800, height = 600)  # Replace with your desired filename
  
  # Create the plot
  plot(cluster_output$ydata, 
       col = cluster_output$y$cluster, 
       xlab = "SIMLR component 1", 
       ylab = "SIMLR component 2", 
        
       main = "SIMILR 2D visualization for BuettnerFlorian")
  
  # Render the plot and close the device (optional for some libraries)
  dev.off()

  end<- Sys.time()
  time = end - start
  
  final <- data.frame(row.names(colData(sce)),cluster_output$y$cluster)
  colnames(final) <- c("cell", "simlr")
  write.csv(final, file=paste0(folder, "/simlr.csv"), quote=FALSE, row.names = FALSE, col.names=TRUE)
  cat(paste0(time, " ", dataset, " simlr"), file="time.txt", append=TRUE, sep = "\n")
  
}