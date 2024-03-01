library(RaceID)
library(SingleCellExperiment)

main <- function(filename_dataset){

  dataset = paste0(strsplit(filename_dataset, "_")[[1]][1], "_", strsplit(filename_dataset, "_")[[1]][2])
  folder = paste0("./", dataset, "/raceid")
  if (!file.exists(folder)) {
    dir.create(folder, recursive = TRUE)
  }
  sce<-readRDS(paste0("./data/sce/", filename_dataset, ".RDS"))
  start<-Sys.time()
  
  sc <- SCseq(counts(sce))
  # Filters cells with less than 2000 total counts
  sc <- filterdata(sc,mintotal=1)
  # Retrieves filtered expression data from the object
  fdata <- getfdata(sc)
  # Computes a pairwise distance matrix using Pearson correlation
  sc <- compdist(sc,metric="pearson")
  # Performs initial clustering of cells
  sc <- clustexp(sc)  
  
  sc <- findoutliers(sc)
  
  sc <- compumap(sc)  # Computes UMAP
  
  plotmap(sc, um=TRUE)  # Plots UMAP
  dev.copy2pdf(file = (paste0(folder, paste0("/pca_clust_raceid_",filename_dataset,".pdf"))))
  
  dev.off()
  
  
  final <- data.frame(colnames(sce), sc@cpart)
  final$RaceID <- sc@cpart
  num_clusters<-length(unique(final$RaceID3))
  colnames(final) <- c("cell", "raceid")
  
  end <- Sys.time()
  time<-end-start
  
  write.csv(final, file=paste0(folder, "/raceid.csv"), quote=FALSE, row.names = FALSE, col.names=TRUE)
  cat(paste0(time, " ", dataset, " raceid"), file="time.txt", append=TRUE, sep = "\n")
  
  
}