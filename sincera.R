library(SINCERA)
library(SingleCellExperiment)

main<-function(filename_dataset){
  dataset = paste0(strsplit(filename_dataset, "_")[[1]][1], "_", strsplit(filename_dataset, "_")[[1]][2])
  folder = paste0("./", dataset, "/sincera")
  if (!file.exists(folder)) {
    dir.create(folder, recursive = TRUE)
  }
  sce<-readRDS(paste0("./data/sce/", filename_dataset, ".RDS"))
  start<-Sys.time()
  # The analysis starts with running the construct function to create an R S4 object, which will hold all the data and analysis results.
  # The function takes expression matrix and cell sample information as input
  sc <- construct(exprmatrix=data.frame(counts(sce)), samplevector = colnames(counts(sce)))
  
  # Identify and remove low quality cells.
  # The key parameters of running this function include: “min.expression”,
  # which specifies the minimum expression value for a gene to be considered as expressed,
  # and “min.genes”, which specifies the lower bound of the number of expressed genes in a cell.
  sc <- filterLowQualityCells(sc, min.expression=1, min.genes=0, do.plot = F)
  
  # Filter out non-expressed genes
  sc <- prefilterGenes(sc, pergroup=FALSE, min.expression=1, min.cells=1, min.samples=1)
  
  # Perform z-score scaling
  sc <- normalization.zscore(sc, pergroup=FALSE)
  
  # do PCA using all genes
  sc <- doPCA(sc, genes=NULL, use.fast = T)
  
  
  
  #In the first iteration, we selected genes for clustering using the following criteria
  # genes expressed in at least 10 cells in at least 6 samples
  min.samples <- 6
  obj <- prefilterGenes(sc, pergroup=FALSE, min.expression=1, min.cells=10, min.samples=min.samples)
  # genes with at least 0.7 specificity in at least 6 samples
  obj <- cluster.geneSelection(obj, method="specificity", pergroup=TRUE, min.samples=min.samples, specifity.thresh=0.7)
  
  # set the selected genes for clustering
  sc <- setGenesForClustering(sc, value=getGenesForClustering(obj))
  
  sc <- cluster.assignment(sc)
  
  sc_pca<-plotRDS(sc, feature.type="pca", pt.size=1, do.label=F, label.size = 0)
  ggsave(paste0(folder, paste0("/pca_clust_sincera_",filename_dataset,".png")), sc_pca, width = 6, height = 4)  # Replace with your desired filename
  
  end<- Sys.time()
  time<-end-start
  
  fin <- data.frame(cell=colnames(sce), sincera = getCellMeta(sc, "CLUSTER"))
  
  write.csv(fin, file=paste0(folder, "/sincera.csv"), quote=FALSE, row.names = FALSE, col.names=TRUE)
  cat(paste0(time, " ", dataset, " sincera"), file="time.txt", append=TRUE, sep = "\n")
  
  
}