#######################
#### Functions #######
#####################
library(ggplot2)

calc_ARI <- function(M, N) {
  
  library("mclust")
  
  if (length(M) != length(N)) {
    print("Error. The length of vectors should be the same.")
    ARI <- NA
  } else {
    ARI <- adjustedRandIndex(M,N)
  }
  return(ARI)
}

calc_AMI <- function(M,N) {
  
  library("aricode")
  
  if (length(M) != length(N)) {
    print("Error. The length of vectors should be the same.")
    AMI <- NA
  } else {
    AMI <- AMI(M,N)
  }
  return(AMI)
}

calc_vi <- function(M, N) {
  library("mcclust")
  
  if (length(M) != length(N)) {
    print("Error. The length of vectors should be the same.")
    VI.measure <- NA
  } else {
    VI.measure <- vi.dist(M, N)
  }
  return(VI.measure)
}
############################
#### Setting ##############
##########################
all_datasets<-c("baron-mouse","baron-mouse_markergenes")
all_clusterings <- c( "seurat","sc3","cidr","sincera","simlr","tscan","raceid")
#all_clusterings <- c("cidr")


for (i in all_datasets) {
  print(i)
  for (j in all_clusterings) {
    source(paste0("./", j, ".R"))
    main(i)
  }
}

 ############################
#### Analysis #############
##########################



ARI_matrix<-matrix("", nrow=length(all_clusterings), ncol=length(all_datasets), dimnames=list(all_clusterings, all_datasets))
AMI_matrix<-matrix("", nrow=length(all_clusterings), ncol=length(all_datasets), dimnames=list(all_clusterings, all_datasets))
VI_matrix<-matrix("", nrow=length(all_clusterings), ncol=length(all_datasets), dimnames=list(all_clusterings, all_datasets))
number_clusters<-matrix("", nrow=length(all_clusterings), ncol=length(all_datasets), dimnames=list(all_clusterings, all_datasets))

for (i in all_clusterings) {
  print(i)
  for (j in all_datasets) {
    dataset = paste0(strsplit(j, "_")[[1]][1], "_", strsplit(j, "_")[[1]][2])
    resultados<-read.csv(paste0("./",dataset,"/",i, "/",i, ".csv"))
    
    #leemos el dataset y sus anotaciones y lo guardamos como cell_type_annotations
    sce<-readRDS(paste0("./data/sce/", j, ".RDS"))
    cell_annotations <- colData(sce)
    cell_type_annotations <- cell_annotations$cell_type1
    
    #leemos la asignación de clustering realizado para ese dataset en el algoritmo i
    cluster_assignation <- resultados[[i]]
    
    # Obtener la fila para el algoritmo i
    row_j <- which(dimnames(ARI_matrix)[[1]] == i)
    
    # Obtener la columna para el dataset j
    column_i <- which(colnames(ARI_matrix) == j)
    
    #calculamos ARI
    ari<-calc_ARI(cell_type_annotations,cluster_assignation)
    ARI_matrix[row_j, column_i] <- ari
    
    #calculamos AMI
    ami<-calc_AMI(cell_type_annotations,cluster_assignation)
    AMI_matrix[row_j, column_i] <- ami
    #calculamos VI
    vi<-calc_vi(cell_type_annotations,cluster_assignation)
    VI_matrix[row_j, column_i] <- vi
    #calcular numero de clusters
    cluster_assignation1<-as.data.frame(cluster_assignation)
    num_clu<-length(unique(cluster_assignation1$cluster_assignation))
    number_clusters[row_j, column_i] <- num_clu
  }
}
write.csv(ARI_matrix, file = "ARI_matrix.csv", row.names = FALSE)
write.csv(AMI_matrix, file = "AMI_matrix.csv", row.names = FALSE)
write.csv(VI_matrix, file = "VI_matrix.csv", row.names = FALSE)

