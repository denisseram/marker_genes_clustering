# Load required libraries
library(data.table)

# Read the data files
sce<-readRDS("./data/sce/baron-mouse.rds")
marker_mouse<-read.csv("./marker_genes/markers_pancreas_mouse.csv")
# Filter genes by matching row names to marker symbols
rownames(sce)<-toupper(rownames(sce))
rownames(sce)
filtered_sce <- sce[rownames(sce) %in% allmarkers$V1, ]

# Print dimensions of the filtered object
dim(filtered_sce)


saveRDS(filtered_sce,  "./data/sce/baron-mouse_markergenes.RDS")

###WITH HUMANDATASET
# Read the data files
sce<-readRDS("./data/sce/baron-human.rds")
marker_mouse<-read.csv("./marker_genes/markers_pancreas_mouse.csv")
# Filter genes by matching row names to marker symbols
rownames(sce)<-toupper(rownames(sce))
rownames(sce)
filtered_sce <- sce[rownames(sce) %in% allmarkers$V1, ]

# Print dimensions of the filtered object
dim(filtered_sce)


saveRDS(filtered_sce,  "./data/sce/baron-mouse_markergenes.RDS")




datos <- read.table("contain_baron.txt")
allmarkers <- read.table("all_markers.txt")
head(allmarkers)


si<-rownames(datos)
uppercase_vector <- toupper(si)

hola<-uppercase_vector %in% allmarkers$V1
count(hola)
