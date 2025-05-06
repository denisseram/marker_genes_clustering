# Single-Cell RNA-seq Clustering with Marker Gene Integration

This repository contains the implementation of a comprehensive methodology to evaluate the impact of marker genes in improving clustering performance for single-cell RNA sequencing (scRNA-seq) data analysis.
## Project Overview

scRNA-seq has transformed our understanding of cellular heterogeneity. A critical step in analyzing scRNA-seq data is identifying distinct cell types through unsupervised clustering. However, the sparse and high-dimensional nature of scRNA-seq data poses significant challenges for accurate clustering.
This project tests the hypothesis that selecting an appropriate subset of genes (marker genes) can enhance the clustering process by reducing noise, ultimately improving the identification of cell types.

## Methodology

![Image](https://github.com/user-attachments/assets/3399b9de-12fe-4f39-b247-5510a15637d7)

### Data Preprocessing:

Download and preprocess scRNA-seq dataset
Obtain marker genes from Cell Marker and PanglaoDB databases


### Feature Selection:

Filter expression matrix to retain only marker genes
Create a second dataset with these filtered features


### Clustering:

Apply seven clustering algorithms to both full and filtered datasets


### Evaluation:

Compare clustering performance using multiple metrics:

Adjusted Rand Index (ARI)
Adjusted Mutual Information (AMI)
Variation of Information (VI)
Number of clusters

## Clustering Algorithms
Seven different clustering algorithms are implemented and evaluated:

Seurat - A widely used R package for scRNA-seq analysis
SC3 - Single-Cell Consensus Clustering
CIDR - Clustering through Imputation and Dimensionality Reduction
SINCERA - SINgle CEll RNA-seq profiling Analysis
SIMLR - Single-cell Interpretation via Multi-kernel LeaRning
TSCAN - Tools for Single-Cell ANalysis
RaceID - Robust analysis of single-cell RNA-seq data

## Repository Structure
```
.
├── data/
│   └── sce/               # Single-cell experiment objects (RDS files)
├── marker_genes/          # CSV files containing marker genes for human and mouse pancreas
├── baron-mouse/           # Results for baron-mouse dataset
│   ├── cidr/
│   ├── raceid/
│   ├── sc3/
│   ├── seurat/
│   ├── simlr/
│   ├── sincera/
│   └── tscan/
├── baron-mouse_markergenes/ # Results for filtered baron-mouse dataset
│   ├── ...                  # Same structure as above
├── *.R                     # R scripts for clustering and analysis
├── *.png                   # Result visualizations
└── time.txt                # Runtime benchmarks
```


## Usage
### Prerequisites


Install required R packages:
```
Rinstall.packages(c("data.table", "cidr", "SingleCellExperiment", "RaceID", "SC3", "scater",
                   "SIMLR", "SINCERA", "dplyr", "ggplot2", "mclust", "aricode", "mcclust"))
```
### Running the Analysis

Filter genes using marker genes:
```
RRscript Filtering_by_markers.R
```
Run all clustering algorithms:

RRscript principal.R

Generate comparison plots:

RRscript save_samesize.R
