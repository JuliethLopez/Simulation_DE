#Context: Assuming a healthy individual, we want to find differentially expressed genes for each cell type
#Objective: Compare differential expression methods performance for single cell
# You can peak some DE methods exploring this article:
# Article - (2022) Recommendations of scRNA-seq Differential Gene Expression Analysis Based on Comprehensive Benchmarking. https://doi.org/10.3390/life12060850
#Tips: Some metrics to compare DEGS in methods are accuracy, precision, recall, F1 score
# Take only overexpressed genes in the cluster of interest
# Filter DEGS by: fold change and pvalue
# See some differences in a fraction of cells expressing a given gene

################################################################################
# Data and libraries
################################################################################

# Set environment
setwd("C:/Users/Administrator/Documents/Jaguar/Simulation_DE")
outdir <- "C:/Users/Administrator/Documents/Jaguar/Simulation_DE/figures"

# libraries
library("scater") #normalization and pca
library("splatter") #simulation
library("Seurat") #single cell pre-processing and processing
library("SeuratData") #data reference PBMC
library("ggplot2") #plotting
library("tidyverse") #database manipulation

# charge data
data <- LoadData("pbmc3k.SeuratData") #not normalized

round(table(data$seurat_annotations)/sum(table(data$seurat_annotations)),2)

# extract counts
data_counts <- as.matrix(data@assays$RNA$counts)
class(data_counts) #matrix
dim(data_counts) #13714  2700

# Cells frecuency
sampled_cells_indices <- sample(1:ncol(data_counts), size = 1000, replace = FALSE)
#Naive CD4 T Memory CD4 T   CD14+ Mono     B        CD8 T 
#0.26         0.18         0.18         0.13         0.10 
#FCGR3A+ Mono   NK           DC     Platelet 
#0.06         0.06         0.01         0.01 

################################################################################
# simulation scRNA data ONE individual with different types of cells
################################################################################
# Article - (2017) Splatter: simulation of single-cell RNA sequencing data. https://doi.org/10.1186/s13059-017-1305-0
# Tutorial - Splatter: https://bioconductor.statistik.tu-dortmund.de/packages/3.13/bioc/vignettes/splatter/inst/doc/splatter.html

set.seed(1)
params <- splatEstimate(subset_counts)
params
getParams(params, names = c("nGenes","nCells", "mean.rate", "mean.shape", "nGroups"))
#$nGenes
#[1] 13714
#
#$nCells
#[1] 1000
#
#$mean.rate
#[1] 14.3284
#
#$mean.shape
#[1] 0.6715326
#
#$nGroups
#[1] 1

sim.groups <- splatSimulateGroups(batchCells=3000, #number of cells per batch
                                  group.prob = c(0.54, 0.24, 0.13, 0.06, 0.01), #number and proportion groups: T, mono, B, NK, DC
                                  de.prob = c(0.08, 0.06, 0.05, 0.05, 0.01), #probability DEGS per group
                                  de.facLoc = c(0.01, 0.2, 0.1, 0.1, 0.2), #define more or less extreme differences between groups
                                  de.facScale = c(0.2, 0.5, 0.2, 0.5, 0.4), #define more or less extreme differences between groups
                                  de.downProb = 0, #probability of subexpressed genes
                                  sparsify = TRUE, #a lot of zeros
                                  verbose = TRUE) #print additional info

# Explore simulation
sim.groups #class SingleCellExperiment
head(sim.groups$Group,20) #group the cell belongs to
head(sim.groups@rowRanges@elementMetadata$GeneMean,20) #Expression level of a gene in a particular group after applying differential expression factors.
sim.groups@rowRanges@elementMetadata$DEFacGroup1 #The differential expression factor for each gene in a particular group (1 is not differentially expressed)

# Comparison mean expression simulation vs real data
comparison <- compareSCEs(list(Splat = sim.groups, Real = as.SingleCellExperiment(data)))
plot_means_comparison <- comparison$Plots$Means
ggsave(filename = paste0(outdir,"/Mean_expression_simulation_vs_realdata.png"), plot = plot_means_comparison)

# Use scater to calculate logcounts and quick pca
sim.groups <- logNormCounts(sim.groups)
# Plot PCA
sim.groups <- runPCA(sim.groups)
pca_groups <- plotPCA(sim.groups, colour_by = "Group") + labs(title="Five groups with different DE probabilities")
ggsave(filename = paste0(outdir,"/pca_groups.png"), plot = pca_groups)

# save simulation
saveRDS(sim.groups, file = "simulation.rds")

################################################################################
# pre-processing with seurat
################################################################################
# Tutorial: https://satijalab.org/seurat/articles/pbmc3k_tutorial

# read rds
sim.groups <- readRDS("simulation.rds")

# count matrix
matrix_counts <- sim.groups@assays@data@listData[["counts"]] #not normalized
dim(matrix_counts)

# CreateSeuratObject
data <- CreateSeuratObject(counts = matrix_counts)

# normalization
data <- NormalizeData(data, normalization.method = "LogNormalize", scale.factor = 10000)

#Identification of highly variable features (feature selection)
data <- FindVariableFeatures(data, selection.method = "vst", nfeatures = 2000)

# Identify the most highly variable genes
top20 <- head(VariableFeatures(data), 20)
VariableFeaturePlot(object = data) #plot most highly variable genes

# Scaling the data for dimension reduction
data <- ScaleData(data)

#Perform linear dimensional reduction
#PCA
data <- RunPCA(data, features = VariableFeatures(object = data))
VizDimLoadings(data, dims = 1:2, reduction = "pca")
DimPlot(data, reduction = "pca")

# Elbow plot to pick enough principal components
ElbowPlot(data)

#Cluster the cells
data <- FindNeighbors(data, dims = 1:5)
data <- FindClusters(data, resolution = 0.5)

# Run UMAP
data <- RunUMAP(data, dims = 1:10)
DimPlot(data, reduction = "umap", pt.size = 0.5) + NoLegend()

################################################################################
# differential expression
################################################################################
# Tutorial: https://satijalab.org/seurat/articles/de_vignette

## differential expression tests
# head(sim.groups$Group,20) #group the cell belongs to
Idents(data_norm_harmony_bulk) <- "---" #idents are the cell types

#Finding differentially expressed features (cluster biomarkers)
markers <- FindAllMarkers(data, test.use = "wilcox", assay="----") #complete the assay
write.csv(markers, file = paste0("/markers_wilcox.csv"), row.names = TRUE) #save markers

# plot gene counts
VlnPlot(data, features = c("Gene1", "Gene2"), slot = "counts", log = TRUE)

## pseudobulk test Deseq2
# agregate data, only when you have more than one individual
# data_bulk <- AggregateExpression(data, assays="---", return.seurat = T, group.by = c("----")) #group by cell type
markers_bulk <- FindAllMarkers(object = data_bulk, test.use = "DESeq2")
write.csv(markers_bulk, file = paste0("/markers_DESEQ.csv"), row.names = TRUE)

# compare the DE p-values between the single-cell level and the pseudobulk level results
markers_bulk <- read.csv(paste0("/markers_DESEQ_azimuth_L1_l2.csv"))
names(markers_bulk) <- paste0(names(markers_bulk), ".bulk")
names(markers_bulk)[names(markers_bulk)=="gene.bulk"] = "gene"
dim(markers_bulk)

markers_all_azimuth_t_L1 <- read.csv(paste0("/markers_all_ttest_azimuth_L1_l2.csv"))
names(markers_all_azimuth_t_L1) <- paste0(names(markers_all_azimuth_t_L1), ".sc")
names(markers_all_azimuth_t_L1)[names(markers_all_azimuth_t_L1)=="gene.sc"] = "gene"
dim(markers_all_azimuth_t_L1)

merge_dat <- merge(markers_all_azimuth_t_L1, markers_bulk, by = "gene")
merge_dat <- merge_dat[order(merge_dat$p_val.bulk), ]

# Number of genes that are marginally significant in both; marginally significant only in bulk; and marginally significant only in single-cell
common <- merge_dat$gene[which(merge_dat$p_val.bulk < 0.05 &
                                 merge_dat$p_val.sc < 0.05)]
only_sc <- merge_dat$gene[which(merge_dat$p_val.bulk > 0.05 &
                                  merge_dat$p_val.sc < 0.05)]
only_bulk <- merge_dat$gene[which(merge_dat$p_val.bulk < 0.05 &
                                    merge_dat$p_val.sc > 0.05)]
