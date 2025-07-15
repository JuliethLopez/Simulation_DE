##### 
rm(list=ls(all=TRUE))
par(mfrow = c(1, 1)) 
plot.new()

## Set environment - Change for your own directory
basedir <- "C:/Users/Administrator/Documents/Jaguar/Simulation_DE/Mast"
setwd("C:/Users/Administrator/Documents/Jaguar/Simulation_DE")
outdir <- "C:/Users/Administrator/Documents/Jaguar/Simulation_DE/Mast/figures"
source("C:/Users/Administrator/Documents/Jaguar/Simulation_DE/Functions.R")

## Libraries
library(Seurat)
library(dplyr)
library(ggplot2)
library(reshape2)
#library(data.table)
library(openxlsx)
library(MAST)

################################################################################
## Data
################################################################################

# Data from simulation
# Here is called data pre-processed in ALL_V1
data <- readRDS("Simulaciones/simulation_2.rds")
data_1_cl <- as.Seurat(data, counts = "counts", data = NULL)
data_1_cl #15000 features across 4000 samples

# change name off assay from "originalexp" to "RNA"
names(data_1_cl@assays) = "RNA"
data_1_cl@active.assay = "RNA"
# change "Groups" name to "Cluster"
colnames(data_1_cl@meta.data)[6] = "Cluster"

# inspect
#head(data_1_cl@assays$RNA$data)
#head(data_1_cl@assays$RNA$counts)
#head(data_1_cl@meta.data)
#head(data_1_cl@assays$RNA@meta.features)

# Preprocessing
data_1_cl <- NormalizeData(data_1_cl, normalization.method = "LogNormalize", scale.factor = 10000)
data_1_cl <- FindVariableFeatures(data_1_cl, selection.method = "vst", nfeatures = 2000)
data_1_cl <- ScaleData(data_1_cl)
data_1_cl <- RunPCA(data_1_cl, features = VariableFeatures(object = data_1_cl))
data_1_cl <- RunUMAP(data_1_cl, dims = 1:10)

################################################################################
# Expected DEGS per group and level
expected_degs_levels_df <- read.csv("DEGs_expected.csv")
all_genes_df <- read.csv("all_genes.csv")
all_genes_df <- all_genes_df %>% select(c(Gene, cluster))

################################################################################
## Differential expression analysis
################################################################################
## mast
Idents(data_1_cl) <- "Cluster"
DE_mast <- DE_analysis_seurat(data_1_cl, assay = "RNA", slot="data", test_type = "MAST")
head(DE_mast)
dim(DE_mast) #1826    7

# save results
write.csv(DE_mast, file = paste0(basedir,"/DE_mast.csv"), row.names = TRUE) #save markers

## Filter of DE results to define the markers
DE_mast = read.csv(paste0(basedir,"/DE_mast.csv"))
markers_mast = DE_classes(DE_mast, min_pct_difference = 0.1, max_p_val_adj = 0.05, type="predicted")

# Data DE results/metrics from simulation
# list for save metrics
results_list <- list()

# Expected DEGS per group and level
results_list[["mast"]] <- DE_metrics(expected_degs_levels_df,markers_mast, all_genes_df, test="mast")
results_list[["mast"]]$Metrics

## Heatmap summary per test
Metrics_heatmap(results_list, outdir, text="_simulation")

## Save the DEGs (given the filters aplpied), the DEGs per cluster and the metrics to an Excel file
wb <- createWorkbook()

for (test in names(results_list)) {
  DEGs <- results_list[[test]]$DEGs
  metrics <- results_list[[test]]$Metrics
  
  # DEGs
  addWorksheet(wb, paste(test, "DEGs"))
  writeData(wb, sheet = paste(test, "DEGs"), DEGs)
  
  # metrics
  addWorksheet(wb, paste(test, "Metrics"))
  writeData(wb, sheet = paste(test, "Metrics"), metrics)
  
  # predicted data by group in all test
  DEGs_per_group <- DEGs %>%
    group_by(cluster) %>% #logfc_level
    summarise(count = n())
  DEGs_per_group <- rbind(DEGs_per_group, data.frame(cluster="all", count = sum(DEGs_per_group$count)))
  
  addWorksheet(wb, paste(test, "DEGs_group"))
  writeData(wb, sheet = paste(test, "DEGs_group"), DEGs_per_group)
}

saveWorkbook(wb, file = paste0(basedir,"/DEGs_and_Metrics_MAST.xlsx"), overwrite = TRUE)
