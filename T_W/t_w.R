##### 
rm(list=ls(all=TRUE))
par(mfrow = c(1, 1)) 
plot.new()

## Set environment - Change for your own directory
basedir <- "C:/Users/Administrator/Documents/Jaguar/Simulation_DE/T_W"
basedir2 <- "C:/Users/Administrator/Documents/Jaguar/Simulation_DE/Wilcoxon_t_Test"
setwd("C:/Users/Administrator/Documents/Jaguar/Simulation_DE")
outdir <- "C:/Users/Administrator/Documents/Jaguar/Simulation_DE/T_W/figures"
source("C:/Users/Administrator/Documents/Jaguar/Simulation_DE/Functions.R")

## Libraries
library(Seurat)
library(dplyr)
library(ggplot2)
library(reshape2)
#library(data.table)
library(openxlsx)

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
## T
Idents(data_1_cl) <- "Cluster"
DE_t <- DE_analysis_seurat(data_1_cl, assay = "RNA", slot="data", test_type = "t")
head(DE_t)
dim(DE_t) #1679    7

DE_t[which(DE_t$gene == "Gene3"),]

## W
DE_w <- DE_analysis_seurat(data_1_cl, assay = "RNA", slot="data", test_type = "wilcox")
head(DE_w)
dim(DE_w) #2714    7

DE_w[which(DE_w$gene == "Gene3"),]

# save results
write.csv(DE_w, file = paste0(basedir,"/DE_W.csv"), row.names = TRUE) #save markers
write.csv(DE_t, file = paste0(basedir,"/DE_t.csv"), row.names = TRUE) #save markers

## Filter of DE results to define the markers
DE_W = read.csv(paste0(basedir,"/DE_w.csv"))
DE_t = read.csv(paste0(basedir,"/DE_t.csv"))

# Cargar los archivos de resultados
DE_t_test <- read.csv(paste0(basedir2,"/Results/DE_results_t.csv"))
DE_Wilcox <- read.csv(paste0(basedir2,"/Results/DE_results_wilcox.csv"))

# Why they are not the same?
DE_t_test$p_val_adj == DE_t$p_val_adj

markers_t <- DE_classes(DE_W, min_pct_difference = 0.1, max_p_val_adj = 0.05, type = "predicted")
head(markers_t)
class(markers_t$cluster)

markers_w <- DE_classes(DE_t, min_pct_difference = 0.1, max_p_val_adj = 0.05, type = "predicted")
head(markers_w)
class(markers_w$cluster)

# Data DE results/metrics from simulation
# list for save metrics
results_list <- list()

# Expected DEGS per group and level
results_list[["t_test"]] <- DE_metrics(expected_degs_levels_df, markers_t, all_genes_df, test="t")
results_list[["wilcoxon_test"]] <- DE_metrics(expected_degs_levels_df, markers_w, all_genes_df, test="w")
results_list[["t_test"]]$Metrics
results_list[["wilcoxon_test"]]$Metrics

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

saveWorkbook(wb, file = paste0(basedir,"/DEGs_and_Metrics.xlsx"), overwrite = TRUE)
# saveWorkbook(wb, file = paste0(basedir,"/Wilcoxon_t_Test/Results/DEGs_and_Metrics.xlsx"), overwrite = TRUE)
