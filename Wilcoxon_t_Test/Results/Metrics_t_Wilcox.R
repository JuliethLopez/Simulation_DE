## Set environment - Change for your own directory
basedir <- "C:/Users/Administrator/Documents/Jaguar/Simulation_DE/Wilcoxon_t_Test"
setwd("C:/Users/Administrator/Documents/Jaguar/Simulation_DE")
outdir <- "C:/Users/Administrator/Documents/Jaguar/Simulation_DE/Wilcoxon_t_Test/Results"
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
DE_t_test <- DE_analysis_seurat(data_1_cl, assay = "RNA", slot="data", test_type = "t")
head(DE_t_test)
dim(DE_t_test) #1679    7

DE_t_test[which(DE_t_test$gene == "Gene3"),]

## W
DE_Wilcox <- DE_analysis_seurat(data_1_cl, assay = "RNA", slot="data", test_type = "wilcox")
head(DE_Wilcox)
dim(DE_Wilcox) #2714    7

DE_Wilcox[which(DE_Wilcox$gene == "Gene3"),]

# save results
write.csv(DE_Wilcox, file = paste0(outdir,"/DE_results_wilcox.csv"), row.names = TRUE) #save markers
write.csv(DE_t_test, file = paste0(outdir,"/DE_results_t.csv"), row.names = TRUE) #save markers

# Cargar los archivos de resultados
DE_t_test <- read.csv(paste0(outdir,"/DE_results_t.csv"))
DE_Wilcox <- read.csv(paste0(outdir,"/DE_results_wilcox.csv"))

# Aplicar la función DE_classes 
markers_t <- DE_classes(DE_t_test, min_pct_difference = 0.1, max_p_val_adj = 0.05, type = "predicted")
head(markers_t)
class(markers_t$cluster)

markers_wilcox <- DE_classes(DE_Wilcox, min_pct_difference = 0.1, max_p_val_adj = 0.05, type = "predicted")
head(markers_wilcox)
class(markers_wilcox$cluster)

# Calcular las métricas para los resultados de t-test y Wilcoxon
results_t <- DE_metrics(expected_degs_levels_df, markers_t, all_genes_df, test="t-test")
results_wilcox <- DE_metrics(expected_degs_levels_df, markers_wilcox, all_genes_df, test="wilcox")

results_list <- list()
results_list[["t_test"]] <- results_t
results_list[["wilcoxon_test"]] <- results_wilcox

# Crear y guardar los heatmaps de las métricas para t-test y Wilcoxon
#outdir <- "resultados" # Asegúrate de tener este directorio creado
Metrics_heatmap(results_list, outdir, text = "")

# Guardar los resultados de las métricas en archivos CSV
# write.csv(results_t$Metrics, paste0(basedir,"/Wilcoxon_t_Test/Results/Metrics_t_test.csv"), row.names = FALSE)
# write.csv(results_wilcox$Metrics, paste0(basedir,"/Wilcoxon_t_Test/Results/Metrics_wilcox.csv"), row.names = FALSE)

# Guardar los genes clasificados para t-test y Wilcoxon
# write.csv(results_t$DEGs,  paste0(outdir,"Wilcoxon_t_Test/DEGs_t_test.csv"), row.names = FALSE)
# write.csv(results_wilcox$DEGs,  paste0(outdir,"/DEGs_wilcox.csv"), row.names = FALSE)

## Save DEGs, number the DEGs per cluster and the metrics in an Excel file
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

saveWorkbook(wb, file = paste0(outdir,"/DEGs_and_Metrics.xlsx"), overwrite = TRUE)

# Revise metrics for both t and wilcoxon.Both are pretty similar but not equal!
head(results_wilcox$Metrics)
head(results_t$Metrics)
