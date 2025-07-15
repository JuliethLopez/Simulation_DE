## Set environment - Change for your own directory
basedir <- "C:/Users/Administrator/Documents/Jaguar/Simulation_DE/DESEQ2"
setwd("C:/Users/Administrator/Documents/Jaguar/Simulation_DE")
outdir <- "C:/Users/Administrator/Documents/Jaguar/Simulation_DE/DESEQ2/figures"
source("C:/Users/Administrator/Documents/Jaguar/Simulation_DE/Functions.R")

## Libraries
library(Seurat)
library(dplyr)
library(ggplot2)
library(reshape2)
#library(data.table)
library(openxlsx)
library(DESeq2)

#### CARGAR DATOS ####

# Cargar el archivo .rds
sim.groups <- readRDS("Simulaciones/simulation_2.rds")

# Verificar que los datos se cargaron correctamente
print(sim.groups)

## estructura
str(sim.groups)
print(sim.groups)

# Ver las anotaciones de los grupos
table(colData(sim.groups)$Group)

# Convertir datos en formato de DESeq2, si no lo están
coldata <- data.frame(Cluster = factor(colData(sim.groups)$Group, levels = c("Group1", "Group2", "Group3", "Group4", "Group5", "Group6")))
rownames(coldata) <- colnames(sim.groups)

################################################################################
# Expected DEGS per group and level
expected_degs_levels_df <- read.csv("DEGs_expected.csv")
all_genes_df <- read.csv("all_genes.csv")
all_genes_df <- all_genes_df %>% select(c(Gene, cluster))

################################################################################
#### SEGUNDA PARTE: ANÁLISIS DE EXPRESIÓN DIFERENCIAL CON DESEQ2 ####
################################################################################
# Crear el conjunto de datos DESeq2 desde la matriz de conteo
deseq2_matrix <- DESeqDataSetFromMatrix(
  countData = assay(sim.groups, "counts"), 
  colData = coldata, 
  design = ~ 0 + Cluster  # Modelo con 6 grupos/clusters
)

DE_deseq2_original <- list()

for (level in levels(deseq2_matrix$Cluster)) {
  print(level)
  # Relevel para comparar cada grupo con los otros
  deseq2_matrix$Cluster <- relevel(deseq2_matrix$Cluster, ref = level)
  dds <- DESeq(deseq2_matrix)
  
  # Realizar el análisis de expresión diferencial comparando el cluster contra todos los demás
  res <- results(dds, contrast = c(1, -1/5, -1/5, -1/5, -1/5, -1/5))
  
  # Guardar resultados en formato de data frame
  results <- as.data.frame(res)
  results$cluster <- level
  results$Gene <- rownames(results)
  
  # Añadir resultados al listado de resultados por grupo
  DE_deseq2_original[[level]] <- results
}

# Combinar los resultados en un solo dataframe
DE_deseq2_groups <- bind_rows(DE_deseq2_original)
which_col = which(colnames(DE_deseq2_groups)=="padj")
colnames(DE_deseq2_groups)[which_col] = "p_val_adj"

# Guardar los resultados en un archivo CSV
write.csv(DE_deseq2_groups, file = paste0(outdir, "/DE_deseq2_original.csv"), row.names = TRUE)

#### TERCERA PARTE: FUNCIONES DE ANÁLISIS Y CÁLCULO DE MÉTRICAS AJUSTADAS A DESEQ2 ####
DE_deseq2_groups = read.csv(paste0(basedir,"/DE_deseq2_original.csv"))

# Filtrar y procesar los resultados de DESeq2 (ajustada a DESeq2)
deseq_results_filtered = DE_classes(DE_deseq2_groups, min_pct_difference = 0.1, max_p_val_adj = 0.05, type="predicted_deseq")

# Calcular las métricas basadas en los resultados filtrados de DESeq2
# list for save metrics
results_list <- list()
# Expected DEGS per group and level
results_list[["deseq2"]] <- DE_metrics(expected_degs_levels_df, deseq_results_filtered, all_genes, test="deseq2")

### HEATMAP ####
Metrics_heatmap(results_list, outdir, text="_simulation")

# Guardar los DEGs (según los filtros aplicados), los DEGs por cluster y las métricas en un archivo CSV
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

saveWorkbook(wb, file = paste0(basedir,"/DEGs_and_Metrics_Deseq2.xlsx"), overwrite = TRUE)


