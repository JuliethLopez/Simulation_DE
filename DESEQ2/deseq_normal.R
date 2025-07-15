###### Instalaciones ###########

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

devtools::install_github('satijalab/seurat-data')

BiocManager::install("splatter") 
BiocManager::install("Seurat") 
BiocManager::install("DESeq2") 
BiocManager::install("scater")
BiocManager::install("SeuratData")
BiocManager::install("Seurat")
BiocManager::install("SingleCellExperiment")
BiocManager::install("scran")
BiocManager::install("splat")
BiocManager::install("Rtools")
BiocManager::install("caret")
BiocManager::install("grid")
BiocManager::install("gridExtra")
BiocManager::install("ggpubr")

BiocManager::install(
  "Oshlack/splatter",
  dependencies = TRUE,
  build_vignettes = TRUE
)

install.packages("patchwork")
install.packages("ggplot2")
install.packages("tidyverse")
install.packages("devtools")
install.packages("https://seurat.nygenome.org/src/contrib/pbmc3k.SeuratData_3.0.0.tar.gz", repos = NULL, type = "source")
install.packages("C:/Users/Usuario/OneDrive/Escritorio/proyectoestadistica/mio/conflicted-main", repos = NULL, type = "source")
install.packages("C:/Users/Usuario/OneDrive/Escritorio/proyectoestadistica/mio/seurat-data-master", repos = NULL, type = "source")

###### llamado de librerias#####
library(caret)

suppressPackageStartupMessages({
  library(splatter)
})

library(dplyr)
library(Seurat)
library(patchwork)
library(BiocFileCache)
library(R.utils)
library(Matrix)
library(SingleCellExperiment)
library(scran)
library("scater") #normalization and pca
library("ggplot2") #plotting
library(patchwork)
library(DESeq2)
library(tibble)
library(grid)
library(gridExtra)
library(ggpubr)
library("tidyverse")
library("SeuratData") #data reference PBMC

outdir <- "C:/Users/Usuario/OneDrive/Escritorio/proyectoestadistica"

#### CARGAR DATOS ####

# Cargar el archivo .rds
sim.groups <- readRDS("C:/Users/Usuario/OneDrive/Escritorio/proyectoestadistica/simulation_2.rds")

# Verificar que los datos se cargaron correctamente
print(sim.groups)

## estructura
str(sim.groups)
print(sim.groups)

# Ver las anotaciones de los grupos
table(colData(sim.groups)$Group)

#### PREPROCESAMIENTO DE DATOS ####

# Preprocesamiento de sim.groups
matrix_counts <- sim.groups@assays@data@listData[["counts"]] # Not normalized
data_1 <- CreateSeuratObject(counts = matrix_counts)
data_1 <- NormalizeData(data_1, normalization.method = "LogNormalize", scale.factor = 10000)
data_1 <- FindVariableFeatures(data_1, selection.method = "vst", nfeatures = 2000)
data_1 <- ScaleData(data_1)
data_1_pca <- RunPCA(data_1, features = VariableFeatures(object = data_1))
ElbowPlot(data_1_pca)
data_1_SNN <- FindNeighbors(data_1_pca, dims = 1:15) # Ajustar el número de dimensiones
data_1_cl <- FindClusters(data_1_SNN, resolution = 0.7) # Ajustar la resolución
data_1_UMAP <- RunUMAP(data_1_cl, dims = 1:12)
DimPlot(data_1_UMAP, label = TRUE, reduction = "umap", pt.size = 0.5) + NoLegend()

# Verificación de la distribución de los grupos
table(Idents(data_1_cl)) # Número de células por clúster
table(colData(sim.groups)$Group) # Verificar la distribución de los grupos simulados

###########################
# Análisis de expresión diferencial usando DESeq2

# Crear el objeto DESeqDataSet a partir de los datos Seurat
library(DESeq2)

dds <- DESeqDataSetFromMatrix(
  countData = matrix_counts, 
  colData = data.frame(cluster = Idents(data_1_cl)), 
  design = ~ cluster
)

# Correr DESeq2
dds <- DESeq(dds)

# Obtener resultados y filtrar por un valor ajustado de p < 0.05
res <- results(dds, contrast = c("cluster", "1", "2"), alpha = 0.05)
res <- res[!is.na(res$padj) & res$padj < 0.05, ]

# Crear la tabla de resumen para DESeq2 y categorizar DEGs según niveles de logfc
summary_over_df_deseq2SIM <- data.frame(res) %>%
  rownames_to_column(var = "Gene") %>%
  mutate(DEFacGroup = log2FoldChange) %>%
  mutate(logfc_level = case_when(
    DEFacGroup >= 1.2 & DEFacGroup < 2 ~ "bajo",
    DEFacGroup >= 2 & DEFacGroup < 3 ~ "medio",
    DEFacGroup >= 3 ~ "alto",
    TRUE ~ "otros"
  ))

# Filtrar DEGs según los niveles de logfc
DEGs_bajo_deseq2SIM <- summary_over_df_deseq2SIM %>% filter(logfc_level == "bajo")
DEGs_medio_deseq2SIM <- summary_over_df_deseq2SIM %>% filter(logfc_level == "medio")
DEGs_alto_deseq2SIM <- summary_over_df_deseq2SIM %>% filter(logfc_level == "alto")

# Contar el número de genes en cada categoría
num_bajo <- nrow(DEGs_bajo_deseq2SIM)
num_medio <- nrow(DEGs_medio_deseq2SIM)
num_alto <- nrow(DEGs_alto_deseq2SIM)

# Imprimir los resultados
cat("Número de genes con logfc bajo:", num_bajo, "\n")
cat("Número de genes con logfc medio:", num_medio, "\n")
cat("Número de genes con logfc alto:", num_alto, "\n")

# Cálculo de métricas de evaluación (manteniendo la función original)
calculate_metrics <- function(predicted, actual) {
  TP <- sum(predicted == 1 & actual == 1)
  TN <- sum(predicted == 0 & actual == 0)
  FP <- sum(predicted == 1 & actual == 0)
  FN <- sum(predicted == 0 & actual == 1)
  
  precision <- TP / (TP + FP)
  recall <- TP / (TP + FN)
  f1_score <- 2 * (precision * recall) / (precision + recall)
  accuracy <- (TP + TN) / (TP + TN + FP + FN)
  fdr <- FP / (FP + TP)
  
  return(list(precision = precision, recall = recall, f1_score = f1_score, accuracy = accuracy, fdr = fdr))
}

# Supongamos que tenemos vectores `actual` y `predicted` para cada nivel de logfc
# Aquí deberías definir los vectores `actual` y `predicted` basados en tus datos reales.
set.seed(42)
actual_bajo <- sample(c(0, 1), size = num_bajo, replace = TRUE)
predicted_bajo <- sample(c(0, 1), size = num_bajo, replace = TRUE)

actual_medio <- sample(c(0, 1), size = num_medio, replace = TRUE)
predicted_medio <- sample(c(0, 1), size = num_medio, replace = TRUE)

actual_alto <- sample(c(0, 1), size = num_alto, replace = TRUE)
predicted_alto <- sample(c(0, 1), size = num_alto, replace = TRUE)

# Calcular métricas para cada nivel
metrics_bajo <- calculate_metrics(predicted_bajo, actual_bajo)
metrics_medio <- calculate_metrics(predicted_medio, actual_medio)
metrics_alto <- calculate_metrics(predicted_alto, actual_alto)

# Crear una tabla con los resultados
metrics_table <- data.frame(
  Level = c("Bajo", "Medio", "Alto"),
  Precision = c(metrics_bajo[1], metrics_medio[1], metrics_alto[1]),
  Recall = c(metrics_bajo[2], metrics_medio[2], metrics_alto[2]),
  F1_Score = c(metrics_bajo[3], metrics_medio[3], metrics_alto[3]),
  Accuracy = c(metrics_bajo[4], metrics_medio[4], metrics_alto[4]),
  FDR = c(metrics_bajo[5], metrics_medio[5], metrics_alto[5])
)

# Mostrar la tabla de métricas
print(metrics_table) 

# Guardar los resultados en Excel
if (!require("openxlsx", quietly = TRUE)) install.packages("openxlsx")
library(openxlsx)

wb_sim <- createWorkbook()

# Añadir una hoja para cada categoría y escribir los datos
addWorksheet(wb_sim, "DEGs Bajo DESeq2")
writeData(wb_sim, sheet = "DEGs Bajo DESeq2", DEGs_bajo_deseq2SIM)

addWorksheet(wb_sim, "DEGs Medio DESeq2")
writeData(wb_sim, sheet = "DEGs Medio DESeq2", DEGs_medio_deseq2SIM)

addWorksheet(wb_sim, "DEGs Alto DESeq2")
writeData(wb_sim, sheet = "DEGs Alto DESeq2", DEGs_alto_deseq2SIM)

saveWorkbook(wb_sim, file = "DEGs_by_logfc_level_DESeq2_SIM.xlsx", overwrite = TRUE)

# Mostrar los primeros resultados de la comparación con DESeq2
head(comparison_result_DESeq2SIM)

###### metricas #####

# Función para calcular las métricas
calculate_metrics <- function(predicted, actual) {
  TP <- sum(predicted == 1 & actual == 1)
  TN <- sum(predicted == 0 & actual == 0)
  FP <- sum(predicted == 1 & actual == 0)
  FN <- sum(predicted == 0 & actual == 1)
  
  precision <- ifelse(TP + FP > 0, TP / (TP + FP), 0)
  recall <- ifelse(TP + FN > 0, TP / (TP + FN), 0)
  f1_score <- ifelse(precision + recall > 0, 2 * (precision * recall) / (precision + recall), 0)
  accuracy <- (TP + TN) / (TP + TN + FP + FN)
  fdr <- ifelse(FP + TP > 0, FP / (FP + TP), 0)
  
  return(c(precision, recall, f1_score, accuracy, fdr))
}

# Supongamos que tenemos vectores `actual` y `predicted` para cada nivel de logfc
# Aquí deberías definir los vectores `actual` y `predicted` basados en tus datos reales.
# Ejemplo de vectores de prueba (reemplaza esto con tus datos reales)
set.seed(42)
actual_bajo <- sample(c(0, 1), size = num_bajo, replace = TRUE)
predicted_bajo <- sample(c(0, 1), size = num_bajo, replace = TRUE)

actual_medio <- sample(c(0, 1), size = num_medio, replace = TRUE)
predicted_medio <- sample(c(0, 1), size = num_medio, replace = TRUE)

actual_alto <- sample(c(0, 1), size = num_alto, replace = TRUE)
predicted_alto <- sample(c(0, 1), size = num_alto, replace = TRUE)

# Calcular métricas para cada nivel
metrics_bajo <- calculate_metrics(predicted_bajo, actual_bajo)
metrics_medio <- calculate_metrics(predicted_medio, actual_medio)
metrics_alto <- calculate_metrics(predicted_alto, actual_alto)

# Crear una tabla con los resultados
metrics_table <- data.frame(
  Level = c("Bajo", "Medio", "Alto"),
  Precision = c(metrics_bajo[1], metrics_medio[1], metrics_alto[1]),
  Recall = c(metrics_bajo[2], metrics_medio[2], metrics_alto[2]),
  F1_Score = c(metrics_bajo[3], metrics_medio[3], metrics_alto[3]),
  Accuracy = c(metrics_bajo[4], metrics_medio[4], metrics_alto[4]),
  FDR = c(metrics_bajo[5], metrics_medio[5], metrics_alto[5])
)

# Mostrar la tabla de métricas
print(metrics_table) 

