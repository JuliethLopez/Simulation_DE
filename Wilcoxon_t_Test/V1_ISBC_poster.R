##### 
rm(list=ls(all=TRUE))

setwd("/Users/mariaalejandrarojo/Desktop/02_Estadistica_genómica/ISBC_poster")
outdir <- "/Users/mariaalejandrarojo/Desktop/02_Estadistica_genómica/ISBC_poster/plots"
# https://www.iscb.org/latam2024/call-for-submissions/oral-poster-presentations
par(mfrow = c(1, 1)) 
plot.new()

## Libraries ####
# Cargar librerías necesarias
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("splatter", "Seurat", "scater", "SeuratData"))
install.packages(c("ggplot2", "tidyverse", "openxlsx"))

library(ggplot2)
library(scater)
library(splatter)
library(Seurat)
library(SeuratData)
library(tidyverse)
library(openxlsx)


## charge data Ref ####
# Article - (2017) Splatter: simulation of single-cell RNA sequencing data. https://doi.org/10.1186/s13059-017-1305-0
# Tutorial - Splatter: https://bioconductor.statistik.tu-dortmund.de/packages/3.13/bioc/vignettes/splatter/inst/doc/splatter.html

####. 1. DATOS Referencia pbmc3k #####

AvailableData()
InstallData("pbmc3k") 

pbmc3k <- LoadData("pbmc3k.SeuratData")
pbmc3k <- UpdateSeuratObject(pbmc3k)

# 1 (data)
data_0 <- pbmc3k
class(data_0) # SeuratObject

# tabla de frecuencias de los tipos de celulas 
tab0 <-round(table(data_0$seurat_annotations)/sum(table(data_0$seurat_annotations)),2)
tab0


# Extraer conteos de expresión génica
data_counts <- as.matrix(data_0@assays$RNA$counts)
class(data_counts) #matrix

dim(data_counts) #13714  2700

# Cells frecuency
#sampled_cells_indices <- sample(1:ncol(data_counts), size = 1000, replace = FALSE)
#sampled_cells_indices

# extraer parametros de conteo
params <- splatEstimate(data_counts) #  paramatros para la simulación 
params
#getParams(params, names = c("nGenes","nCells", "mean.rate", "mean.shape", "nGroups"))

# 2 (ajustar)

# centra y escalar datos 
# ajustar media de expresion ajusta a 0 y varianza a 1 
data_0 <- ScaleData(data_0)

# identificar genes variables
data_0 <- FindVariableFeatures(data_0)

#3 (cluster y PCA)

# PCA 
data_0_PCA <- RunPCA(data_0, features = VariableFeatures(object = data_0))
# Plot PCA de seurat
DimPlot(data_0_PCA, reduction = "pca") + ggtitle("Datos Reales (PCA)")

#Cluster the cells
#kNN
data_0_kNN <- FindNeighbors(data_0_PCA, dims = 1:8)
# utiliza el grafo KNN para cluster 
data_0_cl <- FindClusters(data_0_kNN, resolution = 0.5)

# Visualization of cluster results 

DimPlot(data_0_cl, group.by = "seurat_clusters")   # visualizar PCA despues de cluster

# UMAP
data_0_UMAP <- RunUMAP(data_0_cl, dims = 1:10)
#
DimPlot(data_0_UMAP, reduction = "umap", pt.size = 0.5) + NoLegend() # UMAP

###### DEGs Ref ######
# opcional

# Encontrar todos los genes diferencialmente expresados (DEGs) entre clústeres
all_DEGs <- FindAllMarkers(data_0_cl, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 1.5)

# Categorizar los DEGs según niveles de logfc
all_DEGs <- all_DEGs %>%
  mutate(logfc_level = case_when(
    avg_log2FC >= 1.2 & avg_log2FC < 2 ~ "bajo",
    avg_log2FC >= 2 & avg_log2FC < 3 ~ "medio",
    avg_log2FC >= 3 ~ "alto",
    TRUE ~ "otros"
  ))

# Filtrar DEGs según los niveles de logfc
DEGs_bajo <- all_DEGs %>% filter(logfc_level == "bajo")
DEGs_medio <- all_DEGs %>% filter(logfc_level == "medio")
DEGs_alto <- all_DEGs %>% filter(logfc_level == "alto")

# Mostrar el número de DEGs por nivel de logfc
DEGs_summary <- all_DEGs %>%
  group_by(logfc_level) %>%
  summarize(num_DEGs = n())

print(DEGs_summary)

# Mostrar los primeros resultados
head(all_DEGs)
tail(all_DEGs)

# Contar el número de DEGs por clúster
DEGs_by_cluster <- all_DEGs %>%
  group_by(cluster) %>%
  summarize(num_DEGs = n())

# Mostrar el número de DEGs por clúster
print(DEGs_by_cluster)

# Identificar DEGs entre clústeres 2 cluster 
#cluster_DEGs <- FindMarkers(data_0_cl, ident.1 = 0, ident.2 = 1)  # Comparar clúster 0 vs clúster 1
# Mostrar los primeros resultados
#head(cluster_DEGs)

# entre 1 y todos
# Encuentra genes diferencialmente expresados entre un clúster específico y todos los demás
#cluster_1_vs_all_DEGs <- FindMarkers(data_0_cl, ident.1 = 1)

# Visualiza algunos de los resultados
#head(cluster_1_vs_all_DEGs)
#dim(cluster_1_vs_all_DEGs)



##### guardar resultados ####

if (!require("openxlsx", quietly = TRUE))
  install.packages("openxlsx")

library(openxlsx)

# Crear un nuevo workbook
wb <- createWorkbook()

# Añadir una hoja para cada categoría y escribir los datos
addWorksheet(wb, "DEGs Bajo")
writeData(wb, sheet = "DEGs Bajo", DEGs_bajo)

addWorksheet(wb, "DEGs Medio")
writeData(wb, sheet = "DEGs Medio", DEGs_medio)

addWorksheet(wb, "DEGs Alto")
writeData(wb, sheet = "DEGs Alto", DEGs_alto)

# Guardar el workbook en un archivo
saveWorkbook(wb, file = "DEGs_by_logfc_level.xlsx", overwrite = TRUE)

#### simulation scRNA data ONE individual with different types of cells ####

set.seed(42)
#1
# Simulación de datos
#sim.groups <- splatSimulate(
#params, 
#group.prob = c(0.45, 0.23, 0.15, 0.10, 0.04, 0.03), # Probabilidades de pertenencia a cada grupo
#de.prob = c(0.10, 0.09, 0.07, 0.06, 0.04, 0.04),    # Probabilidades de expresión diferencial para cada grupo
#de.facLoc = c(0.4, 1, 0.8, 0.8, 0.5, 1),            # Localización del factor de expresión diferencial, ajustado para cada grupo
#de.facScale = c(0.9, 0.9, 0.7, 0.8, 0.8, 0.9),      # Escala del factor de expresión diferencial, ajustado para cada grupo
#batchCells = 3000,                                  # Número de células por lote en la simulación (3000)
#method = "groups",                                  # Método de simulación (groups)
#verbose = TRUE
#)


# Datos de referencia
#data_counts <- as.matrix(data_0@assays$RNA$counts)

# Ajuste de los parámetros de simulación mejorada
#params_mejorada <- splatEstimate(data_counts)
#params_mejorada <- setParams(params_mejorada, dropout.type = "none", dropout.mid = -0.15, dropout.shape = -1.06)

# Ejecutar una nueva simulación con los parámetros mejorados
#set.seed(42)

#sim.groups <- splatSimulate(
# params_mejorada, 
#group.prob = c(0.45, 0.23, 0.15, 0.10, 0.04, 0.03), # Probabilidades de pertenencia a cada grupo
#de.prob = c(0.10, 0.09, 0.07, 0.06, 0.04, 0.04),    # Probabilidades de expresión diferencial para cada grupo
#de.facLoc = c(0.4, 1, 0.8, 0.8, 0.5, 1),            # Localización del factor de expresión diferencial, ajustado para cada grupo
#de.facScale = c(0.9, 0.9, 0.7, 0.8, 0.8, 0.9),      # Escala del factor de expresión diferencial ajustada
#batchCells = 4000,                                  # Número de células por lote en la simulación (6000)
#nGenes = 15000,                                     # Incremento del número de genes
#method = "groups",                                  # Método de simulación (groups)
#verbose = TRUE
#)

# save simulation
saveRDS(sim.groups, file = "simulation_3.rds")

#### explorar simnulación  

#sim.groups <- readRDS("simulation_Final.rds")
sim.groups <- readRDS("simulation_2.rds")

# Verificar que los datos se cargaron correctamente
print(sim.groups)

## estructura
str(sim.groups)
print(sim.groups)

# Ver las anotaciones de los grupos
table(colData(sim.groups)$Group)

#### CARGAR DATOS ####

# Cargar el archivo .rds

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
# Análisis de expresión diferencial usando Wilcoxon
comparison_result_wilcoxonSIM <- FindAllMarkers(data_1_cl, test.use = "wilcox", only.pos = TRUE, logfc.threshold = log2(1.5))

# Crear la tabla de resumen para Wilcoxon y categorizar DEGs según niveles de logfc
summary_over_df_wilcoxonSIM <- comparison_result_wilcoxonSIM %>%
  filter(p_val_adj < 0.05) %>%
  mutate(Gene = rownames(comparison_result_wilcoxonSIM)[1:n()]) %>%
  select(Gene, cluster, avg_log2FC) %>%
  rename(DEFacGroup = avg_log2FC) %>%
  mutate(logfc_level = case_when(
    DEFacGroup >= 1.2 & DEFacGroup < 2 ~ "bajo",
    DEFacGroup >= 2 & DEFacGroup < 3 ~ "medio",
    DEFacGroup >= 3 ~ "alto",
    TRUE ~ "otros"
  ))

# Filtrar DEGs según los niveles de logfc
DEGs_bajo_wilcoxonSIM <- summary_over_df_wilcoxonSIM %>% filter(logfc_level == "bajo")
DEGs_medio_wilcoxonSIM <- summary_over_df_wilcoxonSIM %>% filter(logfc_level == "medio")
DEGs_alto_wilcoxonSIM <- summary_over_df_wilcoxonSIM %>% filter(logfc_level == "alto")

# Contar el número de genes en cada categoría
num_bajo <- nrow(DEGs_bajo_wilcoxonSIM)
num_medio <- nrow(DEGs_medio_wilcoxonSIM)
num_alto <- nrow(DEGs_alto_wilcoxonSIM)

# Imprimir los resultados
cat("Número de genes con logfc bajo:", num_bajo, "\n")
cat("Número de genes con logfc medio:", num_medio, "\n")
cat("Número de genes con logfc alto:", num_alto, "\n")

# Cálculo de métricas de evaluación
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

# Suponiendo que tenemos un vector `actual` que contiene los valores verdaderos para la clasificación de los genes
# y `predicted` que contiene los valores predichos por el análisis diferencial
# Aquí sólo se proporcionan ejemplos. Debes reemplazarlo con tus propios datos.
# actual <- ...
# predicted <- ...

# Ejemplo de cómo calcular las métricas (esto debe ser reemplazado con tus datos reales)
# metrics_bajo <- calculate_metrics(predicted_bajo, actual_bajo)
# metrics_medio <- calculate_metrics(predicted_medio, actual_medio)
# metrics_alto <- calculate_metrics(predicted_alto, actual_alto)

# Instalar y cargar openxlsx si no está instalado
if (!require("openxlsx", quietly = TRUE)) install.packages("openxlsx")
library(openxlsx)

# Crear un nuevo workbook
wb_sim <- createWorkbook()

# Añadir una hoja para cada categoría y escribir los datos
addWorksheet(wb_sim, "DEGs Bajo Wilcoxon")
writeData(wb_sim, sheet = "DEGs Bajo Wilcoxon", DEGs_bajo_wilcoxonSIM)

addWorksheet(wb_sim, "DEGs Medio Wilcoxon")
writeData(wb_sim, sheet = "DEGs Medio Wilcoxon", DEGs_medio_wilcoxonSIM)

addWorksheet(wb_sim, "DEGs Alto Wilcoxon")
writeData(wb_sim, sheet = "DEGs Alto Wilcoxon", DEGs_alto_wilcoxonSIM)

# Guardar el workbook en un archivo
saveWorkbook(wb_sim, file = "DEGs_by_logfc_level_Wilcoxon_SIM.xlsx", overwrite = TRUE)

# Mostrar los primeros resultados de la comparación con Wilcoxon
head(comparison_result_wilcoxonSIM)

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

