##### 
rm(list=ls(all=TRUE))

setwd("/Users/mariaalejandrarojo/Desktop/02_Estadistica_genómica/06_proyecto_SG/Proyecto_SingleCell_Est.Gen")
outdir <- "/Users/mariaalejandrarojo/Desktop/02_Estadistica_genómica/06_proyecto_SG/Proyecto_SingleCell_Est.Gen/plots"
par(mfrow = c(1, 1)) 
plot.new()

## Libraries ####

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("splatter") 
BiocManager::install("Seurat") 

#BiocManager::install("DESeq2") 

install.packages("ggplot2")
install.packages("tidyverse")
BiocManager::install("scater")
BiocManager::install("SeuratData")
devtools::install_github('satijalab/seurat-data')

#install.packages("patchwork")
#source("https://bioconductor.org/biocLite.R")
#BiocManager::install("scater", forse = TRUE)

#devtools::install_github('satijalab/seurat-data') ok 
#install.packages("Seurat")

# l
library("ggplot2") #plotting
library("scater") #normalization and pca
library("splatter") #simulation 
library("Seurat") #single cell pre-processing and processing PBMC Dat0
library("SeuratData") #data reference PBMC
library("tidyverse") #database manipulation

library(SeuratData)
library(Matrix)
library(scater)
library(SingleCellExperiment)


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

#### MAST REF #####

# ejecutar primero -->  realizar_comparacion_todos

comparison_result_mastREF <- realizar_comparacion_todos(data_0_cl, test_use = "MAST", only_pos = TRUE, logfc_threshold = log2(1.5))
# Resumen de la cantidad de DEGs por clúster para el método MAST
deg_summary_mastREF <- comparison_result_mastREF %>%
  filter(p_val_adj < 0.05) %>%
  group_by(cluster) %>%
  summarise(Num_DEGs = n())

# Asegurar que todos los clústeres estén representados, incluso si no tienen DEGs significativos
all_clusters <- sort(unique(Idents(data_0_cl)))
deg_summary_mastREF <- merge(data.frame(Cluster = all_clusters), deg_summary_mastREF, by.x = "Cluster", by.y = "cluster", all.x = TRUE)
deg_summary_mastREF$Num_DEGs[is.na(deg_summary_mastREF$Num_DEGs)] <- 0

# Crear la tabla de resumen para MAST
summary_over_df_mast <- comparison_result_mastREF %>%
  filter(p_val_adj < 0.05) %>%
  mutate(Gene = rownames(comparison_result_mastREF)) %>%
  select(Gene, cluster, avg_log2FC) %>%
  rename(DEFacGroup = avg_log2FC)

# Ordenar y mostrar el resumen de DEGs por comparación
deg_summary_mastREF <- deg_summary_mastREF[order(deg_summary_mastREF$Cluster), ]
print(deg_summary_mastREF)
head(comparison_result_mastREF)


####################### SIMULACIÓN ##################################################

#### simulation scRNA data ONE individual with different types of cells ####

set.seed(42)

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
data_counts <- as.matrix(data_0@assays$RNA$counts)

# Ajuste de los parámetros de simulación mejorada
params_mejorada <- splatEstimate(data_counts)
params_mejorada <- setParams(params_mejorada, dropout.type = "none", dropout.mid = -0.15, dropout.shape = -1.06)

# Ejecutar una nueva simulación con los parámetros mejorados
set.seed(42)

sim.groups <- splatSimulate(
  params_mejorada, 
  group.prob = c(0.45, 0.23, 0.15, 0.10, 0.04, 0.03), # Probabilidades de pertenencia a cada grupo
  de.prob = c(0.10, 0.09, 0.07, 0.06, 0.04, 0.04),    # Probabilidades de expresión diferencial para cada grupo
  de.facLoc = c(0.4, 1, 0.8, 0.8, 0.5, 1),            # Localización del factor de expresión diferencial, ajustado para cada grupo
  de.facScale = c(0.9, 0.9, 0.7, 0.8, 0.8, 0.9),      # Escala del factor de expresión diferencial ajustada
  batchCells = 4000,                                  # Número de células por lote en la simulación (6000)
  nGenes = 15000,                                     # Incremento del número de genes
  method = "groups",                                  # Método de simulación (groups)
  verbose = TRUE
)

# save simulation
saveRDS(sim.groups, file = "simulation_Final.rds")


#### explorar simnulación  

## estructura
str(sim.groups)
print(sim.groups)

# Ver las anotaciones de los grupos
table(colData(sim.groups)$Group)

##### Pre- procesamiento ######

matrix_counts <- sim.groups @assays@data@listData[["counts"]] # Not normalized
data_1 <- CreateSeuratObject(counts = matrix_counts)
data_1 <- NormalizeData(data_1, normalization.method = "LogNormalize", scale.factor = 10000)
data_1 <- FindVariableFeatures(data_1, selection.method = "vst", nfeatures = 2000)
data_1 <- ScaleData(data_1)
data_1_pca <- RunPCA(data_1, features = VariableFeatures(object = data_1))
ElbowPlot(data_1_pca)
data_1_SNN <- FindNeighbors(data_1_pca, dims = 1:19) # Ajustar el número de dimensiones
data_1_cl <- FindClusters(data_1_SNN, resolution = 0.8) # Ajustar la resolución
data_1_UMAP <- RunUMAP(data_1_cl, dims = 1:10)
data_1_cl <- AddMetaData(data_1_cl, metadata = Idents(object = data_1_cl), col.name = "Cluster")
DimPlot(data_1_UMAP, label = TRUE, reduction = "umap", pt.size = 0.5) + NoLegend()


table(Idents(data_1_cl)) #numro de celulas por cluster 
table(colData(sim.groups)$Group)# Verificar la distribución de los grupos simulados
#
### 1 ################### DEFacGroup ######

# Obtener los datos de los factores de expresión diferencial para cada grupo
de_fac_groups <- lapply(1:6, function(i) {
  sim.groups@rowRanges@elementMetadata[[paste0("DEFacGroup", i)]]
})

# Calcular los valores mínimos y máximos para cada grupo
min_max_values <- data.frame(
  Group = paste0("DEFacGroup", 1:6),
  Min = sapply(de_fac_groups, min),
  Max = sapply(de_fac_groups, max)
)

# Mostrar los resultados
print(min_max_values)

# Definir genes esperados para cada clúster
expected_degs_filtered <- list()
for (i in 1:6) {
  cluster_name <- paste0("Cluster_", i-1)
  expected_degs_filtered[[cluster_name]] <- rownames(sim.groups)[rowData(sim.groups)[[paste0("DEFacGroup", i)]] > 1.5]
}

# Crear una tabla con el número de DEGs esperados en cada clúster
deg_counts_filtered <- sapply(expected_degs_filtered, length)
deg_counts_df_filtered_1 <- data.frame(
  Cluster = names(deg_counts_filtered),
  Num_DEGs_Overexpressed = deg_counts_filtered
)

# Mostrar la tabla con el número de DEGs esperados
print(deg_counts_df_filtered_1) # DeFacGruop Foldchange > 1.5

# Función para crear tablas de resumen
crear_resumen <- function(expected_degs) {
  gene_names <- c()
  clusters <- c()
  defac_groups <- c()
  
  for (i in 1:6) {
    cluster_name <- paste0("Cluster_", i-1)
    group_name <- paste0("DEFacGroup", i)
    
    genes <- expected_degs[[cluster_name]]
    
    gene_names <- c(gene_names, genes)
    clusters <- c(clusters, rep(cluster_name, length(genes)))
    defac_groups <- c(defac_groups, rowData(sim.groups)[genes, group_name])
  }
  
  data.frame(
    Gene = gene_names,
    Cluster = clusters,
    DEFacGroup = defac_groups
  )
}

# Crear tabla de resumen para DEFacGroup > 1.5
summary_over_df_filtered_1 <- crear_resumen(expected_degs_filtered)

# Mostrar las tablas de resumen
print(head(summary_over_df_filtered_1)) # DeFacGruop Foldchange > 1.5
print(dim(summary_over_df_filtered_1)) # DeFacGruop Foldchange > 1.5



##### comparacion con datos reales ####

comparison <- compareSCEs(list(Splat = sim.groups, Real = as.SingleCellExperiment(data_0))) 
plot_means_comparison <- comparison$Plots$Meansplot_means_comparison <- comparison$Plots$Means
plot_means_comparison 
#ggsave(filename = paste0(outdir,"/Mean_expression_simulation_vs_realdata_1.png"), plot = plot_means_comparison)

if (!requireNamespace("MAST", quietly = TRUE)) {
  BiocManager::install("MAST")
  
}

# Actualizar todos los paquetes instalados
update.packages(ask = FALSE)

# Reinstalar MAST y sus dependencias
BiocManager::install("MAST", type = "source")


# Función para realizar el análisis de DEGs utilizando diferentes métodos de prueba
realizar_comparacion_todos <- function(data, test_use = "wilcox", only_pos = TRUE, min_pct = 0, logfc_threshold = 0) {
  comparison_result <- FindAllMarkers(data, 
                                      only.pos = only_pos, 
                                      min.pct = min_pct, 
                                      logfc.threshold = logfc_threshold, 
                                      test.use = test_use)
  comparison_result$p_val_adj <- p.adjust(comparison_result$p_val, method = "fdr")
  return(comparison_result)
}

# Realizar comparación para DEFacGroup > 1.5 (solo sobreexpresados) usando diferentes métodos de prueba
start_time <- Sys.time()
comparison_result_wilcox <- realizar_comparacion_todos(data_1_cl, test_use = "wilcox", only_pos = TRUE, logfc_threshold = log2(1.5))
comparison_result_t_test <- realizar_comparacion_todos(data_1_cl, test_use = "t", only_pos = TRUE, logfc_threshold = log2(1.5))
comparison_result_deseq2 <- realizar_comparacion_todos(data_1_cl, test_use = "DESeq2", only_pos = TRUE, logfc_threshold = log2(1.5))
comparison_result_mast <- realizar_comparacion_todos(data_1_cl, test_use = "MAST", only_pos = TRUE, logfc_threshold = log2(1.5))
end_time <- Sys.time()

# Resumen de la cantidad de DEGs por clúster para cada caso DEFacGroup > 1.5
deg_summary_wilcox <- comparison_result_wilcox %>%
  filter(p_val_adj < 0.05) %>%
  group_by(cluster) %>%
  summarise(Num_DEGs = n())

deg_summary_t_test <- comparison_result_t_test %>%
  filter(p_val_adj < 0.05) %>%
  group_by(cluster) %>%
  summarise(Num_DEGs = n())

deg_summary_deseq2 <- comparison_result_deseq2 %>%
  filter(p_val_adj < 0.05) %>%
  group_by(cluster) %>%
  summarise(Num_DEGs = n())

deg_summary_mast <- comparison_result_mast %>%
  filter(p_val_adj < 0.05) %>%
  group_by(cluster) %>%
  summarise(Num_DEGs = n())

# Si hay clusters que no tienen DEGs significativos, los incluimos con valor 0
all_clusters <- sort(unique(Idents(data_1_cl)))
deg_summary_wilcox <- merge(data.frame(Cluster = all_clusters), deg_summary_wilcox, by.x = "Cluster", by.y = "cluster", all.x = TRUE)
deg_summary_wilcox$Num_DEGs[is.na(deg_summary_wilcox$Num_DEGs)] <- 0

deg_summary_t_test <- merge(data.frame(Cluster = all_clusters), deg_summary_t_test, by.x = "Cluster", by.y = "cluster", all.x = TRUE)
deg_summary_t_test$Num_DEGs[is.na(deg_summary_t_test$Num_DEGs)] <- 0

deg_summary_deseq2 <- merge(data.frame(Cluster = all_clusters), deg_summary_deseq2, by.x = "Cluster", by.y = "cluster", all.x = TRUE)
deg_summary_deseq2$Num_DEGs[is.na(deg_summary_deseq2$Num_DEGs)] <- 0

deg_summary_mast <- merge(data.frame(Cluster = all_clusters), deg_summary_mast, by.x = "Cluster", by.y = "cluster", all.x = TRUE)
deg_summary_mast$Num_DEGs[is.na(deg_summary_mast$Num_DEGs)] <- 0

# Crear las tablas de resumen para FindAllMarkers
summary_over_df_wilcox <- comparison_result_wilcox %>%
  filter(p_val_adj < 0.05) %>%
  mutate(Gene = rownames(comparison_result_wilcox)) %>%
  select(Gene, cluster, avg_log2FC) %>%
  rename(DEFacGroup = avg_log2FC)

summary_over_df_t_test <- comparison_result_t_test %>%
  filter(p_val_adj < 0.05) %>%
  mutate(Gene = rownames(comparison_result_t_test)) %>%
  select(Gene, cluster, avg_log2FC) %>%
  rename(DEFacGroup = avg_log2FC)

summary_over_df_deseq2 <- comparison_result_deseq2 %>%
  filter(p_val_adj < 0.05) %>%
  mutate(Gene = rownames(comparison_result_deseq2)) %>%
  select(Gene, cluster, avg_log2FC) %>%
  rename(DEFacGroup = avg_log2FC)

summary_over_df_mast <- comparison_result_mast %>%
  filter(p_val_adj < 0.05) %>%
  mutate(Gene = rownames(comparison_result_mast)) %>%
  select(Gene, cluster, avg_log2FC) %>%
  rename(DEFacGroup = avg_log2FC)

# Ordenar y mostrar los resúmenes de DEGs por comparación
deg_summary_wilcox <- deg_summary_wilcox[order(deg_summary_wilcox$Cluster), ]
deg_summary_t_test <- deg_summary_t_test[order(deg_summary_t_test$Cluster), ]
deg_summary_deseq2 <- deg_summary_deseq2[order(deg_summary_deseq2$Cluster), ]
deg_summary_mast <- deg_summary_mast[order(deg_summary_mast$Cluster), ]

print(deg_summary_wilcox) # FindAllMarkers - Wilcoxon
print(deg_summary_t_test) # FindAllMarkers - T-Test
print(deg_summary_deseq2) # FindAllMarkers - DESeq2
print(deg_summary_mast) # FindAllMarkers - MAST

# Mostrar los resultados de DEFacGroup por clúster
print(deg_counts_df_filtered_1) # DeFacGroup Foldchange > 1.5

# Mostrar las tablas de resumen para FindAllMarkers
print(head(summary_over_df_wilcox))
print(dim(summary_over_df_wilcox))

print(head(summary_over_df_t_test))
print(dim(summary_over_df_t_test))

print(head(summary_over_df_deseq2))
print(dim(summary_over_df_deseq2))

print(head(summary_over_df_mast))
print(dim(summary_over_df_mast))

#######################################

### 3 ################### MÉTRICAS #########################


# Cargar las librerías necesarias
library(caret)

# Función para calcular precisión, recall, F1-score, accuracy y FDR
calcular_metricas <- function(esperados, encontrados, all_genes) {
  tp <- length(intersect(esperados, encontrados)) # True Positives
  fp <- length(setdiff(encontrados, esperados))  # False Positives
  fn <- length(setdiff(esperados, encontrados)) # False Negatives
  tn <- length(setdiff(setdiff(all_genes, esperados), encontrados)) # True Negatives
  
  precision <- ifelse(tp + fp == 0, 0, tp / (tp + fp))
  recall <- ifelse(tp + fn == 0, 0, tp / (tp + fn))
  f1_score <- ifelse(precision + recall == 0, 0, 2 * (precision * recall) / (precision + recall))
  accuracy <- ifelse(tp + fp + fn + tn == 0, 0, (tp + tn) / (tp + fp + fn + tn))
  fdr <- ifelse(tp + fp == 0, 0, fp / (tp + fp))
  
  return(list(precision = precision, recall = recall, f1_score = f1_score, accuracy = accuracy, fdr = fdr))
}

# Listas para almacenar los resultados
metricas_findallmarkers <- list()

# Todos los genes
all_genes <- rownames(sim.groups)

# Comparación para DEFacGroup > 1.5

# DEFacGroup esperado
esperados <- summary_over_df_filtered_1$Gene

# Cálculo de métricas para diferentes métodos de prueba
encontrados_wilcox <- summary_over_df_wilcox$Gene
metricas_findallmarkers[["Wilcoxon"]] <- calcular_metricas(esperados, encontrados_wilcox, all_genes)

encontrados_t_test <- summary_over_df_t_test$Gene
metricas_findallmarkers[["T-Test"]] <- calcular_metricas(esperados, encontrados_t_test, all_genes)
encontrados_deseq2 <- summary_over_df_deseq2$Gene
metricas_findallmarkers[["DESeq2"]] <- calcular_metricas(esperados, encontrados_deseq2, all_genes)

encontrados_mast <- summary_over_df_mast$Gene
metricas_findallmarkers[["MAST"]] <- calcular_metricas(esperados, encontrados_mast, all_genes)

# Convertir listas a dataframes para facilitar la visualización
metricas_findallmarkers_df <- do.call(rbind, lapply(metricas_findallmarkers, function(x) unlist(x)))
metricas_findallmarkers_df <- data.frame(Test = names(metricas_findallmarkers), metricas_findallmarkers_df)

# Mostrar las métricas para FindAllMarkers con diferentes métodos de prueba
print("Métricas para FindAllMarkers con diferentes métodos de prueba:")
print(metricas_findallmarkers_df)



############# AUROC ######


library(pROC)
library(PRROC)

library(dplyr)

calcular_auroc <- function(esperados, all_genes, scores) {
  # Crear un vector binario para los genes esperados
  labels <- as.numeric(all_genes %in% esperados)
  
  # Filtrar valores NA e infinitos
  valid_indices <- !is.na(scores) & is.finite(scores)
  labels <- labels[valid_indices]
  scores <- scores[valid_indices]
  
  # Calcular AUROC
  roc_obj <- roc(labels, scores)
  auroc <- auc(roc_obj)
  
  return(auroc)
}

all_genes <- rownames(sim.groups)

# Comparación para DEFacGroup > 1.5
# DEFacGroup esperado
esperados <- summary_over_df_filtered_1$Gene

# Calcular AUROC para Wilcoxon
scores_wilcox <- -log10(comparison_result_wilcox$p_val_adj)
scores_wilcox <- scores_wilcox[match(all_genes, rownames(comparison_result_wilcox))]
auroc_wilcox <- calcular_auroc(esperados, all_genes, scores_wilcox)

# Calcular AUROC para T-Test
scores_t_test <- -log10(comparison_result_t_test$p_val_adj)
scores_t_test <- scores_t_test[match(all_genes, rownames(comparison_result_t_test))]
auroc_t_test <- calcular_auroc(esperados, all_genes, scores_t_test)

# Calcular AUROC para DESeq2
scores_deseq2 <- -log10(comparison_result_deseq2$p_val_adj)
scores_deseq2 <- scores_deseq2[match(all_genes, rownames(comparison_result_deseq2))]
auroc_deseq2 <- calcular_auroc(esperados, all_genes, scores_deseq2)

# Calcular AUROC para MAST
scores_mast <- -log10(comparison_result_mast$p_val_adj)
scores_mast <- scores_mast[match(all_genes, rownames(comparison_result_mast))]
auroc_mast <- calcular_auroc(esperados, all_genes, scores_mast)

# AUROC para cada método
print("AUROC para Wilcoxon:")
print(auroc_wilcox)

print("AUROC para T-Test:")
print(auroc_t_test)

print("AUROC para DESeq2:")
print(auroc_deseq2)

print("AUROC para MAST:")
print(auroc_mast)

############### PRAUC ######


calcular_prauc <- function(esperados, all_genes, scores) {
  # Crear un vector binario para los genes esperados
  labels <- as.numeric(all_genes %in% esperados)
  
  # Filtrar valores NA e infinitos
  valid_indices <- !is.na(scores) & is.finite(scores)
  labels <- labels[valid_indices]
  scores <- scores[valid_indices]
  
  # Calcular PRAUC
  pr_obj <- pr.curve(scores.class0 = scores[labels == 0],
                     scores.class1 = scores[labels == 1],
                     curve = TRUE)
  prauc <- pr_obj$auc.integral
  
  return(prauc)
}


# Comparación para DEFacGroup > 1.5
# DEFacGroup esperado
esperados <- summary_over_df_filtered_1$Gene

# Calcular PRAUC para Wilcoxon
scores_wilcox <- -log10(comparison_result_wilcox$p_val_adj)
scores_wilcox <- scores_wilcox[match(all_genes, rownames(comparison_result_wilcox))]
prauc_wilcox <- calcular_prauc(esperados, all_genes, scores_wilcox)

# Calcular PRAUC para T-Test
scores_t_test <- -log10(comparison_result_t_test$p_val_adj)
scores_t_test <- scores_t_test[match(all_genes, rownames(comparison_result_t_test))]
prauc_t_test <- calcular_prauc(esperados, all_genes, scores_t_test)

# Calcular PRAUC para DESeq2
scores_deseq2 <- -log10(comparison_result_deseq2$p_val_adj)
scores_deseq2 <- scores_deseq2[match(all_genes, rownames(comparison_result_deseq2))]
prauc_deseq2 <- calcular_prauc(esperados, all_genes, scores_deseq2)

# Calcular PRAUC para MAST
scores_mast <- -log10(comparison_result_mast$p_val_adj)
scores_mast <- scores_mast[match(all_genes, rownames(comparison_result_mast))]
prauc_mast <- calcular_prauc(esperados, all_genes, scores_mast)

#  PRAUC para cada método
print("PRAUC para Wilcoxon:")
print(prauc_wilcox)

print("PRAUC para T-Test:")
print(prauc_t_test)

print("PRAUC para DESeq2:")
print(prauc_deseq2)

print("PRAUC para MAST:")
print(prauc_mast)
