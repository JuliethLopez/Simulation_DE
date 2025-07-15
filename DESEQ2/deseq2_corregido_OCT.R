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
library(pROC)
library(reshape2)
library(data.table)
install.packages("C:/Users/Usuario/OneDrive/Escritorio/proyectoestadistica/mio/conflicted-main", repos = NULL, type = "source")
library("tidyverse")
library("SeuratData") #data reference PBMC
# Instalar gplots si no está instalado
if (!requireNamespace("gplots", quietly = TRUE)) {
  install.packages("gplots")
}

library(gplots)

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

# Convertir datos en formato de DESeq2, si no lo están
coldata <- data.frame(Cluster = factor(colData(sim.groups)$Group, levels = c("Group1", "Group2", "Group3", "Group4", "Group5", "Group6")))
rownames(coldata) <- colnames(sim.groups)

#### SEGUNDA PARTE: ANÁLISIS DE EXPRESIÓN DIFERENCIAL CON DESEQ2 ####

# Crear el conjunto de datos DESeq2 desde la matriz de conteo
deseq2_matrix <- DESeqDataSetFromMatrix(
  countData = assay(sim.groups, "counts"), 
  colData = coldata, 
  design = ~ 0 + Cluster  # Modelo con 6 grupos/clusters
)

DE_metacells_deseq2_original <- list()

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
  DE_metacells_deseq2_original[[level]] <- results
}

# Combinar los resultados en un solo dataframe
DE_metacells_deseq2_groups <- bind_rows(DE_metacells_deseq2_original)
colnames(DE_metacells_deseq2_groups, c("log2FoldChange", "padj"), c("DEFacGroup", "p_val_adj"))

# Guardar los resultados en un archivo CSV
write.csv(DE_metacells_deseq2_groups, file = paste0(outdir, "/DE_metacells_deseq2_original.csv"), row.names = TRUE)




#### TERCERA PARTE: FUNCIONES DE ANÁLISIS Y CÁLCULO DE MÉTRICAS AJUSTADAS A DESEQ2 ####

# Definir todos los genes de referencia
all_genes <- unique(DE_metacells_deseq2_groups$Gene)

# Función para analizar expresión diferencial (ajustada a DESeq2)
DE_analysis_deseq2 <- function(deseq_results, max_p_val_adj = 0.05, min_pct_difference = 0.05) {
  summary_over_df <- deseq_results %>%
    mutate(logfc_level = case_when(
      DEFacGroup >= 1.2 & DEFacGroup < 2 ~ "low",
      DEFacGroup >= 2 & DEFacGroup < 3 ~ "medium",
      DEFacGroup >= 3 ~ "high",
      TRUE ~ "otros"
    ))
  
  # Filtrar según el valor ajustado de p
  summary_over_df <- summary_over_df %>%
    filter(p_val_adj < max_p_val_adj)
  
  return(summary_over_df)
}

# Filtrar y procesar los resultados de DESeq2
deseq_results_filtered <- DE_analysis_deseq2(DE_metacells_deseq2_groups)

# Función para calcular el AUROC
calcular_auroc <- function(markers_expected, markers_found) {
  labels_total = c()
  scores_total = c()
  
  for (level in levels(markers_expected$cluster)) {
    founded <- subset(markers_found, cluster == level)
    expected <- subset(markers_expected, cluster == level)
    
    labels <- as.numeric(expected$Gene %in% founded$Gene)
    scores <- expected %>% left_join(founded, by = "Gene") %>%
      select(p_val_adj)
    scores <- -log10(scores[, "p_val_adj"])
    
    no_valid_indices <- is.na(scores) | is.infinite(scores)
    labels <- labels[!no_valid_indices]
    scores <- scores[!no_valid_indices]
    
    labels_total <- c(labels_total, labels)
    scores_total <- c(scores_total, scores)
  }
  
  if (length(unique(labels_total)) > 1) {
    roc_obj <- roc(labels_total, scores_total)
    auroc <- auc(roc_obj)
  } else {
    auroc <- 0
  }
  
  return(auroc)
}

# Función para calcular métricas (precisión, recall, f1, etc.)
calculate_metrics <- function(real, predicted, all_genes) {
  predicted_gene_cluster <- predicted %>% select(c(Gene, cluster))
  TP <- nrow(dplyr::intersect(real, predicted_gene_cluster))
  FP <- nrow(dplyr::setdiff(predicted_gene_cluster, real))
  FN <- nrow(dplyr::setdiff(real, predicted_gene_cluster))
  TN <- nrow(dplyr::setdiff(dplyr::setdiff(all_genes, real$Gene), predicted_gene_cluster$Gene))
  
  precision <- ifelse((TP + FP) > 0, TP / (TP + FP), 0)
  recall <- ifelse((TP + FN) > 0, TP / (TP + FN), 0)
  f1_score <- ifelse((precision + recall) > 0, 2 * (precision * recall) / (precision + recall), 0)
  accuracy <- (TP + TN) / (TP + TN + FP + FN)
  fdr <- ifelse((FP + TP) > 0, FP / (FP + TP), 0)
  auroc <- calcular_auroc(real, predicted)
  
  return(c(precision, recall, f1_score, accuracy, fdr, auroc))
}

# Función para aplicar el cálculo de métricas en diferentes niveles
DE_metrics <- function(expected_df, founded_filtered_df, all_genes) {
  expected_df$cluster <- as.factor(expected_df$cluster)
  
  # Filtrar según niveles de logFC
  expected_low <- expected_df %>% subset(logfc_level == "low") %>% select(c(Gene, cluster))
  expected_medium <- expected_df %>% subset(logfc_level == "medium") %>% select(c(Gene, cluster))
  expected_high <- expected_df %>% subset(logfc_level == "high") %>% select(c(Gene, cluster))
  expected_all <- expected_df %>% select(c(Gene, cluster))
  
  predicted_low <- founded_filtered_df %>% filter(logfc_level == "low") %>% select(c(Gene, cluster, p_val_adj))
  predicted_medium <- founded_filtered_df %>% filter(logfc_level == "medium") %>% select(c(Gene, cluster, p_val_adj))
  predicted_high <- founded_filtered_df %>% filter(logfc_level == "high") %>% select(c(Gene, cluster, p_val_adj))
  predicted_all <- founded_filtered_df %>% select(c(Gene, cluster, p_val_adj))
  
  metrics_low <- calculate_metrics(expected_low, predicted_low, all_genes)
  metrics_medium <- calculate_metrics(expected_medium, predicted_medium, all_genes)
  metrics_high <- calculate_metrics(expected_high, predicted_high, all_genes)
  metrics_all <- calculate_metrics(expected_all, predicted_all, all_genes)
  
  metrics_table <- data.frame(
    Level = c("low", "medium", "high", "all"),
    Precision = c(metrics_low[1], metrics_medium[1], metrics_high[1], metrics_all[1]),
    Recall = c(metrics_low[2], metrics_medium[2], metrics_high[2], metrics_all[2]),
    F1_Score = c(metrics_low[3], metrics_medium[3], metrics_high[3], metrics_all[3]),
    Accuracy = c(metrics_low[4], metrics_medium[4], metrics_high[4], metrics_all[4]),
    FDR = c(metrics_low[5], metrics_medium[5], metrics_high[5], metrics_all[5]),
    AUroc = c(metrics_low[6], metrics_medium[6], metrics_high[6], metrics_all[6])
  )
  
  results_list <- list(DEGs = expected_df, Metrics = metrics_table)
  return(results_list)
}

# Calcular las métricas basadas en los resultados filtrados de DESeq2
metrics_results <- DE_metrics(deseq_results_filtered, deseq_results_filtered, all_genes)

# Guardar las métricas en un archivo CSV
write.csv(metrics_results$Metrics, file = paste0(outdir, "/metrics_deseq2.csv"), row.names = TRUE)


### HEATMAP ####
Metrics_heatmap <- function(results_list, outdir, text = NULL) {
  for (test in names(results_list)) {
    # Transformar la matriz en formato largo (long format)
    df <- reshape2::melt(results_list[[test]]$Metrics)
    colnames(df) <- c("Level", "Metric", "Value")
    
    df$Value[df$Metric == "FDR"] <- 1 - df$Value[df$Metric == "FDR"]
    df$Value_discrete <- ifelse(df$Value < 0.3, 1, ifelse(df$Value >= 0.7, 3, df$Value))
    df$Value_discrete <- ifelse(0.3 <= df$Value & df$Value < 0.7, 2, df$Value_discrete)
    df$Value_discrete <- factor(df$Value_discrete, levels = c("1", "2", "3"))
    df$Test_level <- paste(df$Metric, df$Level, sep = "_")
    
    # Crear el gráfico de heatmap
    heatmap_test_metrics <- ggplot(df, aes(x = Metric, y = Test_level, fill = Value_discrete)) +
      geom_tile() +
      ggtitle("Heatmap test metrics") +
      theme(plot.title = element_text(hjust = 0.5),
            text = element_text(size = 13),
            axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
      coord_fixed() +
      scale_fill_manual(breaks = levels(df$Value_discrete),
                        values = c("#FF99BF", "#FFFF99", "#A5EDFF"),
                        labels = c('Poor', 'Intermediate', 'Good'))
    
    print(heatmap_test_metrics)
    ggsave(filename = paste0(outdir, "/Heatmap_test_metric_", text, "_", test, ".jpg"), plot = heatmap_test_metrics)
  }
}

# Crear el heatmap basado en las métricas calculadas
Metrics_heatmap(list(deseq2 = metrics_results), outdir, text = "DESeq2")
