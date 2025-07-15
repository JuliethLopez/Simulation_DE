# Cargar bibliotecas necesarias
library(Seurat)
library(MAST)
library(dplyr)
library(openxlsx)
library("SeuratData")
library(dplyr)


# Función para calcular métricas de comparación de DEGs
calculate_metrics <- function(reference_genes, simulated_genes) {
  common_genes <- intersect(reference_genes, simulated_genes)
  true_positive <- length(common_genes)
  false_positive <- length(setdiff(simulated_genes, reference_genes))
  false_negative <- length(setdiff(reference_genes, simulated_genes))
  
  precision <- true_positive / (true_positive + false_positive)
  recall <- true_positive / (true_positive + false_negative)
  f1_score <- 2 * (precision * recall) / (precision + recall)
  
  return(list(precision = precision, recall = recall, f1_score = f1_score))
}

# 1. Análisis de expresión diferencial en datos de referencia (pbmc3k)

## charge data Ref ####
####. 1. DATOS Referencia pbmc3k #####
pbmc3k <- LoadData("pbmc3k.SeuratData")
pbmc3k <- UpdateSeuratObject(pbmc3k)

# 1 (data)
data_0 <- pbmc3k

# Preprocesamiento de los datos
pbmc3k <- NormalizeData(pbmc3k)           # Normalización de los datos
pbmc3k <- FindVariableFeatures(pbmc3k)    # Identificación de características variables
pbmc3k <- ScaleData(pbmc3k)               # Escalado de los datos
pbmc3k <- RunPCA(pbmc3k)                  # Análisis de componentes principales

# Agrupamiento de las células
pbmc3k <- FindNeighbors(pbmc3k, dims = 1:10)
pbmc3k <- FindClusters(pbmc3k, resolution = 0.5)

# Corre análisis UMAP o t-SNE para visualización
pbmc3k <- RunUMAP(pbmc3k, dims = 1:10)

# Visualización del agrupamiento en UMAP
DimPlot(pbmc3k, reduction = "umap", label = TRUE)


# Identificar genes diferencialmente expresados en pbmc3k
pbmc3k_markers <- FindAllMarkers(pbmc3k, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
# Ver los primeros 5 genes diferencialmente expresados por cluster
head(pbmc3k_markers)
# Categorizar genes en altos, medios y bajos según logFC
high_fc_genes <- pbmc3k_markers %>% filter(avg_log2FC > 1)
medium_fc_genes <- pbmc3k_markers %>% filter(avg_log2FC <= 1 & avg_log2FC >= 0.5)
low_fc_genes <- pbmc3k_markers %>% filter(avg_log2FC < 0.5)

# Guardar los resultados en archivos Excel
write.xlsx(list("High_logFC" = high_fc_genes, "Medium_logFC" = medium_fc_genes, "Low_logFC" = low_fc_genes), 
           file = "DEGs_by_logfc_level.xlsx")

# Obtener los genes diferencialmente expresados en pbmc3k
reference_genes <- pbmc3k_markers$gene

# 2. Análisis de expresión diferencial en los datos simulados (simulation_2.rds)

# Cargar los datos simulados
# Cargar el archivo simulation_2.rds
simulation_2 <- readRDS("simulation_2.rds")
# Verifica qué tipo de objeto se ha cargado
class(simulation_2)
# Si es una lista o data frame, verifica su contenido
str(simulation_2)
matrix_counts <- simulation_2@assays@data@listData[["counts"]]
# Usar el objeto simulation_2 como matriz de cuentas
simulation_2_seurat <- CreateSeuratObject(counts = matrix_counts)

# Normalizar y escalar los datos
simulation_2_seurat <- NormalizeData(simulation_2_seurat)
simulation_2_seurat <- FindVariableFeatures(simulation_2_seurat)
simulation_2_seurat <- ScaleData(simulation_2_seurat)

# Realizar clustering en la simulación
simulation_2_seurat <- RunPCA(simulation_2_seurat)
simulation_2_seurat <- FindNeighbors(simulation_2_seurat, dims = 1:10)
simulation_2_seurat <- FindClusters(simulation_2_seurat, resolution = 0.5)
simulation_2_seurat <- RunUMAP(simulation_2_seurat, dims = 1:10)

# Realizar análisis de expresión diferencial en la simulación usando MAST
simulation_markers <- FindAllMarkers(simulation_2_seurat, test.use = "MAST", only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# Obtener los genes diferencialmente expresados en la simulación
simulated_genes <- simulation_markers$gene

# 3. Comparar los genes diferencialmente expresados entre pbmc3k y la simulación

# Calcular métricas de comparación
metrics <- calculate_metrics(reference_genes, simulated_genes)

# Mostrar los resultados de las métricas
print(paste("Precisión: ", round(metrics$precision, 3)))
print(paste("Recall: ", round(metrics$recall, 3)))
print(paste("F1-score: ", round(metrics$f1_score, 3)))

# Guardar las métricas en un archivo de texto
write.table(metrics, file = "comparison_metrics.txt", col.names = NA, quote = FALSE)


# Imprime los genes en cada conjunto
print("Genes de referencia:")
print(reference_genes)

print("Genes simulados:")
print(simulated_genes)

####  CONTAR CUANTOS GENES HAY EN AMBOS GRUPOS###############

# Clasificar genes según niveles de expresión
df_genes <- pbmc3k@assays$RNA@data  # o el objeto que estés usando
df_genes <- data.frame(gene = rownames(df_genes), DEFacGroup = apply(df_genes, 1, mean))  # Ejemplo para obtener logFC

df_genes <- df_genes %>%
  mutate(logfc_level = case_when(
    DEFacGroup >= 1.2 & DEFacGroup < 2 ~ "low",
    DEFacGroup >= 2 & DEFacGroup < 3 ~ "medium",
    DEFacGroup >= 3 ~ "high",
    TRUE ~ "not_differential"  # Categoría para genes que no cumplen condiciones
  ))

# Contar genes en cada categoría en el conjunto de referencia
count_reference <- df_genes %>%
  group_by(logfc_level) %>%
  summarise(count = n()) %>%
  deframe()  # Convierte a un vector



# Ver la clase del objeto
class(simulation_2_seurat)
# Acceder a los datos de conteo
df_simulated_genes <- GetAssayData(simulation_2_seurat, slot = "counts")

# Convertir la matriz a dataframe
df_simulated_genes <- as.data.frame(as.matrix(df_simulated_genes))

# Calcular el promedio de expresión para cada gen
df_simulated_genes$DEFacGroup <- rowMeans(df_simulated_genes, na.rm = TRUE)

# Clasificar genes en el conjunto simulado
df_simulated_genes$logfc_level <- case_when(
  df_simulated_genes$DEFacGroup >= 1.2 & df_simulated_genes$DEFacGroup < 2 ~ "low",
  df_simulated_genes$DEFacGroup >= 2 & df_simulated_genes$DEFacGroup < 3 ~ "medium",
  df_simulated_genes$DEFacGroup >= 3 ~ "high",
  TRUE ~ "not_differential"
)

# Contar genes en cada categoría en el conjunto simulado
count_simulated <- df_simulated_genes %>%
  group_by(logfc_level) %>%
  summarise(count = n()) %>%
  ungroup()  # Asegúrate de desagrupar después de resumir

# Mostrar conteos
print(count_simulated)

# Asegúrate de que el conteo simulado esté en el formato correcto
count_simulated_vector <- setNames(as.numeric(count_simulated$count), count_simulated$logfc_level)

# Crear dataframe de métricas
metrics <- data.frame(
  Category = names(count_reference),
  Reference_Count = count_reference,
  Simulated_Count = count_simulated_vector[names(count_reference)],  # Asegúrate de que las categorías coincidan
  Difference = count_simulated_vector[names(count_reference)] - count_reference
)

# Mostrar métricas
print(metrics)

# Guardar las métricas en un archivo de texto
write.table(metrics, file = "comparison_metrics.txt", row.names = FALSE, quote = FALSE)


######CALCULAR METRICAS ESTADISTICAS################
print(count_simulated)

library(dplyr)

# Definir las categorías
categories <- c("high", "medium", "low")

# Inicializar un data frame para almacenar las métricas
metrics <- data.frame(Category = categories, Precision = NA, Recall = NA, F1_Score = NA, Accuracy = NA, FDR = NA)

# Calcular métricas para cada categoría
# Calcular métricas para cada categoría
for (category in categories) {
  
  # Obtener los conteos de referencia
  reference_count <- count_reference[category]
  
  # Obtener el conteo simulado usando dplyr
  simulated_count <- count_simulated %>%
    filter(logfc_level == category) %>%
    pull(count)
  
  # Calcular TP, FP, FN, TN
  TP <- min(simulated_count, reference_count)  # Verdaderos positivos
  FP <- max(0, simulated_count - reference_count)  # Falsos positivos, no pueden ser negativos
  FN <- max(0, reference_count - simulated_count)  # Falsos negativos, no pueden ser negativos
  TN <- sum(count_reference[categories[categories != category]])  # Verdadero negativo
  
  # Calcular las métricas
  precision <- ifelse((TP + FP) > 0, TP / (TP + FP), 0)
  recall <- ifelse((TP + FN) > 0, TP / (TP + FN), 0)
  f1_score <- ifelse((precision + recall) > 0, 2 * (precision * recall) / (precision + recall), 0)
  accuracy <- (TP + TN) / (TP + TN + FP + FN)
  fdr <- ifelse((TP + FP) > 0, FP / (FP + TP), 0)
  
  # Almacenar los resultados
  metrics[metrics$Category == category, "Precision"] <- precision
  metrics[metrics$Category == category, "Recall"] <- recall
  metrics[metrics$Category == category, "F1_Score"] <- f1_score
  metrics[metrics$Category == category, "Accuracy"] <- accuracy
  metrics[metrics$Category == category, "FDR"] <- fdr
}

# Mostrar las métricas
print(metrics)




