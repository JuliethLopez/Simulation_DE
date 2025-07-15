library(SingleCellExperiment)
library(splatter)
library(Seurat)
library(openxlsx)

# 1. Configuración inicial y estimación de parámetros a partir de datos reales
params <- splatEstimate(as.matrix(pbmc3k@assays$RNA$counts))

# 2. Función para realizar una simulación con parámetros específicos
run_simulation <- function(params, group_prob, de_prob, de_facLoc, de_facScale, batchCells, nGenes, dropout_mid, dropout_shape) {
  
  params <- setParams(params, dropout.type = "none", dropout.mid = dropout_mid, dropout.shape = dropout_shape)
  
  sim <- splatSimulate(
    params,
    group.prob = group_prob,
    de.prob = de_prob,
    de.facLoc = de_facLoc,
    de.facScale = de_facScale,
    batchCells = batchCells,
    nGenes = nGenes,
    method = "groups",
    verbose = TRUE
  )
  
  return(sim)
}

# 3. Simulación de referencia con parámetros iniciales
sim1 <- run_simulation(
  params,
  group_prob = c(0.45, 0.23, 0.15, 0.10, 0.04, 0.03),
  de_prob = c(0.10, 0.09, 0.07, 0.06, 0.04, 0.04),
  de_facLoc = c(0.4, 1, 0.8, 0.8, 0.5, 1),
  de_facScale = c(0.9, 0.9, 0.7, 0.8, 0.8, 0.9),
  batchCells = 3000,
  nGenes = 15000,
  dropout_mid = -0.15,
  dropout_shape = -1.06
)

# 4. Simulación alternativa con un mayor número de células y genes
sim2 <- run_simulation(
  params,
  group_prob = c(0.50, 0.20, 0.15, 0.10, 0.03, 0.02),
  de_prob = c(0.12, 0.08, 0.07, 0.06, 0.04, 0.03),
  de_facLoc = c(0.5, 0.9, 0.9, 0.7, 0.6, 1),
  de_facScale = c(0.95, 0.85, 0.75, 0.85, 0.9, 0.95),
  batchCells = 3000,
  nGenes = 14000,
  dropout_mid = -0.2,
  dropout_shape = -1.09
)

sim.groups <- sim2

# 5. Visualización y comparación de las simulaciones

# Preprocesamiento de sim1
sim1_counts <- sim1@assays@data@listData[["counts"]]
data_sim1 <- CreateSeuratObject(counts = sim1_counts)
data_sim1 <- NormalizeData(data_sim1)
data_sim1 <- FindVariableFeatures(data_sim1)
data_sim1 <- ScaleData(data_sim1)
data_sim1 <- RunPCA(data_sim1, features = VariableFeatures(data_sim1))
data_sim1 <- FindNeighbors(data_sim1, dims = 1:10)
data_sim1 <- FindClusters(data_sim1, resolution = 0.8)
data_sim1 <- RunUMAP(data_sim1, dims = 1:10)
DimPlot(data_sim1, reduction = "umap") + ggtitle("Simulación 1 (UMAP)")

# Preprocesamiento de sim2
sim2_counts <- sim2@assays@data@listData[["counts"]]
data_sim2 <- CreateSeuratObject(counts = sim2_counts)
data_sim2 <- NormalizeData(data_sim2)
data_sim2 <- FindVariableFeatures(data_sim2)
data_sim2 <- ScaleData(data_sim2)
data_sim2 <- RunPCA(data_sim2, features = VariableFeatures(data_sim2))
data_sim2 <- FindNeighbors(data_sim2, dims = 1:10)
data_sim2 <- FindClusters(data_sim2, resolution = 0.8)
data_sim2 <- RunUMAP(data_sim2, dims = 1:10)
DimPlot(data_sim2, reduction = "umap") + ggtitle("Simulación 2 (UMAP)")

# 6. Guardar los resultados de la simulación
wb_sim <- createWorkbook()

# Guardar datos de la primera simulación
addWorksheet(wb_sim, "DEGs Sim1")
writeData(wb_sim, sheet = "DEGs Sim1", FindAllMarkers(data_sim1, only.pos = TRUE))

# Guardar datos de la segunda simulación
addWorksheet(wb_sim, "DEGs Sim2")
writeData(wb_sim, sheet = "DEGs Sim2", FindAllMarkers(data_sim2, only.pos = TRUE))

saveWorkbook(wb_sim, file = "Simulations_DEGs.xlsx", overwrite = TRUE)
