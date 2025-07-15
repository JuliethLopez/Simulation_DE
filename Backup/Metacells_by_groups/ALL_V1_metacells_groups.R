## Set environment - Change for your own directory
basedir <- "C:/Users/Administrator/Documents/Jaguar/Simulation_DE/Metacells_by_groups"
setwd("C:/Users/Administrator/Documents/Jaguar/Simulation_DE")
outdir <- "C:/Users/Administrator/Documents/Jaguar/Simulation_DE/Metacells_by_groups/figures_metacells"

#source("C:/Users/Administrator/Documents/Jaguar/Simulation_DE/ALL_V1.R")

# Libraries
library(SuperCell)
library(Seurat)
# If you have Seurat V5 installed, specify that you want to analyze Seurat V4 objects
#if(packageVersion("Seurat") >= 5) {options(Seurat.object.assay.version = "v4"); print("you are using seurat v5 with assay option v4")}
library(SeuratData)
library(tidyverse)
library(reticulate)
library(spatstat.utils)
library(MetacellAnalysisToolkit)

################################################################################
## Data
################################################################################

# Data
# Here is called data pre-processed in ALL_V1
data <- readRDS("simulation_2.rds")
data_1_cl <- as.Seurat(data, counts = "counts", data = NULL)
data_1_cl #15000 features across 4000 samples

# Plot variable features 
data_1_cl <- NormalizeData(data_1_cl, normalization.method = "LogNormalize", scale.factor = 10000)
data_1_cl <- FindVariableFeatures(data_1_cl, selection.method = "vst", nfeatures = 2000)
data_1_cl <- ScaleData(data_1_cl)

hvg <- VariableFeatures(data_1_cl)
length(hvg)
plot1 <- VariableFeaturePlot(data_1_cl)
LabelPoints(plot = plot1, points = hvg[1:20], repel = TRUE)
ggsave(paste0(outdir,"/plot_variable_features_single_cell_data.png"),width = 10,height = 8, dpi = 300)

# Plot PCA (2D representation of scRNA-seq data) colored by cell line
DimPlot(data_1_cl, reduction = "pca", group.by = "Cluster",label=T, raster=F) + ggtitle("PCA_single_cell_by_clusters")
ggsave(paste0(outdir,"/pca_single_cell_data.png"),width = 10,height = 8, dpi = 300)

# Plot UMAP
DimPlot(data_1_cl, reduction = "umap", group.by = "Cluster",label=T, raster=F)+ ggtitle("UMAMP_single_cell_by_clusters")
ggsave(paste0(outdir,"/umap_single_cell_data.png"),width = 10,height = 8, dpi = 300)

################################################################################
## Metacells construction
################################################################################

# Set idents
Idents(data_1_cl)<- "Cluster"
Cluster <- data_1_cl$Cluster

# Get gene expression matrix
GE = GetAssayData(data_1_cl, layer = "counts")
dim(GE) #15000  4000
head(rownames(GE))
head(colnames(GE))

gamma = 20 # Graining level
k.knn = 5 # Number of nearest neighbors to build KNN network

# Compute metacells per groups using SuperCell package
MC <- SCimplify(
  X = GE, # single-cell log-normalized gene expression data
  genes.use = hvg, #highly variable genes
  cell.annotation	= Cluster, #vector with groups of the cells
  #n.var.genes = 2000,
  gamma = gamma,
  k.knn = k.knn
)

length(unique(MC$membership)) #299

#network plot
png(paste0(outdir,"/network_metacell_data_cluster.png"))
supercell_plot(
  MC$graph.supercells, #networks
  #group = MC$Cluster, 
  seed = 1, 
  alpha = -pi/2,
  main  = "Network_metacells"
)
dev.off()

# Compute gene expression of metacells by simply averaging log-normalized gene expression within each metacell
# Alternatively, counts can be averaged (summed up) followed by a lognormalization step (this approach is used in the MetaCell and SEACell algorithms)

if(1){
  MC.counts <- supercell_GE(
    ge = GE,
    groups = MC$membership,
    mode =  "sum"
  )
  MC.ge <- Seurat::LogNormalize(MC.counts, verbose = FALSE)
}

head(rownames(MC.ge))
colnames(MC.ge) = paste0("MetaCell",1:ncol(MC.ge))
head(colnames(MC.ge))
head(MC.counts)
head(MC.ge)

# Annotate metacells to Cells Cluster
MC$Cluster <- supercell_assign(
  cluster = data_1_cl$Cluster,          # single-cell assignment to Cells Cluster 
  supercell_membership = MC$membership,  # single-cell assignment to metacells
  method = "absolute" # available methods are c("jaccard", "relative", "absolute"), function's help() for explanation
)

#network plot metacells
png(paste0(outdir,"/network_metacell_data_by_groups.png"))
supercell_plot(
  MC$graph.supercells, #networks
  group = MC$Cluster,
  seed = 1, 
  alpha = -pi/2,
  main  = "Metacells colored by groups assignment"
)
dev.off()

#network plot single cells
png(paste0(outdir,"/network_single_cell_data_by_groups.png"))
supercell_plot(
  MC$graph.singlecell, 
  group = data_1_cl$Cluster, 
  do.frames = FALSE,
  lay.method = "components",
  seed = 1, 
  alpha = -pi/2,
  main  = "Single cells colored by groups assignment"
)
dev.off()

################################################################################
## Quality Control
################################################################################

#purity
# Compute purity of metacells as :
#  * a proportion of the most abundant cell type withing metacells (`method = `"max_proportion)
#  * an entropy of cell type within metacells (`method = "entropy"`)
method_purity <- c("max_proportion", "entropy")[1]
MC$purity <- supercell_purity(
  clusters = data_1_cl$Cluster,
  supercell_membership = MC$membership, 
  method = method_purity
)

# Metacell purity distribution
summary(MC$purity)
png(paste0(outdir,"/histogram_purity_metacell_data.png"))
hist(MC$purity, main = paste0("Purity of metacells (", method_purity,")"))
dev.off()

png(file=paste0(outdir,"/boxplot_purity.png"),
    width=600, height=350)
boxplot(MC$purity, main = paste0("Purity of metacells (", method_purity,")"))
dev.off() # a function call to save the file

png(file=paste0(outdir,"/boxplot_Clusters_purity.png"),
    width=600, height=350)
boxplot(MC$purity ~ MC$Cluster, main = paste0("Purity of metacells \nin terms of cell groups (", method_purity,")"))
dev.off() # a function call to save the file

#size metacells: Having a homogeneous metacell size distribution is ideal for 
#downstream analyses, since larger metacells will express more genes, which 
#could confound analyses. When heterogeneous size distributions are obtained 
#we recommend weighted downstream analyses 
png(file=paste0(outdir,"/histplot_size_distribution.png"),
    width=600, height=350)
hist(MC$supercell_size, main = "Size distribution", xlab = "Size", breaks = 50)
dev.off() # a function call to save the file

library(vioplot)
png(file=paste0(outdir,"/violinplot_size_distribution.png"),
    width=600, height=350)
vioplot(MC$supercell_size, main = "Size distribution", xlab = "Size")
dev.off() # a function call to save the file

# Other QC metrics
# https://gfellerlab.github.io/MetacellAnalysisTutorial/downstream-analysis.html#weighted-analysis

# compactness
# We compute these metrics in a diffusion map obtained from the pca
membership_df <- MC$membership
diffusion_comp <- get_diffusion_comp(sc.obj = MC, dims = 1:30)

## Computing diffusion maps ...

mc_data$compactness <- mc_compactness(cell.membership = membership_df,
                                      sc.obj = sc_data,
                                      sc.reduction = diffusion_comp)
qc_boxplot(mc.obj = mc_data, qc.metrics = "compactness")

################################################################################
## Transform data for use in seurat
################################################################################

# Standar downstream analysis
# pass supercell to a seurat object.
# Note: since metacells have different size (consist of different number of cells),
# we apply sample-weighted algorithms at most af the steps of the downstream 
# analyses. Thus, when coercing SuperCell to Seurat, we replaced PCA, scaling and 
# kNN graph of Seurat object with those obtained applying sample-weighted version 
# of PCA, scaling or SuperCell graph (i.e., metacell network), respectively. 
# If you then again apply `RunPCA`, `ScaleData`, or `FindNeighbors`, the result 
# will be rewritten, but you will be able to access them with `Embeddings(m.seurat,
# reduction = "pca_weigted")`, `m.seurat@assays$RNA@misc[["scale.data.weighted"]]`, 
# or `m.seurat@graphs$RNA_super_cells`, respectively

MC.seurat <- supercell_2_Seurat(
  SC.GE = MC.counts,
  SC = MC,
  fields = c("Cluster", "purity", "supercell_size"), # elements of MC to save as metacell metadata 
  var.genes = MC$genes.use,
  N.comp = 10
)

MC.seurat.lognorm <- supercell_2_Seurat(
  SC.GE = MC.ge,
  SC = MC,
  fields = c("Cluster", "purity", "supercell_size"),
  var.genes = MC$genes.use,
  N.comp = 10
)

#It sims it is not normalizing the data although it says is done https://github.com/GfellerLab/SIB_workshop/blob/main/workbooks/Workbook_1__cancer_cell_lines.md
head(MC.seurat@assays$RNA$counts,5)
#Gene1  .  3  4 .  9 
head(MC.seurat@assays$RNA$data,5)
#Gene1  .  3  4 .  9
head(MC.seurat@assays$RNA$scale.data,5)
#Gene1 -0.5640929  0.7585165  1.1993864

head(MC.seurat.lognorm@assays$RNA$counts,5)
#Gene1 .         0.2323282 0.19951878
head(MC.seurat.lognorm@assays$RNA$data,5)
#Gene1 .         0.2323282 0.19951878 
head(MC.seurat.lognorm@assays$RNA$scale.data,5)
#Gene1 -0.64866307  0.6260241  0.44601286

#pca and pca_weighted are exactly the same in this point
any(!c(MC.seurat@reductions$pca@cell.embeddings == MC.seurat@reductions$pca_weighted@cell.embeddings))

################################################################################
## Dimensional reduction un-weighted
################################################################################

## Check dimensional reduction

#as all are the same is enought saving only one pca normal for metacells
Idents(MC.seurat) <- "Cluster"
levels(MC.seurat) <- sort(levels(MC.seurat))
DimPlot(MC.seurat, reduction = "pca_seurat")
DimPlot(MC.seurat, reduction = "pca_weighted")
DimPlot(MC.seurat, reduction = "pca") + ggtitle("pca_Clusters_metacells_data")
ggsave(paste0(outdir,"/pca_Clusters_metacells_data.png"),width = 10,height = 8, dpi = 300)

#saving Umap over normal pca, because weighted is not really calculates with supercell_2_seurat as indicated here https://github.com/GfellerLab/SuperCell/blob/master/vignettes/a_SuperCell.Rmd
MC.seurat <- RunUMAP(MC.seurat, reduction = "pca", dims = 1:10)
DimPlot(MC.seurat, reduction = "umap") + ggtitle("umap_Clusters_metacells_data")
ggsave(paste0(outdir,"/umap_Clusters_metacells_data.png"),width = 10,height = 8, dpi = 300)

#now saving pca for MC.seurat.lognorm
Idents(MC.seurat.lognorm) <- "Cluster"
levels(MC.seurat.lognorm) <- sort(levels(MC.seurat.lognorm))
DimPlot(MC.seurat.lognorm, reduction = "pca")
DimPlot(MC.seurat.lognorm, reduction = "pca_seurat")
DimPlot(MC.seurat.lognorm, reduction = "pca_weighted") + ggtitle("pca_Clusters_metacells_data_lognorm")
ggsave(paste0(outdir,"/pca_Clusters_metacells_data_lognorm.png"),width = 10,height = 8, dpi = 300)

#now saving umap for MC.seurat.lognorm
MC.seurat.lognorm <- RunUMAP(MC.seurat.lognorm, reduction = "pca", dims = 1:10)
DimPlot(MC.seurat.lognorm, reduction = "umap") + ggtitle("umap_Clusters_metacells_data_lognorm")
ggsave(paste0(outdir,"/umap_Clusters_metacells_data_lognorm.png"),width = 10,height = 8, dpi = 300)

## Save object

# Seurat object containing normalized metacells gene expression (? NOT) data as well as 
# the first (N.comp) principal components of PCA performed internally using user 
# defined set of genes (by default the genes used for metacells constructions)
if(packageVersion("Seurat") >= 5) {
  MC.seurat[["RNA"]] <- as(object = MC.seurat[["RNA"]], Class = "Assay")
  MC.seurat.lognorm[["RNA"]] <- as(object = MC.seurat.lognorm[["RNA"]], Class = "Assay")
}
saveRDS(MC.seurat, file = file.path(paste0("MC_gamma_", gamma, "_seurat_object.Rds")))
saveRDS(MC.seurat.lognorm, file = file.path(paste0("MC_gamma_", gamma, "_seurat_lognorm.Rds")))

################################################################################
## Dimensional reduction weighted
################################################################################
# read data
MC.seurat <- readRDS("MC_gamma_20_seurat_object.Rds")
MC.seurat.lognorm <- readRDS("MC_gamma_20_seurat_lognorm.Rds")

#normalization and variable features
MC.seurat <- NormalizeData(MC.seurat, normalization.method = "LogNormalize")
MC.seurat <- FindVariableFeatures(MC.seurat)

#calculate pcs weighted
MC_list <- list(N.SC = ncol(MC.seurat),
                supercell_size = MC.seurat$size)
MC_list$PCA <- SuperCell::supercell_prcomp(
  Matrix::t(GetAssayData(MC.seurat, slot = "data")),
  genes.use = VariableFeatures(MC.seurat),  # or a new set of HVG can be computed
  supercell_size = MC_list$supercell_size, # provide this parameter to run sample-weighted version of PCA,
  k = 30
)

#calculate umap
MC_list$UMAP <- supercell_UMAP(
  SC = MC_list,
  PCA_name = "PCA",
  n.comp = 30, n_neighbors = 15, min_dist=0.5
)

## Plot dimensional reduction
# Option 1
supercell_DimPlot(SC = MC_list,
                  groups = MC.seurat@meta.data$Cluster,
                  dim.name = "PCA",
                  title = "UMAP of metacells colored by cell type assignment"
) #+ theme()

# Option 2
#llamar supercell_dimplot, pero con tu nueva versión de plot_theme
supercell_Dimplot_correction <- function(...) {
  theme(aspect.ratio=4/3)
}

supercell_Dimplot_correction(SC = MC_list,
                             groups = MC.seurat@meta.data$Cluster,
                             dim.name = "PCA",
                             title = "UMAP of metacells colored by cell type assignment"
)

# Option 3
# add reduction PCA to MC.seurat
embedings_pca_w <- MC_list$PCA$x #embeddings
dim(embedings_pca_w)
#label the columns to ensure downstream consistency
colnames(embedings_pca_w) <- paste0("PCW_", 1:30)
# store dimensional reduction called 'embedings_pca_w'
MC.seurat[["pca_weighted"]] <- CreateDimReducObject(embeddings = embedings_pca_w, key = "PCW_", assay = "RNA")
MC.seurat@reductions$pca_weighted@stdev <- MC_list$PCA$sdev #stadar deviation
MC.seurat@reductions$pca_weighted@feature.loadings <- MC_list$PCA$rotation #loadings

# add reduction UMAP to MC.seurat
embedings_umap_pca_w <- MC_list$UMAP$layout #embeddings
dim(embedings_umap_pca_w)
#label the columns to ensure downstream consistency
colnames(embedings_umap_pca_w) <- paste0("UMAPW_", 1:2)
# store dimensional reduction called 'embedings_pca_w'
MC.seurat[["umap_pca_weighted"]] <- CreateDimReducObject(embeddings = embedings_umap_pca_w, key = "UMAPW_", assay = "RNA")

#plot reductions
Idents(MC.seurat) <- "Cluster"
levels(MC.seurat) <- sort(levels(MC.seurat))
DimPlot(MC.seurat, reduction = "pca_weighted") + ggtitle("pca_weighted_supercell")
ggsave(paste0(outdir,"/pca_weighted_Clusters_metacells_data.png"),width = 10,height = 8, dpi = 300)

DimPlot(MC.seurat, reduction = "umap_pca_weighted") + ggtitle("umap_pca_weighted_supercell")
ggsave(paste0(outdir,"/umap_pca_weighted_Clusters_metacells_data.png"),width = 10,height = 8, dpi = 300)

MC.seurat <- RunUMAP(MC.seurat, reduction = "pca_weighted", dims = 1:10, reduction.name = "umap_pca_weighted_seurat", reduction.key="UMAPW_SEURAT_")
DimPlot(MC.seurat, reduction = "umap_pca_weighted_seurat") + ggtitle("umap pca_weighted_seurat")
ggsave(paste0(outdir,"/umap_pca_weighted_seurat_Clusters_metacells_data.png"),width = 10,height = 8, dpi = 300)

################################################################################
## Differential expression analysis
################################################################################

# Compute upregulated genes in each cell line (versus other cells)
comparison_result_metacells_deseq2 <- realizar_comparacion_todos(MC.seurat, test_use = "DESeq2", only_pos = TRUE, logfc_threshold = log2(1.5))
write.csv(comparison_result_metacells_deseq2, file = "markers_metacells_deseq2.csv", row.names = TRUE) #save markers

comparison_result_metacells_mast <- realizar_comparacion_todos(MC.seurat.lognorm, test_use = "MAST", only_pos = TRUE, logfc_threshold = log2(1.5))
write.csv(comparison_result_metacells_mast, file = "markers_metacells_mast.csv", row.names = TRUE) #save markers

# DE analysis t-test weighted
head(MC.ge.lognorm)==head(MC.seurat.lognorm@assays$RNA$data) #TRUE!
MC.all.markers <- supercell_FindAllMarkers(
  ge = MC.ge.lognorm, 
  clusters = MC$Cluster, 
  supercell_size = MC$supercell_size,
  only.pos = TRUE, 
  min.pct = 0.25, 
  logfc.threshold = log2(1.5)
)

# Transform the output of `supercell_FindAllMarkers()` to be in the format of the `Seurat::FindAllMarkers()`
MC.all.markers.df <- data.frame()
for(cl in names(MC.all.markers)){
  cur <- MC.all.markers[[cl]]
  cur$cluster <- cl
  cur$gene <- rownames(cur)
  cur$avg_log2FC <- cur$logFC
  MC.all.markers.df <- rbind(MC.all.markers.df, cur)
}
write.csv(MC.all.markers.df, file = "markers_metacells_t_test_weighted.csv", row.names = TRUE) #save markers

# DE analysis wilcoxon-test weighted(???)
comparison_result_metacells_wilcox <- realizar_comparacion_todos(MC.seurat, test_use = "wilcox", only_pos = TRUE, logfc_threshold = log2(1.5))
write.csv(comparison_result_metacells_wilcox, file = "markers_metacells_wilcox.csv", row.names = TRUE) #save markers

#continuar en la linea 163 mirando la expresion diferencial

################################################################################
################################################################################
# Todos los genes
all_genes <- rownames(data_1_cl)

# Definir genes esperados para cada clúster
expected_degs_filtered <- list()
for (i in 1:6) {
  cluster_name <- paste0("Cluster_", i-1)
  expected_degs_filtered[[cluster_name]] <- rownames(data_1_cl)[rowData(data_1_cl)[[paste0("DEFacGroup", i)]] > 1.5]
} #revisar para varios DEFacGroup

# Crear una tabla con el número de DEGs esperados en cada clúster
deg_counts_filtered <- sapply(expected_degs_filtered, length)
deg_counts_df_filtered_1 <- data.frame(
  Cluster = names(deg_counts_filtered),
  Num_DEGs_Overexpressed = deg_counts_filtered
)

# Mostrar la tabla con el número de DEGs esperados
print(deg_counts_df_filtered_1) # DeFacGruop Foldchange > 1.5


# Crear tabla de resumen para DEFacGroup > 1.5
# Funcion modificada, se añadío sim.groups
summary_over_df_filtered_1 <- crear_resumen(expected_degs_filtered, data_1_cl)

# Mostrar las tablas de resumen
print(head(summary_over_df_filtered_1)) # DeFacGruop Foldchange > 1.5
print(dim(summary_over_df_filtered_1)) # 2304    3

##linea 386
basedir <- "C:/Users/Administrator/Documents/Jaguar/Simulation_DE/Metacells"
comparison_result_wilcox <- read.csv(paste0(basedir,"/markers_metacells_wilcox.csv"))
comparison_result_t_test <- read.csv(paste0(basedir,"/markers_metacells_t_test_weighted.csv"))
comparison_result_deseq2 <- read.csv(paste0(basedir,"/markers_metacells_deseq2.csv"))
comparison_result_mast <- read.csv(paste0(basedir,"/markers_metacells_mast.csv"))

colnames(comparison_result_t_test)[3] = "p_val_adj"

a=data_1_cl
a=MC.seurat
a <- AggregateExpression(a, return.seurat = TRUE, group.by = c('Cluster'))
head(a@assays$RNA$counts)

Idents(a)="Cluster"
comparison_result_deseq2_1 <- FindAllMarkers(object = a, test.use = "DESeq2")
comparison_result_deseq2_1

cts <- a@assays$RNA$counts
head(cts)
coldata <- data.frame(Cluster=as.factor(a$Cluster))
head(coldata)
unique(coldata)
class(coldata$Cluster)

comparison_result_deseq2_2 <- DESeqDataSetFromMatrix(countData = cts,
                                                     colData = coldata,
                                                     design = ~ Cluster)

comparison_result_deseq2_2 <- DESeq(comparison_result_deseq2_2)
res <- results(comparison_result_deseq2_2)
res
