## Metacells

## Set environment - Change for your own directory
basedir <- "C:/Users/Administrator/Documents/Jaguar/Simulation_DE/Metacells"
setwd("C:/Users/Administrator/Documents/Jaguar/Simulation_DE")
outdir <- "C:/Users/Administrator/Documents/Jaguar/Simulation_DE/Metacells/figures_metacells"
source("C:/Users/Administrator/Documents/Jaguar/Simulation_DE/Functions.R")

#gc() #libera memoria automáticamente cuando un objeto ya no se utiliza

# Libraries
library(SuperCell)
library(Seurat)
library(SeuratData)
library(tidyverse)
library(dplyr)
library(caret)
library(vioplot)
library(ggplot2)
library(reshape2)
library(MetacellAnalysisToolkit)
library(data.table)
library(DESeq2)
library(openxlsx)
library(pROC)
library(PRROC)
library(dplyr)

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
head(data_1_cl@assays$RNA$data)
head(data_1_cl@assays$RNA$counts)
head(data_1_cl@meta.data)
head(data_1_cl@assays$RNA@meta.features)

# Preprocessing
data_1_cl <- NormalizeData(data_1_cl, normalization.method = "LogNormalize", scale.factor = 10000)
data_1_cl <- FindVariableFeatures(data_1_cl, selection.method = "vst", nfeatures = 2000)
data_1_cl <- ScaleData(data_1_cl)
data_1_cl <- RunPCA(data_1_cl, features = VariableFeatures(object = data_1_cl))
data_1_cl <- RunUMAP(data_1_cl, dims = 1:10)

# Plot variable features
hvg <- VariableFeatures(data_1_cl)
length(hvg) #2000
plot1 <- VariableFeaturePlot(data_1_cl)
LabelPoints(plot = plot1, points = hvg[1:20], repel = TRUE)
ggsave(paste0(outdir,"/plot_variable_features_single_cell_data.png"),width = 10,height = 8, dpi = 300)

# Plot PCA (2D representation of scRNA-seq data) colored by cell line
DimPlot(data_1_cl, reduction = "pca", group.by = "Cluster",label=T, raster=F) + ggtitle("PCA_single_cell_by_clusters")
ggsave(paste0(outdir,"/pca_single_cell_data.png"),width = 10,height = 8, dpi = 300)

# Plot UMAP
DimPlot(data_1_cl, reduction = "umap", group.by = "Cluster",label=F, raster=F) + 
  ggtitle("UMAMP_single_cell_by_clusters") + 
  theme(legend.text = element_text(size = 23),
        #legend.title = element_text(size = 10),
        plot.title = element_text(hjust = 0.5, size=30, face='bold'),
        axis.title.y = element_text(size = 25, face='bold'),
        axis.title.x = element_text(size = 25, face='bold'),
        axis.text.y = element_text(size = 20),
        axis.text.x = element_text(size = 20)
  )
ggsave(paste0(outdir,"/umap_single_cell_data.png"),width = 10,height = 8, dpi = 300)

################################################################################
## Expected DEGs per group and level
## This only take into account the DEFacGroup. Instead we use other option that
## takes into account both the DEFacGroup and the FoldChange
# Here we see cluster names so we'll use it in expected_degs_levels and all_genes
#class(data_1_cl$Cluster) #factor
#levels(data_1_cl$Cluster) #"Group1" until 6
#
## Define expected genes for each cluster
#all_genes <- list()
#expected_degs_levels <- list()
#
#FoldChange = DEFacGroup1 / ((DEFacGroup2 * nCellsGroup2  + DEFacGroup3 * nCellsGroup3) / (nCellsGroup2 + nCellsGroup3))
#
#for (i in 1:6) {
#  group_name <- paste0("Group", i)
#  defac_group = data_1_cl@assays$RNA@meta.features[c("Gene",paste0("DEFacGroup", i))]
#  colnames(defac_group)[2] = "DEFacGroup"
#  defac_group$cluster = group_name
#  colnames(defac_group)[1] = "Gene"
#  all_genes[[group_name]] <- defac_group[c("Gene","cluster")]
#  expected_degs_levels[[group_name]] <- DE_classes(defac_group, type="real")
#}
#names(all_genes)
#head(all_genes[["Group1"]])
#
#names(expected_degs_levels)
#head(expected_degs_levels[["Group1"]])


################################################################################
## MODIFICAR AQUI
################################################################################
# Idea from comment of contributor of splatter
# https://github.com/Oshlack/splatter/issues/57#issuecomment-427703128

# List of all DEFac columns
df = data_1_cl@assays$RNA@meta.features
DEFac_columns = colnames(df)[grepl("^DEFac", colnames(df))]

# number of cells per group
ncells_group = table(data_1_cl$Cluster)
names(ncells_group) = DEFac_columns

# Multiply ncell per DEFACgroup
for(namecol in DEFac_columns){
  df[paste0(namecol,"_ncell")] = df[namecol]*ncells_group[namecol]
}
colnames(df)

# Initialize an empty list to store results
result_list <- list()

# Loop through each DEFac group as First group
for (i in seq_along(DEFac_columns)) {
  
  # Calculate FirstDEFac for the current group
  FirstDEFac_col <- DEFac_columns[i]
  columns_SecondDEFac <- paste0(DEFac_columns[-i],"_ncell")
  ncell_SecondDEFac <- ncells_group[-i]
  sum_ncells_SecondDEFac = sum(ncell_SecondDEFac)
  
  # Compute SecondDEFac as weighted average of Groups
  df_temp <- df %>%
    mutate(cluster = paste0("Group",i),
           FirstDEFac = .[[FirstDEFac_col]],
           SecondDEFac = rowSums(across(all_of(columns_SecondDEFac))) / sum_ncells_SecondDEFac,  # we use across to apply function to all columns indicated and use all_off as helper function to ensure dplyr selects specified columns
           log2FoldChange = log2(FirstDEFac / SecondDEFac)) %>%
    select(Gene, log2FoldChange, DEFAcGroup = FirstDEFac_col, cluster)
  
  # Append the result to the list
  cluster_name = paste0("Group",i)
  result_list[[cluster_name]] <- df_temp
}

# Combine all results into one data frame
expected_degs_levels_df <- do.call(rbind, result_list)
head(expected_degs_levels_df)

# Here we can verify DEFACGROUPS is different to log2FOLDCHANGE
plot(expected_degs_levels_df$DEFAcGroup,2^(expected_degs_levels_df$log2FoldChange))
lines(1:70,1:70, col="blue")

# Revise relation between positive DE genes and FoldChange.
# The values in DEfacGroup which are higher than one (positively differential expressed) has 
# higher values of foldchange than those DEfacGroup values equal to one or lower than one (negatively differential expressed).
# This has sense due to we are simulating cells principally with positive DE genes
expected_degs_levels_copy = expected_degs_levels_df %>%
  mutate(log2FoldChange_1.2 = ifelse(log2FoldChange >= log2(1.2), "red", "black"), #see diference with the threashold we are using in foldchange vs defacgroup higer than 1
         DEFAcGroup_DE = DEFAcGroup>1)

plot(expected_degs_levels_copy$DEFAcGroup_DE, expected_degs_levels_copy$log2FoldChange, col= expected_degs_levels_copy$log2FoldChange_1.2)
table(expected_degs_levels_copy$log2FoldChange_1.2, expected_degs_levels_copy$DEFAcGroup_DE)
# Now, if we compare with the table DEFacGroup and FoldChange, we can see that Foldchange mark
# as differentially expressed some genes that ARE NOT DE in the simulation! This could be done
# due to technical variation, because the biological simulation variation do not mark those
# genes as DE, it means are false positives. So which values should we take?: Values of
# DEfacGroup higher than 1 (TRUE DE) and if we want to know how big is the DE, we should
# use FoldChange

# So, let's filter the DE genes and lets add the metrics low, medium and high
# Define expected genes for each cluster
expected_degs_levels = DE_classes(expected_degs_levels_df, type="real")
table(expected_degs_levels$logfc_level)
View(expected_degs_levels)

all_genes <- expected_degs_levels_df %>%
  select(Gene, cluster)

# View the final result
head(expected_degs_levels)
dim(expected_degs_levels) #3082 4
head(all_genes)
dim(all_genes) #90000 2

#save expected genes
write.csv(expected_degs_levels, file = "DEGs_expected.csv", row.names = TRUE)
write.csv(all_genes, file = "all_genes.csv", row.names = TRUE)

## Crear una tabla con el número de DEGs esperados por nivel en cada clúster
# real data
head(expected_degs_levels)
deg_counts_filtered <- expected_degs_levels %>%
  group_by(cluster,logfc_level) %>%
  summarise(count = n())

View(deg_counts_filtered)

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

# Compute metacells using SuperCell package
MC <- SCimplify(
  X = GE, # single-cell log-normalized gene expression data
  genes.use = hvg, #highly variable genes
  cell.annotation	= Cluster, #vector with groups of the cells for metacells construction
  gamma = gamma,
  k.knn = k.knn
)

#network plot
png(paste0(outdir,"/network_metacell_data.png"))
supercell_plot(
  MC$graph.supercells, #networks
  group = MC$Cluster, 
  seed = 1, 
  alpha = -pi/2,
  main  = "Metacells"
  #title_size = 20
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
  print("Averaged matrix calculated")
}

head(rownames(MC.ge))
colnames(MC.ge) = paste0("MetaCell",1:ncol(MC.ge))
head(colnames(MC.ge))
head(MC.counts)
head(MC.ge)

#############################################################################
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
) + theme(plot.title = element_text(hjust = 0.5, size=30, face='bold'))
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
  main  = "Single cells by groups assignment"
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

png(paste0(outdir,"/histogram_purity_metacell_data.png"),
    width=500, height=500)
# Set graphical parameters using par() for the axis and text sizes
par(cex.main = 2,            # Title size
    font.main = 2,           # Bold title
    cex.lab = 1.5,           # Axis labels size
    font.lab = 2,            # Bold axis labels
    cex.axis = 1.4,          # Axis text size
    las = 1)                 # Rotate axis labels (optional, 1 = horizontal)

# Create the histogram
hist(MC$purity, 
     main = paste0("Purity of metacells (", method_purity, ")"), 
     xlab = "Purity", 
     ylab = "Frequency")
dev.off()

png(file=paste0(outdir,"/boxplot_purity.png"),
    width=500, height=500)
# Set graphical parameters using par() for the axis and text sizes
par(cex.main = 2,            # Title size
    font.main = 2,           # Bold title
    cex.lab = 1.5,           # Axis labels size
    font.lab = 2,            # Bold axis labels
    cex.axis = 1.4,          # Axis text size
    las = 1)                 # Rotate axis labels (optional, 1 = horizontal)

# Create the boxplot
boxplot(MC$purity, 
        main = paste0("Purity of metacells (", method_purity,")"))
dev.off() # a function call to save the file

png(file=paste0(outdir,"/boxplot_Clusters_purity.png"),
    width=500, height=500)
# Set graphical parameters using par() for the axis and text sizes
par(cex.main = 2,            # Title size
    font.main = 2,           # Bold title
    cex.lab = 1.5,           # Axis labels size
    font.lab = 2,            # Bold axis labels
    cex.axis = 1.4,          # Axis text size
    las = 1)                 # Rotate axis labels (optional, 1 = horizontal)
# Plot boxplot
boxplot(MC$purity ~ MC$Cluster, 
        main = paste0("Purity of metacells \nin terms of cell groups (", method_purity,")"))
dev.off() # a function call to save the file

#size metacells: Having a homogeneous metacell size distribution is ideal for 
#downstream analyses, since larger metacells will express more genes, which 
#could confound analyses. When heterogeneous size distributions are obtained 
#we recommend weighted downstream analyses 
png(file=paste0(outdir,"/histplot_size_distribution.png"),
    width=500, height=500)
# Set graphical parameters using par() for the axis and text sizes
par(cex.main = 2,            # Title size
    font.main = 2,           # Bold title
    cex.lab = 1.5,           # Axis labels size
    font.lab = 2,            # Bold axis labels
    cex.axis = 1.4,          # Axis text size
    las = 1)      
# plot hist
hist(MC$supercell_size, main = "Size distribution", xlab = "Size", breaks = 50)
dev.off() # a function call to save the file

png(file=paste0(outdir,"/violinplot_size_distribution.png"),
    width=500, height=500)
# Set graphical parameters using par() for the axis and text sizes
par(cex.main = 2,            # Title size
    font.main = 2,           # Bold title
    cex.lab = 1.5,           # Axis labels size
    font.lab = 2,            # Bold axis labels
    cex.axis = 1.4,          # Axis text size
    las = 1)      
# plot hist 
vioplot(MC$supercell_size, main = "Size distribution", xlab = "Size")
dev.off() # a function call to save the file

# Other QC metrics
# https://gfellerlab.github.io/MetacellAnalysisTutorial/downstream-analysis.html#weighted-analysis
library(reticulate)
library(spatstat.utils)

if(packageVersion("Seurat") >= 5) {options(Seurat.object.assay.version = "v4"); print("you are using seurat v5 with assay option v4")}
library(anndata)
library(MetacellAnalysisToolkit)
library(ggplot2)

# compactness
MC$compactness <- mc_compactness(cell.membership = data.frame(membership = MC$membership), 
                                 sc.obj = data_1_cl,
                                 sc.reduction = "pca", 
                                 dims = 1:30)

png(file=paste0(outdir,"/boxplot_compactness.png"),
    width=500, height=500)
# Set graphical parameters using par() for the axis and text sizes
par(cex.main = 2,            # Title size
    font.main = 2,           # Bold title
    cex.lab = 1.5,           # Axis labels size
    font.lab = 2,            # Bold axis labels
    cex.axis = 1.4,          # Axis text size
    las = 1)                 # Rotate axis labels (optional, 1 = horizontal)
# Plot boxplot
boxplot(MC$compactness, 
        main = paste0("Compactness of metacells"))
dev.off() # a function call to save the file

png(file=paste0(outdir,"/boxplot_Cluster_compactness.png"),
    width=500, height=500)
# Set graphical parameters using par() for the axis and text sizes
par(cex.main = 2,            # Title size
    font.main = 2,           # Bold title
    cex.lab = 1.5,           # Axis labels size
    font.lab = 2,            # Bold axis labels
    cex.axis = 1.4,          # Axis text size
    las = 1)                 # Rotate axis labels (optional, 1 = horizontal)
# Plot boxplot
boxplot(MC$compactness ~ MC$Cluster, 
        main = paste0("Compactness of metacells"))
dev.off() # a function call to save the file

#difussion map: We can also compute the compactness of each metacell using diffusion map components computed based on the PCA axes, as suggested in (Persad et al. 2023). 
diffusion_comp <- get_diffusion_comp(sc.obj = data_1_cl, sc.reduction = "pca", dims = 1:30)

MC$compactness <- mc_compactness(cell.membership = data.frame(membership = MC$membership), 
                                 sc.obj = data_1_cl,
                                 sc.reduction = diffusion_comp, 
                                 dims = 1:ncol(diffusion_comp))
boxplot(MC$compactness ~ MC$Cluster, 
        main = paste0("Compactness of metacells"))

#separation
#The separation of a metacell is the distance to the closest metacell (Persad et al. 2023). The higher the separation value the better.
MC$separation <- mc_separation(cell.membership = data.frame(membership = MC$membership), 
                                    sc.obj = data_1_cl,
                                    sc.reduction = "pca")

png(file=paste0(outdir,"/boxplot_separation.png"),
    width=500, height=500)
# Set graphical parameters using par() for the axis and text sizes
par(cex.main = 2,            # Title size
    font.main = 2,           # Bold title
    cex.lab = 1.5,           # Axis labels size
    font.lab = 2,            # Bold axis labels
    cex.axis = 1.4,          # Axis text size
    las = 1)                 # Rotate axis labels (optional, 1 = horizontal)
# Plot boxplot
boxplot(MC$separation, 
        main = paste0("Separation of metacells"))
dev.off() # a function call to save the file


png(file=paste0(outdir,"/boxplot_Cluster_separation.png"),
    width=500, height=500)
# Set graphical parameters using par() for the axis and text sizes
par(cex.main = 2,            # Title size
    font.main = 2,           # Bold title
    cex.lab = 1.5,           # Axis labels size
    font.lab = 2,            # Bold axis labels
    cex.axis = 1.4,          # Axis text size
    las = 1)                 # Rotate axis labels (optional, 1 = horizontal)
# Plot boxplot
boxplot(MC$separation ~ MC$Cluster, 
        main = paste0("Separation of metacells"))
dev.off() # a function call to save the file

#Note that compactness and separation metrics are correlated, better compactness results in worse separation and vice versa. 
ggplot(data.frame(compactness = log(MC$compactness), separation = log(MC$separation)), 
       aes(x=compactness, y=separation)) + 
  geom_point() +
  geom_smooth(method=lm) + 
  ggpubr::stat_cor(method = "pearson")

#The inner normalized variance (INV) of a metacell is the mean-normalized variance of gene expression within the metacell (Ben-Kiki et al. 2022).
#The lower the INV value the better. Note that it is the only metric that is latent-space independent.
MC$INV <- mc_INV(cell.membership = data.frame(membership = MC$membership),
                 sc.obj = data_1_cl,
                 group.label = "membership")

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
head(MC.seurat@assays$RNA$counts,1)
#Gene1  .  3  4 .  9 
head(MC.seurat@assays$RNA$data,1)
#Gene1  .  3  4 .  9
head(MC.seurat@assays$RNA$scale.data,1)
#Gene1 -0.5640929  0.7585165  1.1993864

head(MC.seurat.lognorm@assays$RNA$counts,1)
#Gene1 .         0.2323282 0.19951878
head(MC.seurat.lognorm@assays$RNA$data,1)
#Gene1 .         0.2323282 0.19951878 
head(MC.seurat.lognorm@assays$RNA$scale.data,1)
#Gene1 -0.64866307  0.6260241  0.44601286

#pca and pca_weighted are exactly the same in this point if calculated with 10 pc
#any(!c(MC.seurat@reductions$pca@cell.embeddings == MC.seurat@reductions$pca_weighted@cell.embeddings))

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
DimPlot(MC.seurat.lognorm, reduction = "pca_weighted") + 
  ggtitle("PCA_metacells_lognorm") + 
  theme(legend.text = element_text(size = 23),
        #legend.title = element_text(size = 10),
        plot.title = element_text(hjust = 0.5, size=30, face='bold'),
        axis.title.y = element_text(size = 25, face='bold'),
        axis.title.x = element_text(size = 25, face='bold'),
        axis.text.y = element_text(size = 20),
        axis.text.x = element_text(size = 20)
  )
ggsave(paste0(outdir,"/pca_Clusters_metacells_data_lognorm.png"),width = 10,height = 8, dpi = 300)

#now saving umap for MC.seurat.lognorm
MC.seurat.lognorm <- RunUMAP(MC.seurat.lognorm, reduction = "pca", dims = 1:10)
DimPlot(MC.seurat.lognorm, reduction = "umap") + 
  ggtitle("UMAP_metacells_lognorm") + 
  theme(legend.text = element_text(size = 23),
        #legend.title = element_text(size = 10),
        plot.title = element_text(hjust = 0.5, size=30, face='bold'),
        axis.title.y = element_text(size = 25, face='bold'),
        axis.title.x = element_text(size = 25, face='bold'),
        axis.text.y = element_text(size = 20),
        axis.text.x = element_text(size = 20)
  )
ggsave(paste0(outdir,"/umap_Clusters_metacells_data_lognorm.png"),width = 10,height = 8, dpi = 300)

# Seurat object containing normalized metacells gene expression (? NOT) data as well as 
# the first (N.comp) principal components of PCA performed internally using user 
# defined set of genes (by default the genes used for metacells constructions)
if(packageVersion("Seurat") >= 5) {
  MC.seurat[["RNA"]] <- as(object = MC.seurat[["RNA"]], Class = "Assay")
  MC.seurat.lognorm[["RNA"]] <- as(object = MC.seurat.lognorm[["RNA"]], Class = "Assay")
}

# Add gene metadata, so we can use DEFacGroup from simulation
MC.seurat@assays$RNA@meta.features = cbind(MC.seurat@assays$RNA@meta.features, data@rowRanges@elementMetadata)
#colnames(MC.seurat@assays$RNA@meta.features)[1] = "gene"
MC.seurat.lognorm@assays$RNA@meta.features = cbind(MC.seurat.lognorm@assays$RNA@meta.features, data@rowRanges@elementMetadata)
#colnames(MC.seurat.lognorm@assays$RNA@meta.features)[1] = "gene"

# save as Rds
saveRDS(MC.seurat, file = file.path(paste0(basedir, "/MC_gamma_", gamma, "_seurat_object.Rds")))
saveRDS(MC.seurat.lognorm, file = file.path(paste0(basedir, "/MC_gamma_", gamma, "_seurat_lognorm.Rds")))

################################################################################
## Dimensional reduction weighted
################################################################################
# read data
MC.seurat <- readRDS(paste0(basedir,"/MC_gamma_20_seurat_object.Rds"))
MC.seurat.lognorm <- readRDS(paste0(basedir,"/MC_gamma_20_seurat_lognorm.Rds"))

# MC.seurat is DIFFERENT to MC.seurat.lognormalized (only check, not relevant)
any(!c(MC.seurat@assays$RNA$data == MC.seurat.lognorm@assays$RNA$data)[[1]])

dim(MC.seurat) #15000  200
head(MC.seurat@assays$RNA@meta.features)

#normalization and variable features
MC.seurat <- NormalizeData(MC.seurat, normalization.method = "LogNormalize")
MC.seurat <- FindVariableFeatures(MC.seurat)
MC.seurat <- ScaleData(MC.seurat)

# MC.seurat normalized is EQUAL to MC.seurat.lognormalized (only check, not relevant)
any(!c(MC.seurat@assays$RNA$data == MC.seurat.lognorm@assays$RNA$data)[[1]])

#calculate PCA weighted
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
DimPlot(MC.seurat, reduction = "pca_weighted") + 
  ggtitle("PCA_weighted_supercell") +
  theme(legend.text = element_text(size = 23),
        #legend.title = element_text(size = 10),
        plot.title = element_text(hjust = 0.5, size=30, face='bold'),
        axis.title.y = element_text(size = 25, face='bold'),
        axis.title.x = element_text(size = 25, face='bold'),
        axis.text.y = element_text(size = 20),
        axis.text.x = element_text(size = 20)
  )
ggsave(paste0(outdir,"/pca_weighted_Clusters_metacells_data.png"),width = 10,height = 8, dpi = 300)

DimPlot(MC.seurat, reduction = "umap_pca_weighted") + 
  ggtitle("UMAP_PCA_weighted_supercell") + 
  theme(legend.text = element_text(size = 23),
        #legend.title = element_text(size = 10),
        plot.title = element_text(hjust = 0.5, size=30, face='bold'),
        axis.title.y = element_text(size = 25, face='bold'),
        axis.title.x = element_text(size = 25, face='bold'),
        axis.text.y = element_text(size = 20),
        axis.text.x = element_text(size = 20)
  )
ggsave(paste0(outdir,"/umap_pca_weighted_Clusters_metacells_data.png"),width = 10,height = 8, dpi = 300)

MC.seurat <- RunUMAP(MC.seurat, reduction = "pca_weighted", dims = 1:10, reduction.name = "umap_pca_weighted_seurat", reduction.key="UMAPW_SEURAT_")
DimPlot(MC.seurat, reduction = "umap_pca_weighted_seurat") + 
  ggtitle("UMAP_PCA_weighted_seurat") + 
  theme(legend.text = element_text(size = 23),
        #legend.title = element_text(size = 10),
        plot.title = element_text(hjust = 0.5, size=30, face='bold'),
        axis.title.y = element_text(size = 25, face='bold'),
        axis.title.x = element_text(size = 25, face='bold'),
        axis.text.y = element_text(size = 20),
        axis.text.x = element_text(size = 20)
  )
ggsave(paste0(outdir,"/umap_pca_weighted_seurat_Clusters_metacells_data.png"),width = 10,height = 8, dpi = 300)

#plot umap with sized
data <- cbind(MC.seurat@reductions$umap_pca_weighted_seurat@cell.embeddings, MC.seurat$size)
colnames(data)[1:2] <- c("umap_1", "umap_2")
colnames(data)[3] <- c("size")#, "cluster")
ggplot(data,aes(x = umap_1, y = umap_2, color = MC.seurat$Cluster)) + 
  geom_point(aes(size=size)) + 
  ggtitle("UMAP_PCA_weighted_seurat") +
  labs(
    color = "Cluster") +
  theme_classic() +
  theme(legend.text = element_text(size = 18),
        legend.title = element_text(size = 20),
        plot.title = element_text(hjust = 0.5, size=30, face='bold'),
        axis.title.y = element_text(size = 25, face='bold'),
        axis.title.x = element_text(size = 25, face='bold'),
        axis.text.y = element_text(size = 20),
        axis.text.x = element_text(size = 20)
  ) +
  guides(color = guide_legend(override.aes = list(size = 4)))  # Ajustar tamaño de los puntos de la leyenda 'color'

ggsave(paste0(outdir,"/umap_pca_weighted_seurat_Clusters_size_metacells_data.png"),width = 10,height = 8, dpi = 300)

# save as Rds
saveRDS(MC.seurat, file = file.path(paste0(basedir, "/MC_gamma_", gamma, "_seurat_object_normalized.Rds")))

################################################################################
## Differential expression analysis
################################################################################

# read data
MC.seurat <- readRDS(paste0(basedir,"/MC_gamma_20_seurat_object_normalized.Rds"))
MC.seurat.lognorm <- readRDS(paste0(basedir,"/MC_gamma_20_seurat_lognorm.Rds"))
MC.seurat$Cluster <- as.factor(MC.seurat$Cluster)

## Calculate differential expression
# MC.seurat normalized is equal to MC.seurat.lognormalized (only check, not relevant)
any(!c(MC.seurat@assays$RNA$data == MC.seurat.lognorm@assays$RNA$data)[[1]])

##
# DE analysis t-test weighted
MC.all.markers <- supercell_FindAllMarkers(
  ge = MC.seurat@assays$RNA$data, 
  clusters = MC.seurat$Cluster, 
  supercell_size = MC.seurat$supercell_size,
  only.pos = TRUE, 
  min.pct = 0.01, 
  logfc.threshold = log2(1.2)
)

# Transform the output of `supercell_FindAllMarkers()` to be in the format of the `Seurat::FindAllMarkers()`
MC.all.markers.df <- data.frame()
for(cl in names(MC.all.markers)){
  cur <- MC.all.markers[[cl]]
  cur$cluster <- cl
  cur$Gene <- rownames(cur)
  #cur$avg_log2FC <- cur$logFC
  MC.all.markers.df <- rbind(MC.all.markers.df, cur)
}

# Change variable names to be the same as in seurat
colnames(MC.all.markers.df)[c(1,2,5)] = c("p_val","p_val_adj","avg_log2FC")
head(MC.all.markers.df)
dim(MC.all.markers.df) #1143  9
unique(MC.all.markers.df$cluster) #"Group1"  until 6
class(MC.all.markers.df$cluster) #character
# save results
write.csv(MC.all.markers.df, file = paste0(basedir,"/DE_metacells_t_test_weighted.csv"), row.names = TRUE) #save markers

##
# DE analysis t-test un-weighted
Idents(MC.seurat)
MC.seurat$Cluster
data_1_cl$Cluster
DE_metacells_t_unweighted = DE_analysis_seurat(MC.seurat, test_type = "t", assay="RNA", slot="data")
dim(DE_metacells_t_unweighted) #5517  7
levels(DE_metacells_t_unweighted$cluster) #"Group1" until 6
class(DE_metacells_t_unweighted$cluster) #factor
# save resuts
write.csv(DE_metacells_t_unweighted, file = paste0(basedir,"/DE_metacells_t_test_unweighted.csv"), row.names = TRUE) #save markers

## deseq2
# Compute upregulated genes in each cell line (versus other cells)
DE_metacells_deseq2 <- DE_analysis_seurat(MC.seurat, assay = "RNA", slot="counts", test_type="DESeq2")
# save results
write.csv(DE_metacells_deseq2, file = paste0(basedir,"/DE_metacells_deseq2.csv"), row.names = TRUE) #save markers

## deseq2 original package
coldata <- data.frame(Cluster=factor(MC.seurat$Cluster, levels = c("Group1", "Group2", "Group3", "Group4", "Group5", "Group6")))
rownames(coldata) <- colnames(MC.seurat)
class(coldata$Cluster) #factor

deseq2_matrix <- DESeqDataSetFromMatrix(
  countData = MC.seurat@assays$RNA$counts, 
  colData = coldata, 
  design = ~ 0+Cluster
)

DE_metacells_deseq2_original <- list()
for (level in levels(deseq2_matrix$Cluster)){
  print(level)
  # Relevel
  deseq2_matrix$Cluster <- relevel(deseq2_matrix$Cluster, ref = level)
  # Run deseq
  dds <- DESeq(deseq2_matrix)
  # Compare cluster selected against all others
  res <- results(dds, contrast = c(1, -1/5, -1/5, -1/5, -1/5, -1/5))
  # Rename and save in list
  results <- as.data.frame(res)
  results$cluster = level
  results$Gene = rownames(results)
  DE_metacells_deseq2_original[[level]] <- results
}
resultsNames(dds)
DE_metacells_deseq2_groups <- bind_rows(DE_metacells_deseq2_original)
setnames(DE_metacells_deseq2_groups, c("padj"), c("p_val_adj"))

# save results
write.csv(DE_metacells_deseq2_groups, file = paste0(basedir,"/DE_metacells_deseq2_original.csv"), row.names = TRUE) #save markers

## mast. It performs better with normalized and scaled data
DE_metacells_mast1 <- DE_analysis_seurat(MC.seurat, assay = "RNA", slot="counts", test_type = "MAST")
DE_metacells_mast2 <- DE_analysis_seurat(MC.seurat, assay = "RNA", slot="data", test_type = "MAST")
head(DE_metacells_mast1)
dim(DE_metacells_mast1)
head(DE_metacells_mast2)
dim(DE_metacells_mast2)
head(MC.seurat@assays$RNA$data)
head(MC.seurat@assays$RNA$scale.data)
# save results
write.csv(DE_metacells_mast2, file = paste0(basedir,"/DE_metacells_mast.csv"), row.names = TRUE) #save markers

## wilcoxon 
DE_metacells_wilcox <- DE_analysis_seurat(MC.seurat, assay = "RNA", slot="data", test_type = "wilcox")
# save results
write.csv(DE_metacells_wilcox, file = paste0(basedir,"/DE_metacells_wilcox.csv"), row.names = TRUE) #save markers

################################################################################
## Filter of DE results to define the markers
################################################################################

DE_metacells_t_unweighted = read.csv(paste0(basedir,"/DE_metacells_t_test_unweighted.csv"))
DE_metacells_t_weighted = read.csv(paste0(basedir,"/DE_metacells_t_test_weighted.csv"))
DE_wilcox = read.csv(paste0(basedir,"/DE_metacells_wilcox.csv"))
DE_mast = read.csv(paste0(basedir,"/DE_metacells_mast.csv"))
DE_deseq2 = read.csv(paste0(basedir,"/DE_metacells_deseq2.csv"))
DE_deseq2_original = read.csv(paste0(basedir,"/DE_metacells_deseq2_original.csv"))
# Expected DEGS per group and level
expected_degs_levels_df <- read.csv("DEGs_expected.csv")
all_genes_df <- read.csv("all_genes.csv")
all_genes_df <- all_genes_df %>% select(c(Gene, cluster))

# Select DEGs and divide them by classes
markers_metacells_t_unweighted = DE_classes(DE_metacells_t_unweighted, min_pct_difference = 0.1, max_p_val_adj = 0.05, type="predicted")
markers_metacells_t_weighted = DE_classes(DE_metacells_t_weighted, min_pct_difference = 0.1, max_p_val_adj = 0.05, type="predicted")
markers_metacells_wilcox = DE_classes(DE_wilcox, min_pct_difference = 0.1, max_p_val_adj = 0.05, type="predicted")
markers_metacells_mast = DE_classes(DE_mast, min_pct_difference = 0.1, max_p_val_adj = 0.05, type="predicted")
markers_metacells_deseq2_original = DE_classes(DE_deseq2_original, max_p_val_adj = 0.05, type="predicted_deseq")

head(DE_deseq2)
# as there are not differentially expressed genes, we do not need to filter anything

################################################################################
# Data DE results/metrics from simulation

# list for save metrics
results_list <- list()

# Calculate metrics
results_list[["t_weighted"]] <- DE_metrics(expected_degs_levels_df,markers_metacells_t_weighted, all_genes_df, test="t_weighted")
results_list[["t_unweighted"]] <- DE_metrics(expected_degs_levels_df,markers_metacells_t_unweighted, all_genes_df, test="t_unweighted")
results_list[["mast"]] <- DE_metrics(expected_degs_levels_df,markers_metacells_mast, all_genes_df, test="mast")
results_list[["wilcox"]] <- DE_metrics(expected_degs_levels_df,markers_metacells_wilcox, all_genes_df, test="wilcox")
results_list[["deseq2"]] <- DE_metrics(expected_degs_levels_df,markers_metacells_deseq2_original, all_genes_df, test="deseq2")

## Heatmap summary per test
Metrics_heatmap(results_list, outdir, text="metacells")

## Save the DEGs (given the filters applied), the DEGs per cluster and the metrics to an Excel file
wb <- createWorkbook()

for (test in names(results_list)) {
  DEGs <- results_list[[test]]$DEGs
  metrics <- results_list[[test]]$Metrics
  
  # DEGs
  addWorksheet(wb, paste(test, "DEGs", sep="_"))
  writeData(wb, sheet = paste(test, "DEGs", sep="_"), DEGs)
  
  # metrics
  addWorksheet(wb, paste(test, "Metrics", sep="_"))
  writeData(wb, sheet = paste(test, "Metrics", sep="_"), metrics)

  # predicted data by group in all test
  DEGs_per_group <- DEGs %>%
    group_by(cluster) %>% #logfc_level
    summarise(count = n())
  DEGs_per_group <- rbind(DEGs_per_group, data.frame(cluster="all", count = sum(DEGs_per_group$count)))
  
  addWorksheet(wb, paste(test, "DEGs_group", sep="_"))
  writeData(wb, sheet = paste(test, "DEGs_group", sep="_"), DEGs_per_group)
}

saveWorkbook(wb, file = paste0(basedir,"/DEGs_and_Metrics_metacells_all_tests.xlsx"), overwrite = TRUE)

################################################################################
## Generate graphic of metrics from all tests in one heatmap

## Read data
wb = loadWorkbook(paste0(basedir,"/DEGs_and_Metrics_metacells_all_tests.xlsx"))
# explore date
t_weighted_DEGs = read.xlsx(wb, sheet = "t_weighted_DEGs")
head(t_weighted_DEGs)
t_weighted_Metrics = read.xlsx(wb, sheet = "t_weighted_Metrics")
print(t_weighted_Metrics)

## Heatmap all tests in one
# only select pages with "metrics" and only select metrics for row "all"
data_all = data.frame()

for (name in names(wb)){
  if (grepl("*Metrics$", name)){
    data <- read.xlsx(wb, sheet = name)
    data_all = rbind(data_all,data %>% filter(Level == "all"))
  }
}

# fuction for plot heatmap
Metrics_heatmap_all = function(results_list, outdir, title){
  # Transformar la matriz en long format
  df = reshape2::melt(results_list)
  colnames(df) = c("Test", "Level", "Metric", "Value")
  df$Value[df$Metric=="FDR"] = 1-df$Value[df$Metric=="FDR"]
  df$Value_discrete = ifelse(df$Value<0.3,1,ifelse(df$Value>=0.7,3,df$Value))
  df$Value_discrete = ifelse(0.3<=df$Value & df$Value<0.7,2,df$Value_discrete)
  df$Value_discrete <- factor(df$Value_discrete, levels = c("1","2","3"))
  df$Test_level = paste(df$Test, df$Level, sep="_")
  
  # Gráfico
  heatmap_test_metrics = ggplot(df, aes(x = Metric, y = Test_level, fill = Value_discrete)) +
    geom_tile() +
    ggtitle(title) +
    theme(legend.text = element_text(size = 23),
          legend.title = element_blank(),
          plot.title = element_text(hjust = 0.5, size=30, face='bold'),
          axis.title.y = element_text(size = 25, face='bold'),
          axis.title.x = element_text(size = 25, face='bold'),
          axis.text.y = element_text(size = 20),
          axis.text.x = element_text(size = 20, angle = 90, vjust = 0.5, hjust=1),
          text = element_text(size=13)) +
    coord_fixed() +
    scale_fill_manual(breaks = levels(df$Value_discrete),
                      values = c("#FF99BF", "#FFFF99", "#A5EDFF"),
                      labels=c('Poor','Intermediate','Good'))
  print(heatmap_test_metrics)
  ggsave(filename = paste0(outdir,"/", sub(" ","_",title), ".jpg"), plot = heatmap_test_metrics,
         width = 10,height = 8, dpi = 300)
}

# plot heatmap
Metrics_heatmap_all(data_all, outdir, "Heatmap metacells tests mestrics")
