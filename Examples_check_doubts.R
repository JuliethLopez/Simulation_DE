################################################################################
# Tests and examples (not necessary - only backup)
################################################################################
# Little example in metacells about include clusters with genes
test_low_expected = expected_degs_filtered_df %>% subset(logfc_level == "bajo") %>% select(c(Gene, cluster))
test_low_predicted = markers_metacells_t_weighted %>% filter(logfc_level == "bajo")  %>% select(c(Gene, cluster))

intersec_results = intersect(test_low_predicted$Gene, test_low_expected$Gene) #two genes : "Gene4689"  "Gene14355"
test_low_expected[test_low_expected$Gene %in% intersec_results,]
#                        Gene cluster
#Gene4689...4689    Gene4689  Group1
#Gene14355...44355 Gene14355  Group3
test_low_predicted[test_low_predicted$Gene %in% intersec_results,]
#         Gene cluster
#8   Gene4689  Group2
#14 Gene14355  Group5
dplyr::intersect(test_low_predicted, test_low_expected) # no coincidence because were found in different groups

################################################################################
# DE analysis wilcoxon-test weighted(???)
#install.packages("wgsea", dependencies = TRUE) #it was removed from CRAN
#install.packages("cNORM", dependencies = TRUE) #not install correctly
#devtools::install_github("WLenhard/cNORM") #not install correctly
#wilcoxon(MC.seurat, MC.seurat$Cluster, weights = MC.seurat$size, binsize = 0.05)

#little example for wilcoxon with weights
weights=c(5,5,20,10,10)
ge=matrix(runif(25), 5,5)
colnames(ge)=c("C","C","C","O","O")
rownames(ge)=paste0(rep("gene",5),1:5)

# One-sample wilcoxon test
wilcox.test(ge[1,1:3],ge[1,4:5]) #no se rechaza la hipotesis nula, no hay diferencias en la expresion del gen1

# One-sample wilcoxon test
library(dplyr)
library(tidyverse)
library(rstatix)
library(ggpubr)
ge_melt=reshape2::melt(ge)
stat.test <- ge_melt[ge_melt$Var1=="gene1",] %>% 
  wilcox_test(value ~ Var2) %>%
  add_significance()
stat.test

# Apply weigts by replicating the data based on the weights
# CORRECT because gives more importances to certain values, traduced in more ranks for sum
C_weighted <- rep(ge[1,1:3], times = weights[1:3])  # First 3 weights for C group
O_weighted <- rep(ge[1,4:5], times = weights[4:5])  # Last 2 weights for O group
wilcox.test(C_weighted,O_weighted, exact = FALSE, correct = TRUE) #se rechaza la hipotesis nula, si hay diferencias en la expresión del gen1

# Run the Wilcoxon rank-sum test on the weighted data
# INCORRECT in this way because ranks keep the same
wilcox.test(ge[1,1:3],ge[1,4:5]) #no se rechaza la hipotesis nula, no hay diferencias en la expresion del gen1
wilcox.test(ge[1,1:3]*weights[1:3],ge[1,4:5]*weights[4:5]) 

################################################################################
#DESEQ duda, revision
DE_wilcox <- read.csv(paste0(basedir,"/markers_metacells_wilcox.csv"))
DE_t_test <- read.csv(paste0(basedir,"/markers_metacells_t_test_weighted.csv"))
DE_deseq2 <- read.csv(paste0(basedir,"/markers_metacells_deseq2.csv"))
DE_mast <- read.csv(paste0(basedir,"/markers_metacells_mast.csv"))

colnames(DE_t_test)[3] = "p_val_adj"

a=data_1_cl
a=MC.seurat
a <- AggregateExpression(a, return.seurat = TRUE, group.by = c('Cluster'))
head(a@assays$RNA$counts)

Idents(a)="Cluster"
DE_deseq2_1 <- FindAllMarkers(object = a, test.use = "DESeq2")
DE_deseq2_1

cts <- a@assays$RNA$counts
head(cts)
coldata <- data.frame(Cluster=as.factor(a$Cluster))
head(coldata)
unique(coldata)
class(coldata$Cluster)

DE_deseq2_2 <- DESeqDataSetFromMatrix(countData = cts,
                                      colData = coldata,
                                      design = ~ Cluster)

DE_deseq2_2 <- DESeq(DE_deseq2_2)
res <- results(DE_deseq2_2)
res
