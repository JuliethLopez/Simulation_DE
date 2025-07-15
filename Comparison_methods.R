## Set environment - Change for your own directory
setwd("C:/Users/Administrator/Documents/Jaguar/Simulation_DE")
outdir <- "C:/Users/Administrator/Documents/Jaguar/Simulation_DE"
source("C:/Users/Administrator/Documents/Jaguar/Simulation_DE/Functions.R")

## libraries
library(tidyverse)
library(dplyr)
library(caret)
library(ggplot2)
library(reshape2)
library(data.table)
library(openxlsx)
library(pROC)
library(PRROC)
library(dplyr)

## Read data
wb = loadWorkbook("Metacells/DEGs_and_Metrics_metacells_all_tests.xlsx")
wb_mast = loadWorkbook("Mast/DEGs_and_Metrics_MAST.xlsx")
wb_t_wilcox = loadWorkbook("Wilcoxon_t_Test/Results/DEGs_and_Metrics.xlsx")
wb_deseq2_groups = loadWorkbook("DESEQ2/DEGs_and_Metrics_Deseq2.xlsx")

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
    data$Test = paste0("MC_",data$Test)
    data_all = rbind(data_all,data %>% filter(Level == "all"))
  }
}

for (name in names(wb_mast)){
  if (grepl("*Metrics$", name)){
    data <- read.xlsx(wb_mast, sheet = name)
    data_all = rbind(data_all,data %>% filter(Level == "all"))
  }
}

for (name in names(wb_t_wilcox)){
  if (grepl("*Metrics$", name)){
    data <- read.xlsx(wb_t_wilcox, sheet = name)
    data_all = rbind(data_all,data %>% filter(Level == "all"))
  }
}

for (name in names(wb_deseq2_groups)){
  if (grepl("*Metrics$", name)){
    data <- read.xlsx(wb_deseq2_groups, sheet = name)
    data_all = rbind(data_all,data %>% filter(Level == "all"))
  }
}

# some adjustment of form
data_all <- data_all %>% select(-AUroc)
data_all$Test[data_all$Test=="t-test"] = "SC_t_test"
data_all$Test[data_all$Test=="wilcox"] = "SC_wilcoxon_test"
data_all$Test[data_all$Test=="deseq2"] = "SC_deseq2"
data_all$Test[data_all$Test=="mast"] = "SC_mast"
data_all$Test[data_all$Test=="MC_wilcox"] = "MC_wilcoxon"

# fuction for plot heatmap
Metrics_heatmap_all = function(results_list, outdir, title){
  # Transformar la matriz en long format
  df = reshape2::melt(results_list)
  #df$Test = as.factor(df$Test)
  #levels(df$Test) = c(data_all$Test[grep("*_MC$",data_all$Test)], data_all$Test[!grepl("*_MC$",data_all$Test)])
  colnames(df) = c("Test", "Level", "Metric", "Value")
  df$Value[df$Metric=="FDR"] = 1-df$Value[df$Metric=="FDR"]
  df$Value_discrete = ifelse(df$Value<0.3,1,ifelse(df$Value>=0.7,3,df$Value))
  df$Value_discrete = ifelse(0.3<=df$Value & df$Value<0.7,2,df$Value_discrete)
  df$Value_discrete <- factor(df$Value_discrete, levels = c("1","2","3"))
  
  # Gráfico
  heatmap_test_metrics = ggplot(df, aes(x = Metric, y = Test, fill = Value_discrete)) +
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
Metrics_heatmap_all(data_all, outdir, "Heatmap tests metrics")

