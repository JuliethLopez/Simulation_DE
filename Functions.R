# Function to perform differential expression analysis
# min.pct= only test genes that are detected in a minimum fraction of min.pct cells in either of the two populations
DE_analysis_seurat <- function(seurat_obj, test_type = "wilcox", only.pos=TRUE, min_pct = 0.01, logfc_threshold = log2(1.2), assay=NULL, slot=NULL) { 
  # DE performance
  DE_seurat_corrected <- FindAllMarkers(seurat_obj,
                                        assay = assay,
                                        slot = slot,
                                        min_pct = min_pct,
                                        test.use = test_type, #t, negbinom, poisson, LR, MAST, DESeq2, wilcox, wilcox_limma, bimod, roc
                                        only.pos = only.pos,
                                        logfc.threshold = logfc_threshold)
  print("Differential expression calculated")
  
  #correction of p value
  DE_seurat_corrected$p_val_adj <- p.adjust(DE_seurat_corrected$p_val, method = "fdr")
  print("Correction of p value calculated calculated")
  return(DE_seurat_corrected)
}

# function and filters to identify genes per group with expression levels of interest
DE_classes <- function(data_DE, max_p_val_adj = 0.05, min_pct_difference = 0.05, type = "predicted"){
  # division of DE levels over real data
  if (type == "real"){
  summary_over_df <- data_DE %>%
    filter(DEFAcGroup > 1) %>%
    mutate(logfc_level = case_when( # indicate how much it's expressed a gene in one cluster against others
             log2FoldChange >= log2(1.2) & log2FoldChange < log2(2) ~ "low",
             log2FoldChange >= log2(2) & log2FoldChange < log2(3) ~ "medium",
             log2FoldChange >= log2(3) ~ "high",
             TRUE ~ "tooLow"),
           cluster = as.factor(cluster))}
  
  # division of DE levels over predicted data
  if (type == "predicted"){
  summary_over_df <- data_DE %>%
    mutate(logfc_level = case_when(
             avg_log2FC >= log2(1.2) & avg_log2FC < log2(2) ~ "low",
             avg_log2FC >= log2(2) & avg_log2FC < log2(3) ~ "medium",
             avg_log2FC >= log2(3) ~ "high",
             TRUE ~ "tooLow"),
           cluster = as.factor(cluster)) %>%
    # filter over predicted data for find markers
    mutate(pct_difference = pct.1-pct.2) %>%
    filter(p_val_adj < max_p_val_adj) %>%
    filter(pct_difference > min_pct_difference)}
  
  # division of DE levels over predicted deseq data
  if (type == "predicted_deseq"){
    summary_over_df <- data_DE %>%
      mutate(logfc_level = case_when(
             log2FoldChange >= log2(1.2) & log2FoldChange < log2(2) ~ "low",
             log2FoldChange >= log2(2) & log2FoldChange < log2(3) ~ "medium",
             log2FoldChange >= log2(3) ~ "high",
             TRUE ~ "tooLow"),
           cluster = as.factor(cluster)) %>%
      # filter over predicted data for find markers in deseq results
      filter(p_val_adj < max_p_val_adj)}
  
  # locate gene column for rename it
  index_column = grep("gene",tolower(colnames(summary_over_df)))
  colnames(summary_over_df)[index_column] = "Gene"
  
  return(summary_over_df)
}

calcular_auroc <- function(markers_expected, markers_found) {
  labels_total = c()
  scores_total = c()
  for (level in levels(markers_expected$cluster)){
    founded = subset(markers_found,cluster == level)
    expected = subset(markers_expected,cluster == level)
    # Crear un vector binario para las filas (genes esperados por grupo)
    labels <- as.numeric(expected$Gene %in% founded$Gene)
    
    # Completar valores NA e infinitos como 1
    scores <- expected %>% left_join(founded, by="Gene") %>%
      select(p_val_adj)
    scores <- -log10(scores[,"p_val_adj"])
  
    no_valid_indices <- is.na(scores) | is.infinite(scores)
    
    labels = labels[!no_valid_indices]
    scores = scores[!no_valid_indices] 
    
    labels_total <- c(labels_total, labels)
    scores_total <- c(scores_total, scores)
  }
  
  # Calcular AUROC
  if (length(unique(labels_total))>1){
    roc_obj <- roc(labels_total, scores_total)
    auroc <- auc(roc_obj)
  } else {
    auroc=0
  }
  
  return(auroc)
}

# Funtion to calculate metrics
calculate_metrics <- function(real, predicted, all_genes) {
  predicted_gene_cluster = predicted %>% select(c(Gene, cluster))
  TP <- nrow(dplyr::intersect(real, predicted_gene_cluster)) # True Positives
  FP <- nrow(dplyr::setdiff(predicted_gene_cluster, real))  # False Positives
  FN <- nrow(dplyr::setdiff(real, predicted_gene_cluster)) # False Negatives
  TN <- nrow(dplyr::setdiff(dplyr::setdiff(all_genes, real), predicted_gene_cluster)) # True Negatives
  
  precision <- ifelse((TP + FP) > 0, TP / (TP + FP), 0)
  recall <- ifelse((TP + FN) > 0, TP / (TP + FN), 0)
  f1_score <- ifelse((precision + recall) > 0, 2 * (precision * recall) / (precision + recall), 0)
  accuracy <- (TP + TN) / (TP + TN + FP + FN)
  fdr <- ifelse((FP + TP) > 0, FP / (FP + TP), 0)
  auroc <- calcular_auroc(real, predicted)
  
  return(c(precision, recall, f1_score, accuracy, fdr, auroc))
}

# function to calculate metrics per level
DE_metrics = function(expected_df, founded_filtered_df, all_genes, test="wilcox"){
  # adjust tyoe data
  expected_df$cluster = as.factor(expected_df$cluster)
  
  # Filter DEGs based on actual logfc levels
  expected_low <- expected_df %>% subset(logfc_level == "low") %>% select(c(Gene, cluster))
  expected_medium <- expected_df %>% subset(logfc_level == "medium") %>% select(c(Gene, cluster))
  expected_high <- expected_df %>% subset(logfc_level == "high") %>% select(c(Gene, cluster))
  expected_all <- expected_df %>% subset(logfc_level != "tooLow") %>% select(c(Gene, cluster))
  
  # Count the number of genes in each category
  num_all_r <- nrow(expected_all)
  print(paste("There are expected", num_all_r, "DEG in total"))
  num_low_r <- nrow(expected_low)
  print(paste("There are expected", num_low_r, "DEG in low level"))
  num_medium_r <- nrow(expected_medium)
  print(paste("There are expected", num_medium_r, "DEG in mid level"))
  num_high_r <- nrow(expected_high)
  print(paste("There are expected", num_high_r, "DEG in high level"))
  
  # Filtrar DEGs según los niveles de logfc predicted
  predicted_low <- founded_filtered_df %>% filter(logfc_level == "low")  %>% select(c(Gene, cluster, p_val_adj))
  predicted_medium <- founded_filtered_df %>% filter(logfc_level == "medium")  %>% select(c(Gene, cluster, p_val_adj))
  predicted_high <- founded_filtered_df %>% filter(logfc_level =="high")  %>% select(c(Gene, cluster, p_val_adj))
  predicted_all <- founded_filtered_df %>% select(c(Gene, cluster, p_val_adj))
  
  # Count the number of genes in each category
  num_all <- nrow(predicted_all)
  print(paste("There were found", num_all, "DEG in total"))
  num_low <- nrow(predicted_low)
  print(paste("There were found", num_low, "DEG in low level"))
  num_medium <- nrow(predicted_medium)
  print(paste("There were found", num_medium, "DEG in mid level"))
  num_high <- nrow(predicted_high)
  print(paste("There were found", num_high, "DEG in high level"))
  
  # Calculate metrics
  print("Calculating which DEGs found are in expected")
  metrics_low <- calculate_metrics(expected_low, predicted_low, all_genes)
  metrics_medium <- calculate_metrics(expected_medium, predicted_medium, all_genes)
  metrics_high <- calculate_metrics(expected_high, predicted_high, all_genes)
  metrics_all <- calculate_metrics(expected_all, predicted_all, all_genes)
  print("calculating metrics DONE")
  
  # Create a table with the results of the metrics
  metrics_table <- data.frame(
    Test = test,
    Level = c("low", "medium", "high", "all"),
    Precision = c(metrics_low[1], metrics_medium[1], metrics_high[1], metrics_all[1]),
    Recall = c(metrics_low[2], metrics_medium[2], metrics_high[2], metrics_all[2]),
    F1_Score = c(metrics_low[3], metrics_medium[3], metrics_high[3], metrics_all[3]),
    Accuracy = c(metrics_low[4], metrics_medium[4], metrics_high[4], metrics_all[4]),
    FDR = c(metrics_low[5], metrics_medium[5], metrics_high[5], metrics_all[5]),
    AUroc = c(metrics_low[6], metrics_medium[6], metrics_high[6], metrics_all[6])
  )
  
  results_list <- list(DEGs = founded_filtered_df, Metrics = metrics_table)
  return(results_list)
}

# Heatmap for all tests passes through DE_metrics
Metrics_heatmap = function(results_list, outdir, text = NULL){
  for (test in names(results_list)) {
    # Transformar la matriz en long format
    print(class(results_list[[test]]$Metrics))
    df = reshape2::melt(results_list[[test]]$Metrics)
    colnames(df) = c("Test", "Level", "Metric", "Value")
    df$Value[df$Metric=="FDR"] = 1-df$Value[df$Metric=="FDR"]
    df$Value_discrete = ifelse(df$Value<0.3,1,ifelse(df$Value>=0.7,3,df$Value))
    df$Value_discrete = ifelse(0.3<=df$Value & df$Value<0.7,2,df$Value_discrete)
    df$Value_discrete <- factor(df$Value_discrete, levels = c("1","2","3"))
    df$Test_level = paste(df$Test, df$Level, sep="_")
    
    # Gráfico
    heatmap_test_metrics = ggplot(df, aes(x = Metric, y = Test_level, fill = Value_discrete)) +
      geom_tile() +
      ggtitle("Heatmap test metrics") +
      theme(plot.title = element_text(hjust = 0.5),
            text = element_text(size=13),
            axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
      coord_fixed() +
      scale_fill_manual(breaks = levels(df$Value_discrete),
                        values = c("#FF99BF", "#FFFF99", "#A5EDFF"),
                        labels=c('Poor','Intermediate','Good'))
    print(heatmap_test_metrics)
    ggsave(filename = paste0(outdir,"/Heatmap_test_mestric",text,"_",test,".jpg"), plot = heatmap_test_metrics)
  }
}

