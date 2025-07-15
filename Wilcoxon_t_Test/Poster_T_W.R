# Función para realizar el análisis diferencial de expresión
perform_DE_analysis <- function(seurat_obj, test_type = "wilcox", logfc_threshold = log2(1.5)) {
  comparison_result <- FindAllMarkers(seurat_obj, test.use = test_type, only.pos = TRUE, logfc.threshold = logfc_threshold)
  
  summary_over_df <- comparison_result %>%
    filter(p_val_adj < 0.05) %>%
    mutate(Gene = rownames(comparison_result)[1:n()]) %>%
    select(Gene, cluster, avg_log2FC) %>%
    rename(DEFacGroup = avg_log2FC) %>%
    mutate(logfc_level = case_when(
      DEFacGroup >= 1.2 & DEFacGroup < 2 ~ "bajo",
      DEFacGroup >= 2 & DEFacGroup < 3 ~ "medio",
      DEFacGroup >= 3 ~ "alto",
      TRUE ~ "otros"
    ))
  
  return(summary_over_df)
}

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

# Análisis de expresión diferencial para ambos métodos
tests <- c("wilcox", "t")
results_list <- list()

for (test in tests) {
  # Ejecutar el análisis diferencial
  summary_df <- perform_DE_analysis(data_1_cl, test_type = test)
  
  # Filtrar DEGs según los niveles de logfc
  DEGs_bajo <- summary_df %>% filter(logfc_level == "bajo")
  DEGs_medio <- summary_df %>% filter(logfc_level == "medio")
  DEGs_alto <- summary_df %>% filter(logfc_level == "alto")
  
  # Contar el número de genes en cada categoría
  num_bajo <- nrow(DEGs_bajo)
  num_medio <- nrow(DEGs_medio)
  num_alto <- nrow(DEGs_alto)
  
  # Calcular las métricas (usa tus datos reales para `actual` y `predicted`)
  set.seed(42)
  actual_bajo <- sample(c(0, 1), size = num_bajo, replace = TRUE)
  predicted_bajo <- sample(c(0, 1), size = num_bajo, replace = TRUE)
  
  actual_medio <- sample(c(0, 1), size = num_medio, replace = TRUE)
  predicted_medio <- sample(c(0, 1), size = num_medio, replace = TRUE)
  
  actual_alto <- sample(c(0, 1), size = num_alto, replace = TRUE)
  predicted_alto <- sample(c(0, 1), size = num_alto, replace = TRUE)
  
  metrics_bajo <- calculate_metrics(predicted_bajo, actual_bajo)
  metrics_medio <- calculate_metrics(predicted_medio, actual_medio)
  metrics_alto <- calculate_metrics(predicted_alto, actual_alto)
  
  # Crear una tabla con los resultados de las métricas
  metrics_table <- data.frame(
    Test = test,
    Level = c("Bajo", "Medio", "Alto"),
    Precision = c(metrics_bajo[1], metrics_medio[1], metrics_alto[1]),
    Recall = c(metrics_bajo[2], metrics_medio[2], metrics_alto[2]),
    F1_Score = c(metrics_bajo[3], metrics_medio[3], metrics_alto[3]),
    Accuracy = c(metrics_bajo[4], metrics_medio[4], metrics_alto[4]),
    FDR = c(metrics_bajo[5], metrics_medio[5], metrics_alto[5])
  )
  
  results_list[[test]] <- list(DEGs = summary_df, Metrics = metrics_table)
}

# Guardar los resultados en un archivo Excel
wb <- createWorkbook()

for (test in names(results_list)) {
  DEGs <- results_list[[test]]$DEGs
  metrics <- results_list[[test]]$Metrics
  
  addWorksheet(wb, paste(test, "DEGs"))
  writeData(wb, sheet = paste(test, "DEGs"), DEGs)
  
  addWorksheet(wb, paste(test, "Metrics"))
  writeData(wb, sheet = paste(test, "Metrics"), metrics)
}

saveWorkbook(wb, file = "DEGs_and_Metrics_Wilcoxon_T_test.xlsx", overwrite = TRUE)

# Mostrar resultados
for (test in names(results_list)) {
  cat("Resultados del test:", test, "\n")
  print(head(results_list[[test]]$DEGs))
  print(results_list[[test]]$Metrics)
}

