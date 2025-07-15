Resumen de lo realizado con Metacélulas:

- En el archivo "Metacells.R" inicio con la simulación "simulation_2.rds". Esa simulación la convierto a tipo SEURAT para hacer el preprocesamiento con el pipeline de Seurat.
- Luego, calculo los DEGs esperados "DEGs_expected.csv", esos sirven para la comparación de métodos normal como para la comparación de métodos con metacélulas.
- A continuación realizo la construcción de metacélulas, junto con las métricas de su construcción. Además realizo el preprocesamiento, el archivo donde se puede encontrar es "MC_gamma_20_seurat_object.Rds".
- Luego, realicé el análisis de expresión diferencial, para el cuál utilicé las funciones de "Functions.R" y todos los resultados están guardados en formato .csv, por ejemplo: "DE_metacells_t_test_unweighted.csv".
- Finalmente, calculé las métricas y realicé los gráficos utilizando las funciones de "Functions.R", estás métricas estan guardadas en el archivo "DEGs_and_Metrics_metacells_all_tests.xlsx" y las gráficas en la carpeta "figures_metacells".