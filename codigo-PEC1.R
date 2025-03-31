## -------------------------------------------------------------------------------------------------------------------------

library(POMA)
library(dplyr)

## -------------------------------------------------------------------------------------------------------------------------
# Copiar la URL y descargar el dataset
url<- "https://raw.githubusercontent.com/nutrimetabolomics/metaboData/refs/heads/main/Datasets/2024-Cachexia/human_cachexia.csv"
download.file(url, destfile = "human_cachexia.csv", method = "auto",quiet = TRUE)
# Guardar en un dataframe para poder trabajar con el
cachexia<-read.csv("human_cachexia.csv", sep =",")
# Guardar datos del ID del paciente y del grupo al que pertenecen
metadata<-data.frame(cachexia[,c(1,2)])
# Guardar los valores de todos los metabolitos
features<-data.frame(cachexia[,-c(1,2)])
#creamos el SummarizedExperiment
object<-PomaCreateObject(metada = metadata,
                         features = features)
#Lo guardamos como Rda
save(object, file = "cachexia_summarized_experiment.Rda")


## -------------------------------------------------------------------------------------------------------------------------
library(ggplot2, quiet = TRUE) # para
library(patchwork, quiet = TRUE) #para añadir los títulos
#boxplot para observar los datos en crudo
#eliminamos el nombre de los features para que se vea más claro
grafica_boxplot<-PomaBoxplots(object, x = "features",theme_params = list(axistext = "y")) + plot_annotation(
  title = "Boxplot de datos crudos", theme = theme(plot.title = element_text(hjust = 0.5, size = 13)))



## -------------------------------------------------------------------------------------------------------------------------
#normalización pareto: transformación logarítmica y escalado de los datos 
normalized<-PomaNorm(object, method="log_pareto")
#boxplot para observar los datos normalizados
boxplot_norm<-PomaBoxplots(normalized, x = "features",theme_params = list(axistext = "y"))+ plot_annotation(
  title = "Boxplot de los datos normalizados", theme = theme(plot.title = element_text(hjust = 0.5, size = 13)))


## -------------------------------------------------------------------------------------------------------------------------
#buscar outliers en nuestros datos normalizados 
poma_outliers<-PomaOutliers(normalized)
#representación gráfica 
grafica_outliers<-poma_outliers$polygon_plot + plot_annotation(
  title = "Búsqueda de outliers", theme = theme(plot.title = element_text(hjust = 0.5, size = 13)))
#visualizacion de los outliers (son 3)
outliers<-poma_outliers$outliers
## -------------------------------------------------------------------------------------------------------------------------
#guardamos el conjunto de datos limpio para trabajar con el 
poma_processed<-poma_outliers$data


## ----include = FALSE------------------------------------------------------------------------------------------------------
#pca de los datos procesados, sin escalar (tiene escalado pareto)
#centramos los datos, pero no escalamos 
poma_pca<-PomaPCA(poma_processed, center = TRUE, scale = FALSE, ellipse = TRUE)
#representación gráfica bidimensional de la PCA
grafica_pca<-poma_pca$factors_plot + plot_annotation(
  title = "Representación de la PCA", theme = theme(plot.title = element_text(hjust = 0.5, size = 13)))
#vamos a escoger los 5 metabolitos que más contribuyen a explicar la varianza total
load<-as.data.frame(poma_pca$loadings)
load5<-head(load[order(load[,"PC1"],decreasing=TRUE),],5)


## -------------------------------------------------------------------------------------------------------------------------
#test mann-whitney para observar diferencias signficativas entre grupos para cada metabolito 
#especifico la variable grupo, ya que la primera es el ID del paciente
poma_univariate<-PomaUnivariate(poma_processed, method = "mann", covs = "Muscle.loss")
#guardamos los resultados de los que tienen una diferencia signficativa 
dif_sig<-poma_univariate[poma_univariate$adj_pvalue<0.05,]
#guardamos los datos de los metabolitos signficativos
#para el analisis multivariante posterior
datos_sig<-poma_processed[dif_sig$feature,]


## -------------------------------------------------------------------------------------------------------------------------
#Volcano plot para observar las diferencias individuales 
grafica_volcano<-poma_univariate %>% select(feature, fold_change,pvalue) %>%
  PomaVolcano (pval_cutoff = 0.05, labels = TRUE)

## -------------------------------------------------------------------------------------------------------------------------
#sPLSA-DA para intentar maximizar la separación de los grupos 
#Utilizamos los metabolitos con diferencias significativas
pls<-PomaPLS(datos_sig, method ="splsda", y ="Muscle.loss", ellipse = TRUE)
#grafica de las cargas
grafica_pls_features<-pls$selected_features_plot
#gráfico de factores 
grafica_pls_factors<-pls$factors_plot
# Guardar los 5 metabolitos que mejor separan los grupos 
pls_features<-as.data.frame(pls$selected_features)
feat5<-head(pls_features[order(abs(pls_features[,"value"]),decreasing=TRUE),],5)
