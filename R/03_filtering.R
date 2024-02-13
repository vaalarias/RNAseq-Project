### Procesamiento de los datos y formateo

## Procesar la información como la provee SRA, expandirla, separar los atributos e incorporarlos al data frame. . 
rse_gene_SRP127581 <- expand_sra_attributes(rse_gene_SRP127581)

colData(rse_gene_SRP127581)[
  ,
  grepl("^sra_attribute", colnames(colData(rse_gene_SRP127581)))
]

## Cambiar el tipo de dato a factor 
rse_gene_SRP127581$sra_attribute.time <- factor(rse_gene_SRP127581$sra_attribute.time)
rse_gene_SRP127581$sra_attribute.source_name <- factor(rse_gene_SRP127581$sra_attribute.source_name)
rse_gene_SRP127581$sra_attribute.tissue <- factor(rse_gene_SRP127581$sra_attribute.tissue)

### Filtrado 

## proporción de lecturas a genes
rse_gene_SRP127581$assigned_gene_prop <- rse_gene_SRP127581$recount_qc.gene_fc_count_all.assigned / rse_gene_SRP127581$recount_qc.gene_fc_count_all.total

## Resumen por grupo de tiempo
with(colData(rse_gene_SRP127581), tapply(assigned_gene_prop, sra_attribute.time, summary))

## Visalizar calidad de manera gráfica
with(colData(rse_gene_SRP127581), plot(assigned_gene_prop, sra_attribute.time))
abline(v = 0.35, col = "deepskyblue4", lwd = 4, lty = "solid")
## Histrograma
hist(rse_gene_SRP127581$assigned_gene_prop, col="gray")
abline(v=0.35,col="purple", lwd=7, lty = "dashed")

## copia de seguridad
rse_gene_SRP127581_copia<-rse_gene_SRP127581

## Elección de un valor de corte y eliminación de muestras con proporción menor a 0.35
table(rse_gene_SRP127581$assigned_gene_prop < 0.35)

## Eliminar las muestras con proporción menor a 0.35
rse_gene_SRP127581 <- rse_gene_SRP127581[, rse_gene_SRP127581$assigned_gene_prop > 0.35]

## Calcular niveles medios de expresión y eliminar genes poco significativos
gene_means <- rowMeans(assay(rse_gene_SRP127581, "counts"))
summary(gene_means)

# Eliminar genes poco significativos 
rse_gene_SRP127581 <- rse_gene_SRP127581[gene_means > 0.1, ]

# Dimensiones finales
dim(rse_gene_SRP127581)
# Porcentaje de genes retenidos.
round(nrow(rse_gene_SRP127581) / nrow(rse_gene_SRP127581_copia) * 100, 2)

####rse_gene_SRP127581<-rse_gene_SRP127581_copia
