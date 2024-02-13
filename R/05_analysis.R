### Diferential Expression Analysis
# Boxplot
ggplot(as.data.frame(colData(rse_gene_SRP127581)), aes(y = assigned_gene_prop, x = sra_attribute.time )) +
  geom_boxplot() +
  theme_bw(base_size = 20) +
  ylab("Assigned Gene Prop") +
  xlab("Time")

### iSEE(rse_gene_SRP127581)


## Modelo Estadistico las lecturas de las cuentas a logaritmo base 2 de las cuentas por millón.
design <- model.matrix(~0 + sra_attribute.time + assigned_gene_prop, data = colData(rse_gene_SRP127581))
#design
## Estimacion de la varianza promedio
vGene <- voom(dge, design, plot = TRUE)

# Aplicar eBayes para obtener estadísticas bayesianas
fit_bayes <- eBayes(lmFit(vGene))
plotMA(fit_bayes, coef = 4, col = "steelblue")
results <- topTable(fit_bayes, coef = "sra_attribute.time12:00", number = nrow(rse_gene_SRP127581), sort.by = "none")
results

# visualizacion
volcanoplot(fit_bayes, coef = "sra_attribute.time12:00", highlight = 5, names = results$gene_name, col = "steelblue", hl.col="hotpink4")

top30_genes <- head(results[order(abs(results$logFC), decreasing = TRUE), ], 30)$gene_id

# exprs_top30 <- vGene$E[intersect(top30_genes, rownames(vGene$E)), ]
exprs_heatmap <-vGene$E[intersect(top30_genes, rownames(vGene$E)), ]
nombres <- rownames(results)
#nombres
rownames(exprs_heatmap) <- results$gene_name[match(rownames(exprs_heatmap), nombres)]
df <- as.data.frame(colData(rse_gene_SRP127581)[, c("sra_attribute.time","assigned_gene_prop" )])
colnames(df) <- c("time", "gene_prop")
df
pheatmap(
  exprs_heatmap,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = TRUE,
  show_colnames = FALSE,
  annotation_col = df
)