## Normalizaciooon con edgeR
dge <- DGEList(
  counts = assay(rse_gene_SRP127581, "counts"),
  genes = rowData(rse_gene_SRP127581)
)
dge <- calcNormFactors(dge)
