## Descarga del estudio 
## estudio a utilizar SRP127581
## Proyectos disponibles
human_projects <- available_projects()
## Subset
rse_gene_SRP127581 <- create_rse(
  subset(
    human_projects,
    project == "SRP127581" & project_type == "data_sources"
  )
)
## Convertir las cuentas crudas a cuentas por lectura
assay(rse_gene_SRP127581, "counts") <- compute_read_counts(rse_gene_SRP127581)
## Ver 454 muestras
#rse_gene_SRP127581$sra.sample_attributes

