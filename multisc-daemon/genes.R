#!/usr/bin/env Rscript

library(SeuratObject)

write_genes <- function(rds = readRDS(file = "data.rds"), path = "genes.json") {
  genes <- rownames(rds@assays$RNA@counts)
  jsonlite::write_json(genes, path)
  unique(genes)
}
