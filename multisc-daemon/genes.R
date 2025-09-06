#!/usr/bin/env Rscript

library(SeuratObject)
rds <- readRDS(file = "./data.rds", refhook = NULL)
genes <- rownames(rds@assays$RNA@counts)
jsonlite::write_json(genes, "./genes.json")
