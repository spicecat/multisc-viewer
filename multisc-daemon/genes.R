library(Seurat)
library(jsonlite)

rds <- readRDS(file = "./data.rds", refhook = NULL)
genes <- rownames(rds@assays$RNA@counts)

write_json(genes, "./genes.json")
