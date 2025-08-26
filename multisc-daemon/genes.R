library("Seurat")
library("rjson")

rds <- readRDS(file = "./data.rds", refhook = NULL)
genes <- rownames(rds@assays$RNA@counts)

write(rjson::toJSON(as.list(genes)), "./genes.json")
