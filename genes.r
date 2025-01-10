library(Seurat)
library(rjson)

rds <- readRDS(file = "./data.rds", refhook = NULL)

genes <- list()
for (gene in rownames(rds@assays$RNA@counts)) {
	genes <- append(genes, gene)
}

write(toJSON(genes), "./genes.json")
