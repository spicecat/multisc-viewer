
start <- Sys.time()
library(Seurat)
end <- Sys.time()
importSeurat <- end - start

start <- Sys.time()
library(sctransform)
#library(Matrix)
library(RColorBrewer) # colorRampPalette(), brewer.pal
library(ggplot2) # qplot(), position_nudge(), geom_text()
library(cowplot) # for plot_grid
library(gplots) # for heatmap2
library(dplyr) # for mutate, top_n# Setup the Seurat Object
end <- Sys.time()
importOther <- end - start

args = commandArgs(trailingOnly = TRUE)

if (length(args) < 5) {
	stop("Invalid args, usage: Rscript plot.r <gene> <groupBy> <splitBy> <clustering file> <violin file>")
}

start <- Sys.time()
rds <- readRDS(file="./data.rds", refhook=NULL)
end <- Sys.time()
rdsRead <- end - start

if (args[2] == 'CellType') {
	colors <- readRDS(file="./colors.rds", refhook=NULL)
} else {
	colors <- readRDS(file="../genotype_colors.rds", refhook=NULL)
}

# c("coral3", "deepskyblue3", "goldenrod2", "green3", "red3", "pink3")

rds[["integrated"]] <- NULL
rds[["SCT"]] <- NULL

start <- Sys.time()
png(args[4],
    width = 5*300,        # 5 x 300 pixels
    height = 4*300,
    res = 300,            # 300 pixels per inch
    pointsize = 12)        # smaller font size
DimPlot(rds, reduction="umap", label=FALSE, group.by=args[2], cols=colors)
end <- Sys.time()
clusterPlot <- end - start

marker_genes <- c(args[1])

start <- Sys.time()
png(args[5], width=13 * 300, height=13 * 300, res=300, pointsize=4)  
VlnPlot(
	rds,
	assay="RNA",
	features=marker_genes,
	pt.size=0,
	split.by=args[3],
	group.by=args[2],
	cols=colors,
	ncol=4
)
end <- Sys.time()
violinPlot <- end - start

times <- c(importSeurat, importOther, rdsRead, clusterPlot, violinPlot)
write(times, sep=" ", file="../times.dump", append=TRUE)