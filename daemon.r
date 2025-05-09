library(Seurat)
library(sctransform)
# library(Matrix)
library(RColorBrewer) # colorRampPalette(), brewer.pal
library(ggplot2) # qplot(), position_nudge(), geom_text()
library(cowplot) # for plot_grid
library(gplots) # for heatmap2
library(dplyr) # for mutate, top_n# Setup the Seurat Object
library(parallel)

stdin <- file("stdin", "r")

datasets <- list()

ln <- readLines(stdin, n = 1)
while (ln != "quit") {
  cmd <- strsplit(ln, " ")[[1]]

  if (length(cmd) > 1) {
    opid <- cmd[1]

    if (cmd[2] == "load") {
      ds <- cmd[3]

      datasets[[ds]] <- readRDS(
        file = sprintf("./datasets/%s/data.rds", ds), refhook = NULL
      )

      write(sprintf("ack %s", opid), stdout())
    } else if (cmd[2] == "render") {
      ds <- cmd[3]
      gene <- cmd[4]
      group_by <- cmd[5]
      split_by <- cmd[6]

      if (group_by == "Genotype") {
        colors <- readRDS(
          file = sprintf("./datasets/%s/genotype.colors.rds", ds), refhook = NULL
        )
      } else {
        colors <- readRDS(
          file = sprintf("./datasets/%s/cluster.colors.rds", ds), refhook = NULL
        )
      }

      png(sprintf("./datasets/%s/%s_umap.png", ds, opid),
        width = 5 * 300, # 5 x 300 pixels
        height = 4 * 300,
        res = 300, # 300 pixels per inch
        pointsize = 12
      ) # smaller font size
      umap <- DimPlot(
        datasets[[ds]],
        reduction = "umap", label = FALSE, group.by = group_by, cols = colors
      )
      print(umap)
      dev.off()

      png(
        sprintf("./datasets/%s/%s_vln.png", ds, opid),
        width = (dim(table(datasets[[ds]]@meta.data$Genotype)) + 2) * 300, height = 3 * 300, res = 300, pointsize = 4
      )
      vln <- VlnPlot(
        datasets[[ds]],
        assay = "RNA",
        features = c(gene),
        pt.size = 0,
        split.by = group_by,
        group.by = split_by,
        cols = colors,
        ncol = 1
      )
      print(vln)
      dev.off()

	  png(sprintf("./datasets/%s/%s_feat.png", ds, opid),     
		width = 4 *300,        # 5 x 300 pixels
		height = 4*300,
		res = 300,            # 300 pixels per inch
		pointsize = 5)        # smaller font size
	  
	  feature <- FeaturePlot(datasets[[ds]], features = c(gene), min.cutoff = "q5", max.cutoff = "q95")
	  print(feature)
	  dev.off()

      write(sprintf("ack %s", opid), stdout())
    } else if (cmd[2] == "unload") {
      ds <- cmd[3]

      datasets[[ds]] <- NULL

      write(sprintf("ack %s", opid), stdout())
    } else {
      write("Malformed command", stderr())
    }
  } else {
    write("Malformed command", stderr())
  }


  ln <- readLines(stdin, n = 1)
}
