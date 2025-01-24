
library(Seurat)
library(sctransform)
#library(Matrix)
library(RColorBrewer) # colorRampPalette(), brewer.pal
library(ggplot2) # qplot(), position_nudge(), geom_text()
library(cowplot) # for plot_grid
library(gplots) # for heatmap2
library(dplyr) # for mutate, top_n# Setup the Seurat Object
library(parallel)

stdin <- file('stdin', 'r')

GENOTYPE_COLORS <- readRDS(file='./datasets/genotype_colors.rds', refhook=NULL)
datasets <- list()

ln <- readLines(stdin, n=1)
while (ln != 'quit') {
	cmd <- strsplit(ln, ' ')[[1]]

	if (length(cmd) > 1) {
		opid <- cmd[1]
		
		if (cmd[2] == 'load') {
			ds <- cmd[3]

			datasets[[ds]] <- readRDS(file=sprintf('./datasets/%s/data.rds', ds), refhook=NULL)

			write(sprintf('ack %s', opid), stdout())
		} else if (cmd[2] == 'render') {
			ds <- cmd[3]
			gene <- cmd[4]
			groupBy <- cmd[5]
			splitBy <- cmd[6]

			if (groupBy == 'Genotype') {
				colors <- GENOTYPE_COLORS
			} else {
				colors <- readRDS(file=sprintf('./datasets/%s/colors.rds', ds), refhook=NULL)
			}
			
			png(sprintf('./datasets/%s/%s_umap.png', ds, opid),
				width = 5*300,        # 5 x 300 pixels
				height = 4*300,
				res = 300,            # 300 pixels per inch
				pointsize = 12)        # smaller font size
			umap <- DimPlot(datasets[[ds]], reduction="umap", label=FALSE, group.by=groupBy, cols=colors)
			print(umap)
			dev.off()

			png(sprintf("./datasets/%s/%s_vln.png", ds, opid), width=13 * 300, height=13 * 300, res=300, pointsize=4)  
			vln <- VlnPlot(
				datasets[[ds]],
				assay="RNA",
				features=c(gene),
				pt.size=0,
				split.by=splitBy,
				group.by=groupBy,
				cols=colors,
				ncol=4
			)
			print(vln)
			dev.off()

			write(sprintf('ack %s', opid), stdout())
		} else if (cmd[2] == 'unload') {
			ds <- cmd[3]

			datasets[[ds]] <- NULL

			write(sprintf('ack %s', opid), stdout())
		} else {
			write('Malformed command', stderr())
		}
	} else {
		write('Malformed command', stderr())
	}
	

	ln <- readLines(stdin, n=1)
}

