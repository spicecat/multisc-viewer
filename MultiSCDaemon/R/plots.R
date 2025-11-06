genotype_color_file <- "genotype.colors.rds"
cluster_color_file <- "cluster.colors.rds"
plot_file <- "plot.png"

get_colors <- function(ds_path, group_by) {
  color_file <- ifelse(
    group_by == "Genotype",
    genotype_color_file,
    cluster_color_file
  )
  readRDS(file.path(ds_path, color_file))
}

save_plot <- function(plot_object, path, ...) {
  temp_path <- tempfile(tmpdir = dirname(path))
  png(temp_path, ...)
  print(plot_object)
  dev.off()
  file.rename(temp_path, path)
}


render_plot <- function(ds_path, ds, gene, pt, group_by, split_by) {
  data <- .ds_data[[ds]]
  colors <- get_colors(ds_path, group_by)
  plot_id <- paste(ds, gene, pt,
    paste(group_by, split_by, sep = ":"),
    sep = "/"
  )
  plot_path <- file.path(
    .config$plots_dir,
    plot_id,
    plot_file
  )
  dir.create(dirname(plot_path), recursive = TRUE, showWarnings = FALSE)
  print(paste("Rendering plot:", plot_id))
  if (pt == "umap") {
    p <- Seurat::DimPlot(
      data,
      reduction = "umap",
      label = FALSE,
      group.by = group_by,
      cols = colors
    )
    save_plot(p, plot_path,
      width = 5 * 300, height = 4 * 300, res = 300, pointsize = 12
    )
  } else if (pt == "vln") {
    p <- Seurat::VlnPlot(data,
      assay = "RNA", features = gene, pt.size = 0,
      split.by = group_by, group.by = split_by, cols = colors, ncol = 1
    )
    save_plot(p, file.path(plot_path),
      width = 1200, height = 900, res = 300
    )
  } else if (pt == "feat") {
    p <- Seurat::FeaturePlot(data,
      features = gene, min.cutoff = "q5", max.cutoff = "q95"
    )
    save_plot(p, file.path(plot_path),
      width = 1200, height = 1200, res = 300
    )
  } else {
    stop(paste("Unknown plot type:", pt))
  }
  plot_id
}

plots <- function(
  ds, gene, pt = list("umap", "vln", "feat"),
  groupBy = "CellType", splitBy = "Genotype"
) {
  ds_index <- get_ds_index()
  unlist(lapply(ds, function(d) {
    try({
      if (!d %in% ls(.ds_data)) {
        .ds_data[[d]] <- get_ds_data(ds_index[[d]])
      }
    })
    unlist(lapply(gene, function(g) {
      lapply(pt, function(p) {
        render_plot(ds_index[[d]], d, g, p, groupBy, splitBy)
      })
    }))
  }))
}
