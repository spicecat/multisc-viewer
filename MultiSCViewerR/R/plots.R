#' @importFrom ggplot2 ggsave

color_file <- "colors.json"
plot_file <- "plot.png"

#' Save Plot to File
#'
#' Saves a ggplot object as PNG.
#'
#' @param plot_object A ggplot object to save
#' @param path Destination file path for the PNG output
#' @param ... Additional arguments passed to [ggplot2::ggsave()]
#' @keywords internal
save_plot <- function(plot_object, path, ...) {
  temp_path <- tempfile(tmpdir = dirname(path))
  ggplot2::ggsave(
    temp_path, plot_object,
    device = "png", units = "px", ...,
    create.dir = TRUE
  )
  file.rename(temp_path, path)
}

#' Render Plot for Dataset and Gene
#'
#' Creates visualization plots (UMAP, violin, or feature plots) for specified
#' dataset and gene combinations.
#'
#' @param ds Dataset ID to visualize
#' @param gene Gene name to plot
#' @param pt Plot type: "umap" (UMAP projection), "vln" (violin),
#'   or "feat" (feature plot)
#' @param group_by Column name for grouping/coloring cells
#' @param split_by Column name for splitting the plot
#' @return Plot ID string (format: "ds/gene/pt/groupBy:splitBy"),
#'   or `NULL` if rendering fails
#' @keywords internal
render_plot <- function(ds, gene, pt, group_by, split_by) {
  plot_id <- paste(ds, gene, pt,
    paste(group_by, split_by, sep = ":"),
    sep = "/"
  )
  plot_path <- file.path(
    .env$plots_dir,
    plot_id,
    plot_file
  )
  print(paste("Rendering plot:", plot_id, plot_path))
  data <- .ds_data[[ds]]
  colors <- jsonlite::fromJSON(file.path(
    .env$ds_index[[ds]], color_file
  ))[[group_by]]

  if (pt == "umap") {
    p <- Seurat::DimPlot(
      data,
      reduction = "umap",
      label = FALSE,
      group.by = group_by,
      cols = colors
    )
    save_plot(p, plot_path,
      width = 5 * 300, height = 4 * 300, dpi = 300, pointsize = 12
    )
  } else if (pt == "vln") {
    p <- Seurat::VlnPlot(data,
      assay = "RNA", features = gene, pt.size = 0,
      split.by = group_by, group.by = split_by, cols = colors, ncol = 1
    )
    save_plot(p, file.path(plot_path),
      width = 1200, height = 900, dpi = 300
    )
  } else if (pt == "feat") {
    p <- Seurat::FeaturePlot(data,
      features = gene, min.cutoff = "q5", max.cutoff = "q95"
    )
    save_plot(p, file.path(plot_path),
      width = 1200, height = 1200, dpi = 300
    )
  } else {
    return()
  }
  plot_id
}

#' Generate Plots for Datasets
#'
#' Generates plots based on query parameters. Loads datasets and renders
#' visualizations for each gene and plot type combination.
#'
#' @param query List with required elements:
#'   - `ds`: Character vector of dataset IDs
#'   - `gene`: Character vector of gene names
#'   - `pt`: Character vector of plot types ("umap", "vln", "feat")
#'   - `groupBy`: Grouping column name
#'   - `splitBy`: Splitting column name
#' @return Character vector of successfully rendered plot IDs,
#'   or `NULL` values for failed plots
#' @keywords internal
plots <- function(query) {
  load_ds(query)
  unlist(lapply(query$ds, function(d) {
    unlist(lapply(query$gene, function(g) {
      unlist(lapply(query$pt, function(p) {
        tryCatch(
          render_plot(d, g, p, query$groupBy, query$splitBy),
          error = function(e) NULL
        )
      }))
    }))
  }))
}
