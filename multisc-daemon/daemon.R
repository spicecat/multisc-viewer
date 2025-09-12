suppressPackageStartupMessages({
  library(SeuratObject)
})

# --- Configuration ---
datasets_dir <- Sys.getenv("DATASETS_DIR", "../data/datasets")
datasets_meta <- Sys.getenv("DATASETS_META", "meta.json")
publications_dir <- Sys.getenv("PUBLICATIONS_DIR", "../data/publications")
publications_meta <- Sys.getenv("PUBLICATIONS_META", "meta.json")
plot_dir <- Sys.getenv("PLOT_DIR", "../data/plots")
plot_file <- Sys.getenv("PLOT_FILE", "plot.png")
data_file <- Sys.getenv("DATA_FILE", "data.rds")
genotype_color_file <- Sys.getenv("GENOTYPE_COLOR_FILE", "genotype.colors.rds")
cluster_color_file <- Sys.getenv("CLUSTER_COLOR_FILE", "cluster.colors.rds")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

# --- Global State ---
dataset_data <- list()
dataset_status <- list()

# --- Helpers ---
get_datasets <- function() {
  n <- names(dataset_data)
  if (is.null(n)) {
    list(datasets = list())
  } else {
    list(datasets = n)
  }
}

load_dataset <- function(ds) {
  if (!ds %in% names(dataset_data)) {
    dataset_data[[ds]] <<- readRDS(file.path(datasets_dir, ds, data_file))
    print(paste("Loaded dataset:", ds))
  }
  TRUE
}

get_colors <- function(ds, group_by) {
  color_file <- ifelse(
    group_by == "Genotype",
    genotype_color_file,
    cluster_color_file
  )
  readRDS(file.path(datasets_dir, ds, color_file))
}

save_plot <- function(plot_object, path, ...) {
  temp_path <- tempfile(tmpdir = dirname(path))
  png(temp_path, ...)
  print(plot_object)
  dev.off()
  file.rename(temp_path, path)
}

render_plot <- function(ds, gene, group_by, split_by, pt) {
  dataset <- dataset_data[[ds]]
  colors <- get_colors(ds, group_by)
  plot_id <- paste(ds, gene,
    paste(group_by, split_by, sep = ":"), pt,
    sep = "/"
  )
  plot_path <- file.path(
    plot_dir,
    plot_id,
    plot_file
  )
  dir.create(dirname(plot_path), recursive = TRUE, showWarnings = FALSE)
  print(paste("Rendering plot:", plot_id))
  if (pt == "umap") {
    umap_plot <- Seurat::DimPlot(
      dataset,
      reduction = "umap",
      label = FALSE,
      group.by = group_by,
      cols = colors
    )
    save_plot(umap_plot, plot_path,
      width = 5 * 300, height = 4 * 300, res = 300, pointsize = 12
    )
  } else if (pt == "vln") {
    p <- Seurat::VlnPlot(dataset,
      assay = "RNA", features = gene, pt.size = 0,
      split.by = group_by, group.by = split_by, cols = colors, ncol = 1
    )
    save_plot(p, file.path(plot_path),
      width = 1200, height = 900, res = 300
    )
  } else if (pt == "feat") {
    p <- Seurat::FeaturePlot(dataset,
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

#* @apiTitle MultiSC Viewer Daemon

#* Log request information
#* @filter logger
function(req) {
  cat(
    as.character(Sys.time()), "-", req$REQUEST_METHOD, req$PATH_INFO, "-",
    req$HTTP_USER_AGENT, "@", req$REMOTE_ADDR, "\n"
  )
  plumber::forward()
}

#* Get datasets
#* @get /datasets
function() {
  jsonlite::fromJSON(file.path(datasets_dir, datasets_meta))
}

#* Get publications
#* @get /publications
function() {
  jsonlite::fromJSON(file.path(publications_dir, publications_meta))
}

#* Get loaded datasets
#* @get /status
#* @get /loaded
function() {
  get_datasets()
}

#* Unload a dataset
#* @post /unload
#* @param datasets:[str]* Dataset names to unload
function(datasets) {
  loaded_datasets <<- dataset_data[!names(dataset_data) %in% datasets]
  get_datasets()
}

#* Load datasets
#* @serializer unboxedJSON
#* @post /load
#* @param datasets:[str]* Dataset names to load
function(datasets) {
  setNames(lapply(datasets, function(ds) {
    tryCatch(
      {
        load_dataset(ds)
      },
      error = function(e) {
        FALSE
      }
    )
  }), datasets)
}

#* Render plots for datasets
#* @serializer unboxedJSON
#* @post /render
#* @param datasets:[str]* Dataset names
#* @param genes:[str]* Genes to plot
#* @param groupBy:str Grouping variable
#* @param splitBy:str Splitting variable
#* @param plotTypes:[str] Plot types to render (umap, vln, feat)
function(
    datasets, genes, groupBy = "CellType", splitBy = "Genotype",
    plotTypes = list("umap", "vln", "feat")) {
  lapply(datasets, load_dataset)
  unlist(lapply(datasets, function(ds) {
    unlist(lapply(genes, function(gene) {
      lapply(plotTypes, function(pt) {
        render_plot(ds, gene, groupBy, splitBy, pt)
      })
    }))
  }))
}
