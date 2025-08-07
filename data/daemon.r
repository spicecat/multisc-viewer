suppressPackageStartupMessages({
  library(plumber)
  library(Seurat)
  library(ggplot2)
  library(future)
  library(future.apply)
  library(promises)
})

# --- Plumber Configuration ---
# workers <- strtoi(Sys.getenv("WORKERS", availableCores()))
# cat(paste("plumber using", workers, "workers\n"))
plan(multisession)
# future::plan(future::multisession(workers = workers))

globals_max_size <- strtoi(Sys.getenv(
  "GLOBALS_MAX_SIZE",
  8
)) * 1024^3 # GB
options(future.globals.maxSize = globals_max_size)

# --- Data Configuration ---
datasets_dir <- Sys.getenv("DATASETS_DIR", "datasets")
cache_dir <- Sys.getenv("PLOT_CACHE_DIR", "plots")
if (!dir.exists(cache_dir)) {
  dir.create(cache_dir)
}
data_file <- Sys.getenv("DATA_FILE", "data.rds")
genotype_color_file <- Sys.getenv(
  "GENOTYPE_COLOR_FILE",
  "genotype.colors.rds"
)
cluster_color_file <- Sys.getenv(
  "CLUSTER_COLOR_FILE",
  "cluster.colors.rds"
)

# --- Global State ---
loaded_datasets <- list()

# --- Helpers ---
get_datasets <- function() {
  list(datasets = names(loaded_datasets))
}

load_dataset <- function(ds) {
  if (!ds %in% names(loaded_datasets)) {
    loaded_datasets[[ds]] <<- readRDS(file.path(datasets_dir, ds, data_file))
    ds
  }
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
  png(path, ...)
  print(plot_object)
  dev.off()
}

#* @apiTitle MultiSC Plotting Daemon

#* Log request information
#* @filter logger
function(req) {
  cat(
    as.character(Sys.time()), "-", req$REQUEST_METHOD, req$PATH_INFO, "-",
    req$HTTP_USER_AGENT, "@", req$REMOTE_ADDR, "\n"
  )
  forward()
}

#* Get server status
#* @get /status
function() {
  future_promise({
    get_datasets()
  })
}

#* Load datasets
#* @post /load
#* @param datasets:[str]* Dataset names to load
function(datasets) {
  lapply(datasets, load_dataset)
  get_datasets()
}

#* Unload a dataset
#* @post /unload
#* @param datasets:[str]* Dataset names to unload
function(datasets) {
  loaded_datasets <- loaded_datasets[!names(loaded_datasets) %in% datasets]
  get_datasets()
}

#* Render plots for a dataset
#* @post /render
#* @param ds:str* Dataset name
#* @param gene:str* Gene to plot
#* @param groupBy:str* Grouping variable
#* @param splitBy:str* Splitting variable
#* @param plots:[str]* Plots to render (umap, vln, feat)
function(ds, gene, groupBy, splitBy, plots = list("umap", "vln", "feat")) {
  load_dataset(ds)
  future({
    dataset <- loaded_datasets[[ds]]
    cache_key <- paste(ds, gene, groupBy, splitBy, sep = ":")
    plot_dir <- file.path(cache_dir, cache_key)
    if (!dir.exists(plot_dir)) {
      dir.create(plot_dir, recursive = TRUE)
    }
    colors <- get_colors(ds, groupBy)

    lapply(
      plots,
      function(plot_type) {
        if (plot_type == "umap") {
          # UMAP plot
          umap_plot <- DimPlot(
            dataset,
            reduction = "umap", label = FALSE, group.by = groupBy, cols = colors
          )
          save_plot(umap_plot, file.path(plot_dir, "umap.png"),
            width = 5 * 300, height = 4 * 300, res = 300, pointsize = 12
          )
        } else if (plot_type == "vln") {
          # Violin plot
          vln_plot <- VlnPlot(
            dataset,
            assay = "RNA",
            features = c(gene),
            pt.size = 0,
            split.by = groupBy,
            group.by = splitBy,
            cols = colors,
            ncol = 1
          )
          save_plot(vln_plot, file.path(plot_dir, "vln.png"),
            width = (length(unique(dataset@meta.data[[groupBy]])) + 2) * 300,
            height = 3 * 300, res = 300, pointsize = 4
          )
        } else if (plot_type == "feat") {
          # Feature plot
          feature_plot <- FeaturePlot(dataset,
            features = c(gene),
            min.cutoff = "q5", max.cutoff = "q95"
          )
          save_plot(feature_plot, file.path(plot_dir, "feat.png"),
            width = 4 * 300, height = 4 * 300, res = 300, pointsize = 5
          )
        }
      }
    )

    list(key = cache_key)
  })
}
