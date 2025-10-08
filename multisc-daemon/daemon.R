suppressPackageStartupMessages({
  library(SeuratObject)
  source("genes.R")
})

# --- Configuration ---
ds_dir <- Sys.getenv("DATASETS_DIR", "../data/datasets")
ds_index_file <- "index.json"
ds_meta_file <- "metadata.json"
ds_data_file <- "data.rds"
ds_genes_file <- "genes.json"
ds_degs_file <- "degs.json"
genotype_color_file <- "genotype.colors.rds"
cluster_color_file <- "cluster.colors.rds"

pub_dir <- Sys.getenv("PUBLICATIONS_DIR", "../data/publications")
pub_index_file <- "index.json"
pub_meta_file <- "metadata.json"

plots_dir <- Sys.getenv("PLOTS_DIR", "../data/plots")
plot_file <- "plot.png"
dir.create(plots_dir, FALSE, TRUE)

# --- Global State ---
ds_data <- list()

# --- Helpers ---
get_ds_path <- function() {
  jsonlite::fromJSON(file.path(ds_dir, ds_index_file))
}

get_loaded <- function() {
  n <- names(ds_data)
  if (is.null(n)) {
    list(datasets = list())
  } else {
    list(datasets = n)
  }
}

read_ds_data <- function(ds_path) {
  readRDS(file.path(ds_dir, ds_path, ds_data_file))
}

get_colors <- function(ds_path, group_by) {
  color_file <- ifelse(
    group_by == "Genotype",
    genotype_color_file,
    cluster_color_file
  )
  readRDS(file.path(ds_dir, ds_path, color_file))
}

save_plot <- function(plot_object, path, ...) {
  temp_path <- tempfile(tmpdir = dirname(path))
  png(temp_path, ...)
  print(plot_object)
  dev.off()
  file.rename(temp_path, path)
}

render_plot <- function(ds_path, ds, gene, pt, group_by, split_by) {
  data <- ds_data[[ds]]
  colors <- get_colors(ds_path, group_by)
  plot_id <- paste(ds, gene, pt,
    paste(group_by, split_by, sep = ":"),
    sep = "/"
  )
  plot_path <- file.path(
    plots_dir,
    plot_id,
    plot_file
  )
  dir.create(dirname(plot_path), recursive = TRUE, showWarnings = FALSE)
  print(paste("Rendering plot:", plot_id))
  if (pt == "umap") {
    umap_plot <- Seurat::DimPlot(
      data,
      reduction = "umap",
      label = FALSE,
      group.by = group_by,
      cols = colors
    )
    save_plot(umap_plot, plot_path,
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

#* @apiTitle MultiSC-Daemon API

#* Log request information
#* @filter logger
function(req) {
  cat(
    as.character(Sys.time()), "-", req$REQUEST_METHOD, req$PATH_INFO, "-",
    req$HTTP_USER_AGENT, "@", req$REMOTE_ADDR, "\n"
  )
  plumber::forward()
}

#* Check daemon health.
#* @serializer unboxedJSON
#* @get /health
function() {
  list(status = "ok")
}

#* Get loaded datasets on daemon.
#* @get /loaded
function() {
  get_loaded()
}

#* Unload datasets on daemon.
#* @post /unload
#* @param ds:[str] Dataset ids. Unset to unload all datasets.
function(ds) {
  if (missing(ds)) ds <- names(ds_data)
  ds_data <<- ds_data[!names(ds_data) %in% ds]
  get_loaded()
}

#* Load datasets on daemon.
#* @serializer unboxedJSON
#* @post /load
#* @param ds:[str]* Dataset ids
load_ds <- function(ds) {
  ds_path <- get_ds_path()
  lapply(ds, function(d) {
    try({
      if (!d %in% names(ds_data)) ds_data[[d]] <<- read_ds_data(ds_path[[d]])
    })
  })
  get_loaded()
}

#* Get datasets metadata.
#* @serializer unboxedJSON
#* @get /datasets
#* @param ds:[str] Dataset ids. Unset to get all datasets.
function(ds) {
  ds_path <- get_ds_path()
  if (missing(ds)) ds <- names(ds_path)
  metadata <- setNames(lapply(ds, function(d) {
    tryCatch(
      {
        meta <- jsonlite::fromJSON(file.path(
          ds_dir, ds_path[[d]], ds_meta_file
        ))
        meta$size <- file.info(file.path(
          ds_dir, ds_path[[d]], ds_data_file
        ))$size
        meta
      },
      error = function(e) NULL
    )
  }), ds)
  metadata[!vapply(metadata, is.null, logical(1))]
}

#* Get publications metadata.
#* @serializer unboxedJSON
#* @get /publications
#* @param pub:[str] Publication ids. Unset to get all publications.
function(pub) {
  pub_path <- jsonlite::fromJSON(file.path(pub_dir, pub_index_file))
  if (missing(pub)) pub <- names(pub_path)
  metadata <- setNames(lapply(pub, function(p) {
    tryCatch(
      jsonlite::fromJSON(file.path(pub_dir, pub_path[[p]], pub_meta_file)),
      error = function(e) NULL
    )
  }), pub)
  metadata[!vapply(metadata, is.null, logical(1))]
}

#* Get genes for datasets.
#* @get /genes
#* @param ds:[str]* Dataset ids
function(ds) {
  ds_path <- get_ds_path()
  setNames(lapply(ds, function(d) {
    tryCatch(
      jsonlite::fromJSON(file.path(ds_dir, ds_path[[d]], ds_genes_file)),
      error = function(e) {
        tryCatch(
          write_genes(read_ds_data(ds_path[[d]])),
          error = function(e) {
            list()
          }
        )
      }
    )
  }), ds)
}

#* Get differentially expressed genes for datasets.
#* @get /degs
#* @param ds:[str]* Dataset ids
function(ds) {
  ds_path <- get_ds_path()
  setNames(lapply(ds, function(d) {
    tryCatch(
      unique(jsonlite::fromJSON(file.path(ds_dir, ds_path[[d]], ds_degs_file))),
      error = function(e) {
        list()
      }
    )
  }), ds)
}

#* Generate plots for datasets.
#* @serializer unboxedJSON
#* @post /plots
#* @param datasets:[str]* Dataset ids
#* @param gene:[str]* Genes
#* @param plotTypes:[str] Plots to render
#* @param groupBy:str Grouping variable
#* @param splitBy:str Splitting variable
function(
    ds, gene, pt = list("umap", "vln", "feat"),
    groupBy = "CellType", splitBy = "Genotype") {
  ds_path <- get_ds_path()
  unlist(lapply(ds, function(d) {
    try({
      if (!d %in% names(ds_data)) ds_data[[d]] <<- read_ds_data(ds_path[[d]])
    })
    unlist(lapply(gene, function(g) {
      lapply(pt, function(p) {
        render_plot(ds_path[[d]], d, g, p, groupBy, splitBy)
      })
    }))
  }))
}
