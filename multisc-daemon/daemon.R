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
get_ds_index <- function() {
  jsonlite::fromJSON(file.path(ds_dir, ds_index_file))
}

ds_loaded <- function() {
  n <- names(ds_data)
  if (is.null(n)) {
    list(datasets = list())
  } else {
    list(datasets = n)
  }
}

get_ds_data <- function(ds_path) {
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
  temp_index <- tempfile(tmpdir = dirname(path))
  png(temp_index, ...)
  print(plot_object)
  dev.off()
  file.rename(temp_index, path)
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
  ds_loaded()
}

#* Unload datasets on daemon.
#* @post /unload
#* @param ds:[str] Dataset ids. Unset to unload all datasets.
function(ds) {
  if (missing(ds)) ds <- names(ds_data)
  ds_data <<- ds_data[!names(ds_data) %in% ds]
  ds_loaded()
}

#* Load datasets on daemon.
#* @serializer unboxedJSON
#* @post /load
#* @param ds:[str]* Dataset ids
function(ds) {
  ds_index <- get_ds_index()
  lapply(ds, function(d) {
    try({
      if (!d %in% names(ds_data)) ds_data[[d]] <<- get_ds_data(ds_index[[d]])
    })
  })
  ds_loaded()
}

#* Get datasets metadata.
#* @serializer unboxedJSON
#* @get /datasets
#* @param ds:[str] Dataset ids. Unset to get all datasets.
function(ds) {
  ds_index <- get_ds_index()
  if (missing(ds)) ds <- names(ds_index)
  metadata <- setNames(lapply(ds, function(d) {
    tryCatch(
      {
        meta <- jsonlite::fromJSON(file.path(
          ds_dir, ds_index[[d]], ds_meta_file
        ))
        meta$size <- file.info(file.path(
          ds_dir, ds_index[[d]], ds_data_file
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
  pub_index <- jsonlite::fromJSON(file.path(pub_dir, pub_index_file))
  if (missing(pub)) pub <- names(pub_index)
  metadata <- setNames(lapply(pub, function(p) {
    tryCatch(
      jsonlite::fromJSON(file.path(pub_dir, pub_index[[p]], pub_meta_file)),
      error = function(e) NULL
    )
  }), pub)
  metadata[!vapply(metadata, is.null, logical(1))]
}

#* Get genes for datasets.
#* @serializer unboxedJSON
#* @get /genes
#* @param ds:[str]* Dataset ids
function(ds) {
  ds_index <- get_ds_index()
  setNames(lapply(ds, function(d) {
    tryCatch(
      jsonlite::fromJSON(file.path(ds_dir, ds_index[[d]], ds_genes_file)),
      error = function(e) {
        tryCatch(
          write_genes(get_ds_data(ds_index[[d]])),
          error = function(e) list()
        )
      }
    )
  }), ds)
}

#* Get differentially expressed genes for datasets.
#* @serializer unboxedJSON
#* @get /degs
#* @param ds:[str]* Dataset ids
function(ds) {
  ds_index <- get_ds_index()
  degs <- setNames(lapply(ds, function(d) {
    tryCatch(
      jsonlite::fromJSON(file.path(ds_dir, ds_index[[d]], ds_degs_file)),
      error = function(e) NULL
    )
  }), ds)
  degs[!vapply(degs, is.null, logical(1))]
}

#* Get aggregated gene rows and summary metadata for datasets.
#* Computes gene rows from DEGs and returns lightweight metadata to avoid transferring large arrays.
#* @serializer unboxedJSON
#* @get /gene_rows
#* @param ds:[str]* Dataset ids
function(ds) {
  ds_index <- get_ds_index()

  # Aggregation maps
  rows_map <- new.env(parent = emptyenv())
  degs_meta <- list()
  genes_count <- list()

  lapply(ds, function(d) {
    # Read DEGs for dataset and build meta and rows
    degs <- tryCatch(
      jsonlite::fromJSON(file.path(ds_dir, ds_index[[d]], ds_degs_file)),
      error = function(e) NULL
    )
    if (!is.null(degs)) {
      # per-dataset meta container
      degs_meta[[d]] <<- list()
      # iterate each DEG group
      lapply(names(degs), function(deg_id) {
        ds_deg <- degs[[deg_id]]
        # record deg meta: name and gene count (without returning full list)
        degs_meta[[d]][[deg_id]] <<- list(
          `_id` = deg_id,
          name = if (!is.null(ds_deg$name)) ds_deg$name else NULL,
          geneCount = length(ds_deg$gene)
        )
        # aggregate rows by gene id
        lapply(ds_deg$gene, function(g) {
          if (!exists(g, envir = rows_map, inherits = FALSE)) {
            assign(g, list(
              `_id` = g,
              datasets = as.list(character()),
              degs = as.list(character())
            ), envir = rows_map)
          }
          row <- get(g, envir = rows_map, inherits = FALSE)
          # append dataset id if not present
          if (!(d %in% row$datasets)) row$datasets <- append(row$datasets, d)
          # append deg id if not present
          if (!(deg_id %in% row$degs)) row$degs <- append(row$degs, deg_id)
          assign(g, row, envir = rows_map)
          NULL
        })
        NULL
      })
    }

    # Read genes to return only counts
    genes <- tryCatch(
      jsonlite::fromJSON(file.path(ds_dir, ds_index[[d]], ds_genes_file)),
      error = function(e) {
        # fallback to computing when genes.json missing
        tryCatch(
          write_genes(get_ds_data(ds_index[[d]])),
          error = function(e) list()
        )
      }
    )
    genes_count[[d]] <<- length(genes)
    NULL
  })

  # Convert rows map to list and sort by degs count desc then id
  rows <- as.list(mget(ls(envir = rows_map), envir = rows_map))
  rows <- rows[order(sapply(rows, function(x) length(x$degs)), decreasing = TRUE)]
  # tie-breaker: _id ascending
  if (length(rows) > 1) {
    # stable sort by name where degs count equal
    deg_counts <- sapply(rows, function(x) length(x$degs))
    same <- duplicated(deg_counts) | duplicated(deg_counts, fromLast = TRUE)
    if (any(same)) {
      rows_same <- rows[same]
      order_same <- order(vapply(rows_same, function(x) x$`_id`, character(1)))
      rows[same] <- rows_same[order_same]
    }
  }

  list(
    rows = rows,
    degsMeta = degs_meta,
    genesCount = genes_count
  )
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
  groupBy = "CellType", splitBy = "Genotype"
) {
  ds_index <- get_ds_index()
  unlist(lapply(ds, function(d) {
    try({
      if (!d %in% names(ds_data)) ds_data[[d]] <<- get_ds_data(ds_index[[d]])
    })
    unlist(lapply(gene, function(g) {
      lapply(pt, function(p) {
        render_plot(ds_index[[d]], d, g, p, groupBy, splitBy)
      })
    }))
  }))
}
