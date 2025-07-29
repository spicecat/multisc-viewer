suppressPackageStartupMessages({
  library(plumber)
  library(Seurat)
  library(ggplot2)
  library(future)
  library(future.apply)
})

plan(multisession)

# --- Cache Directory ---
cache_dir <- Sys.getenv("PLOT_CACHE_DIR", unset = "plots")
if (!dir.exists(cache_dir)) {
  dir.create(cache_dir, recursive = TRUE)
}

# --- Helper Functions ---
get_colors <- function(ds, group_by) {
  color_file <- if (group_by == "Genotype") {
    "genotype.colors.rds"
  } else {
    "cluster.colors.rds"
  }
  readRDS(file.path("datasets", ds, color_file), refhook = NULL)
}

save_plot <- function(plot_object, path, ...) {
  png(path, ...)
  print(plot_object)
  dev.off()
}

# --- Global State ---
datasets <- new.env()

#* @apiTitle MultiSC Plotting Daemon

#* Log some information about the request
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
  list(datasets = ls(datasets))
}

#* Load datasets
#* @post /load
#* @param ds Dataset names to load
function(ds, res) {
  results <- lapply(ds, function(current_ds) {
    if (exists(ds, envir = datasets)) {
      return(list(status = "already loaded"))
    }
    tryCatch(
      {
        datasets[[ds]] <- readRDS(
          file = file.path("datasets", ds, "data.rds"), refhook = NULL
        )
        return(list(status = "loaded"))
      },
      error = function(e) {
        res$status <- 500
        list(error = paste("Failed to load dataset:", e$message))
      }
    )
  })
}

#* Unload a dataset
#* @post /unload
#* @param ds Dataset names to unload
function(ds) {
  if (exists(ds, envir = datasets)) {
    rm(list = ds, envir = datasets)
    return(list(status = "unloaded"))
  }
  list(status = "not found")
}

#* Render plots for a dataset
#* @post /render
#* @param ds Dataset names
#* @param gene Gene to plot
#* @param groupBy Grouping variable
#* @param splitBy Splitting variable
function(ds, gene, groupBy, splitBy, res) {
  results <- future_lapply(ds, function(current_ds) {
    if (!exists(current_ds, envir = datasets)) {
      res$status <- 404
      return(list(error = "Dataset not loaded"))
    }
    dataset <- datasets[[current_ds]]

    cache_key <- paste(current_ds, gene, groupBy, splitBy, sep = ":")
    plot_dir <- file.path(cache_dir, cache_key)
    if (!dir.exists(plot_dir)) {
      dir.create(plot_dir, recursive = TRUE)
    }
    colors <- get_colors(current_ds, groupBy)

    tryCatch(
      {
        # UMAP plot
        umap_plot <- DimPlot(
          dataset,
          reduction = "umap", label = FALSE, group.by = groupBy, cols = colors
        )
        save_plot(umap_plot, file.path(plot_dir, "umap.png"),
          width = 5 * 300, height = 4 * 300, res = 300, pointsize = 12
        )

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

        # Feature plot
        feature_plot <- FeaturePlot(dataset,
          features = c(gene),
          min.cutoff = "q5", max.cutoff = "q95"
        )
        save_plot(feature_plot, file.path(plot_dir, "feat.png"),
          width = 4 * 300, height = 4 * 300, res = 300, pointsize = 5
        )

        return(list(key = cache_key))
      },
      error = function(e) {
        res$status <- 500
        list(error = paste("Failed to render plot:", e$message))
      }
    )
  })
}
