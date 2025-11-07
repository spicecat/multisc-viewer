#' @importFrom stats setNames

ds_index_file <- "index.json"
ds_meta_file <- "metadata.json"
ds_data_file <- "data.rds"

pub_index_file <- "index.json"
pub_meta_file <- "metadata.json"

.config <- new.env(parent = emptyenv())

initialize_config <- function() {
  .config$ds_dir <- getOption("MultiSCDaemon.DATASETS_DIR")
  .config$pub_dir <- getOption("MultiSCDaemon.PUBLICATIONS_DIR")
  .config$plots_dir <- getOption("MultiSCDaemon.PLOTS_DIR")
  dir.create(.config$plots_dir, FALSE, TRUE)
}

get_ds_index <- function() {
  lapply(jsonlite::fromJSON(file.path(
    .config$ds_dir, ds_index_file
  )), function(path) {
    normalizePath(
      if (R.utils::isAbsolutePath(path)) {
        path
      } else {
        file.path(.config$ds_dir, path)
      }
    )
  })
}

get_pub_index <- function() {
  lapply(jsonlite::fromJSON(file.path(
    .config$pub_dir, pub_index_file
  )), function(path) {
    normalizePath(
      if (R.utils::isAbsolutePath(path)) {
        path
      } else {
        file.path(.config$pub_dir, path)
      }
    )
  })
}

datasets <- function(ds) {
  ds_index <- get_ds_index()
  if (missing(ds)) ds <- names(ds_index)
  ds_genes <- genes(ds)
  ds_degs <- degs(ds)
  metadata <- setNames(lapply(ds, function(d) {
    ds_deg <- ds_degs[[d]]
    lapply(names(ds_deg), function(deg) {
      ds_deg[[deg]]$gene <<- length(ds_deg[[deg]]$gene)
    })
    tryCatch(
      {
        path <- file.path(ds_index[[d]], ds_meta_file)
        meta <- jsonlite::fromJSON(path)
        meta$gene <- length(ds_genes[[d]])
        meta$deg <- ds_deg
        meta$size <- file.info(path)$size
        meta
      },
      error = function(e) NULL
    )
  }), ds)
  metadata[!vapply(metadata, is.null, logical(1))]
}

publications <- function(pub) {
  pub_index <- get_pub_index()
  if (missing(pub)) pub <- names(pub_index)
  metadata <- setNames(lapply(pub, function(p) {
    tryCatch(
      jsonlite::fromJSON(file.path(pub_index[[p]], pub_meta_file)),
      error = function(e) NULL
    )
  }), pub)
  metadata[!vapply(metadata, is.null, logical(1))]
}
