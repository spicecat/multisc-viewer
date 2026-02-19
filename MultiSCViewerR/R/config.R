#' @keywords internal
.env <- new.env(parent = emptyenv())

#' Initialize Configuration
#'
#' Sets up package-level configuration including data directory,
#' database file path, and indices for datasets and publications.
#'
#' @keywords internal
init <- function() {
  .env$data_dir <- normalizePath(
    getOption("MultiSCViewer.DATA_DIR", "../data"),
    mustWork = TRUE
  )
  print(paste("Using data directory:", .env$data_dir))

  .env$db_file <- file.path(.env$data_dir, "metadata.duckdb")
  .env$plots_dir <- file.path(.env$data_dir, "plots")
  dir.create(.env$plots_dir, FALSE, TRUE)

  ds_paths <- list.files(file.path(.env$data_dir, "datasets"),
    pattern = "metadata.json$",
    recursive = TRUE,
    full.names = TRUE
  )
  .env$ds_index <- setNames(
    dirname(ds_paths),
    sapply(ds_paths, function(f) jsonlite::fromJSON(f)$id)
  )

  pub_paths <- list.files(file.path(.env$data_dir, "publications"),
    pattern = "metadata.json$",
    recursive = TRUE,
    full.names = TRUE
  )
  .env$pub_index <- setNames(
    dirname(pub_paths),
    sapply(pub_paths, function(f) jsonlite::fromJSON(f)$id)
  )
}
