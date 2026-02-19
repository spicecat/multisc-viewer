#' @keywords internal
.ds_data <- new.env(parent = emptyenv())

#' Get Loaded Datasets
#'
#' Returns IDs of currently loaded datasets.
#'
#' @return Character vector of dataset IDs
#' @keywords internal
loaded <- function() {
  ls(.ds_data)
}

#' Unload Datasets
#'
#' Unloads specified datasets. If none specified, unloads all datasets.
#'
#' @param query List containing optional `ds` element with dataset IDs to unload
#' @return Character vector of remaining loaded dataset IDs
#' @keywords internal
unload_ds <- function(query) {
  if (length(query$ds) == 0) query$ds <- ls(.ds_data)
  rm(list = query$ds, envir = .ds_data)
  ls(.ds_data)
}

#' Load Datasets
#'
#' Loads Seurat objects from QS2 files and caches them in memory.
#'
#' @param query List containing `ds` element with dataset IDs to load
#' @return Character vector of currently loaded dataset IDs
#' @keywords internal
load_ds <- function(query) {
  lapply(query$ds, function(d) {
    try({
      if (!d %in% ls(.ds_data)) {
        .ds_data[[d]] <- qs2::qs_read(
          file.path(.env$ds_index[[d]], "data.qs2")
        )
      }
    })
  })
  ls(.ds_data)
}
