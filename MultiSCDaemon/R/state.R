.ds_data <- new.env(parent = emptyenv())

ds_loaded <- function() {
  n <- ls(.ds_data)
  if (is.null(n)) {
    list(datasets = list())
  } else {
    list(datasets = n)
  }
}

get_ds_data <- function(ds_path) {
  readRDS(file.path(ds_path, ds_data_file))
}

loaded <- function() {
  ds_loaded()
}

unload <- function(ds) {
  if (missing(ds)) ds <- ls(.ds_data)
  rm(list = ds, envir = .ds_data)
  ds_loaded()
}

load <- function(ds) {
  ds_index <- get_ds_index()
  lapply(ds, function(d) {
    try({
      if (!d %in% ls(.ds_data)) {
        .ds_data[[d]] <- get_ds_data(ds_index[[d]])
      }
    })
  })
  ds_loaded()
}
