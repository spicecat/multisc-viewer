#' @import data.table
#' @importFrom utils stack

ds_genes_file <- "genes.json"
ds_degs_file <- "degs.json"

#' Write genes to JSON file
#'
#' Writes the gene names from a Seurat object to a JSON file.
#' @param rds A Seurat object. Default reads from "data.rds".
#' @param path The file path to write the JSON file to. Default is "genes.json".
#' @return A character vector of unique gene names.
#' @export
write_genes <- function(rds = readRDS(file = "data.rds"), path = "genes.json") {
  genes <- rownames(rds@assays$RNA@counts)
  jsonlite::write_json(genes, path)
  unique(genes)
}

genes <- function(ds) {
  ds_index <- get_ds_index()
  setNames(lapply(ds, function(d) {
    tryCatch(
      jsonlite::fromJSON(file.path(ds_index[[d]], ds_genes_file)),
      error = function(e) {
        tryCatch(
          write_genes(get_ds_data(ds_index[[d]])),
          error = function(e) list()
        )
      }
    )
  }), ds)
}

degs <- function(ds) {
  ds_index <- get_ds_index()
  metadata <- setNames(lapply(ds, function(d) {
    tryCatch(
      jsonlite::fromJSON(file.path(ds_index[[d]], ds_degs_file)),
      error = function(e) NULL
    )
  }), ds)
  metadata[!vapply(metadata, is.null, logical(1))]
}

genes_rows <- function(ds) {
  ds_genes <- genes(ds)
  ds_degs <- degs(ds)

  if (length(genes) > 0) {
    genes_dt <- as.data.table(stack(ds_genes))
    setnames(genes_dt, c("values", "ind"), c("gene", "datasets"))
    genes_dt <- genes_dt[, .(datasets = list(datasets)), by = "gene"]
  } else {
    return(list())
  }

  if (length(ds_degs) > 0) {
    all_degs <- lapply(names(ds_degs), function(ds) {
      lapply(names(ds_degs[[ds]]), function(deg) {
        data.table(
          gene = ds_degs[[ds]][[deg]]$gene,
          degs = deg
        )
      })
    })
    degs_dt <- rbindlist(unlist(all_degs, recursive = FALSE))
    degs_dt <- degs_dt[, .(degs = list(degs)), by = "gene"]
  } else {
    degs_dt <- data.table(gene = character(0), degs = list())
  }

  combined_dt <- merge(genes_dt, degs_dt, by = "gene", all = TRUE)
  combined_dt[, "len_degs" := lengths(degs)]
  combined_dt[, "len_datasets" := lengths(datasets)]
  setorder(combined_dt, -"len_degs", -"len_datasets", "gene")
  genes_data <- apply(combined_dt, 1, function(row) {
    list(
      `_id` = jsonlite::unbox(row$gene),
      datasets = row$datasets,
      degs = ifelse(row$degs == "", list(), row$degs)
    )
  })

  unname(genes_data)
}
