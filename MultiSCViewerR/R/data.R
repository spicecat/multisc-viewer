#' Get Datasets Index
#'
#' Retrieves dataset IDs and sizes.
#'
#' @return Data frame with dataset IDs and sizes
#' @keywords internal
get_datasets_index <- function() {
  con <- connect_db()
  q <- "SELECT id, size FROM datasets"
  DBI::dbGetQuery(con, q)
}

#' Get Datasets Metadata
#'
#' Retrieves dataset metadata.
#'
#' @param query List with optional `ds` element for dataset IDs
#' @return Data frame with dataset metadata
#' @keywords internal
get_datasets <- function(query) {
  con <- connect_db()
  q <- "SELECT * EXCLUDE (gene, size)
REPLACE (
  list_transform(
    deg,
    x -> struct_pack(id := x.id, name := x.name, genes := len(x.gene))
  ) AS deg
),
len(gene) AS genes
FROM datasets"
  if (length(query$ds) > 0) {
    q <- paste0(
      q, " WHERE id IN (",
      paste(rep("?", length(query$ds)), collapse = ", "),
      ")"
    )
  }
  DBI::dbGetQuery(con, q, params = c(query$ds))
}

#' Get Publications Metadata
#'
#' Retrieves publication metadata.
#'
#' @param query List with optional `pub` element for publication IDs
#' @return Data frame with publication metadata
#' @keywords internal
get_publications <- function(query) {
  con <- connect_db()
  q <- "SELECT * FROM publications"
  if (length(query$pub) > 0) {
    q <- paste0(
      q, " WHERE id IN (",
      paste(rep("?", length(query$pub)), collapse = ", "),
      ")"
    )
  }
  DBI::dbGetQuery(con, q, params = c(query$pub))
}

#' Get Genes for Datasets
#'
#' Retrieves gene lists for datasets.
#'
#' @param query List with `ds` element containing dataset IDs
#' @return Data frame with gene information
#' @keywords internal
get_genes <- function(query) {
  con <- connect_db()
  q <- paste0(
    "SELECT id, gene FROM datasets WHERE id IN (",
    paste(rep("?", length(query$ds)), collapse = ", "),
    ")"
  )
  DBI::dbGetQuery(con, q, params = c(query$ds))
}

#' Get Differentially Expressed Genes for Datasets
#'
#' Retrieves DEG information for datasets.
#'
#' @param query List with `ds` element containing dataset IDs
#' @return Data frame with DEG information
#' @keywords internal
get_degs <- function(query) {
  con <- connect_db()
  q <- paste0(
    "SELECT id, deg FROM datasets WHERE id IN (",
    paste(rep("?", length(query$ds)), collapse = ", "),
    ")"
  )
  DBI::dbGetQuery(con, q, params = c(query$ds))
}


#' Query Differentially Expressed Genes Across Datasets
#'
#' Searches genes across datasets and returns paginated results.
#'
#' @param query List containing `ds` (dataset IDs), `gene` (search pattern),
#'   `limit` (results per page), and `offset` (pagination offset)
#' @return Data frame with gene rows sorted by DEG/gene frequency
#' @keywords internal
get_genes_rows <- function(query) {
  con <- connect_db()
  q <- paste0(
    "WITH ds as (SELECT * FROM datasets WHERE id in (",
    paste(rep("?", length(query$ds)), collapse = ", "),
    ")), gene_ds AS (SELECT unnest(gene) as gene, id as ds_id FROM ds),
gene_deg AS (SELECT unnest(deg.gene) as gene, deg.id as deg_id
  FROM (SELECT unnest(deg) as deg FROM ds)),
genes AS (SELECT gene FROM gene_ds UNION SELECT gene FROM gene_deg),
gene_rows AS (SELECT  g.gene AS id,
  (SELECT list(ds_id) FROM gene_ds WHERE gene = g.gene) as gene,
  (SELECT COALESCE(list(deg_id), []) FROM gene_deg WHERE gene = g.gene) as deg
  FROM genes g)
SELECT * FROM gene_rows WHERE regexp_matches(id, ?)
ORDER BY len(deg) DESC, len(gene) DESC, id ASC"
  )
  if (!is.null(query$limit)) q <- paste(q, "LIMIT ?")
  if (!is.null(query$offset)) q <- paste(q, "OFFSET ?")

  pattern <- paste(query$gene, collapse = "|")
  params <- c(query$ds, pattern, query$limit, query$offset)
  DBI::dbGetQuery(con, q, params = params)
}
