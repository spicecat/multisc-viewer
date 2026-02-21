author_schema <- "STRUCT(name VARCHAR)[]"
ontology_term_schema <- "STRUCT(
  displayName VARCHAR,
  identifier VARCHAR,
  name VARCHAR,
  url VARCHAR
)[]"
dataset_schema <- c(
  id = "VARCHAR PRIMARY KEY",
  name = "VARCHAR",
  displayName = "VARCHAR",
  description = "VARCHAR",
  identifier = "VARCHAR",
  url = "VARCHAR",
  publicationDate = "DATE",
  author = author_schema,
  citation = "STRUCT(identifier VARCHAR, name VARCHAR, url VARCHAR)[]",
  species = ontology_term_schema,
  healthCondition = ontology_term_schema,
  tissue = ontology_term_schema,
  cellType = ontology_term_schema,
  defaultGene = "VARCHAR",
  gene = "VARCHAR[]",
  deg = "STRUCT(id VARCHAR, name VARCHAR, gene VARCHAR[])[]",
  size = "UBIGINT"
)
publication_schema <- c(
  id = "VARCHAR PRIMARY KEY",
  name = "VARCHAR",
  description = "VARCHAR",
  identifier = "VARCHAR",
  url = "VARCHAR",
  journalName = "VARCHAR",
  publicationDate = "DATE",
  author = author_schema,
  datasets = "VARCHAR[]"
)

duckdb_columns <- function(schema) {
  clean_types <- gsub(" PRIMARY KEY", "", schema)
  parts <- paste0("'", names(clean_types), "': '", clean_types, "'")
  paste0("{", paste(parts, collapse = ", "), "}")
}

convert_rds_to_qs2 <- function(dir) {
  qs2_path <- file.path(dir, "data.qs2")
  if (file.exists(qs2_path)) {
    return(qs2_path)
  }

  rds_path <- file.path(dir, "data.rds")
  if (file.exists(rds_path)) {
    print(paste("Converting RDS to QS2:", dir))
    data <- readRDS(rds_path)
    qs2::qs_save(data, qs2_path)
    qs2_path
  }
}

insert_datasets <- function(con) {
  ds_cols <- duckdb_columns(dataset_schema)
  for (dir in .env$ds_index) {
    print(paste("Adding:", dir))
    meta <- jsonlite::fromJSON(file.path(dir, "metadata.json"))

    # Convert RDS to QS2 if necessary
    qs_path <- convert_rds_to_qs2(dir)

    # Add gene, deg, size to metadata
    data <- qs2::qs_read(qs_path)
    meta$gene <- rownames(data[["RNA"]])
    if (file.exists(file.path(dir, "degs.json"))) {
      meta$deg <- jsonlite::fromJSON(file.path(dir, "degs.json"))
    }
    meta$size <- as.numeric(object.size(data))

    ds_fields <- intersect(names(dataset_schema), names(meta))
    temp_meta_path <- tempfile(pattern = "temp_meta_", fileext = ".json")
    jsonlite::write_json(
      meta[ds_fields],
      temp_meta_path,
      auto_unbox = TRUE,
      pretty = TRUE
    )
    DBI::dbExecute(
      con,
      sprintf(
        "INSERT INTO datasets SELECT * FROM read_json(?, columns = %s)",
        ds_cols
      ),
      params = list(temp_meta_path)
    )
  }
}

insert_publications <- function(con) {
  pub_cols <- duckdb_columns(publication_schema)
  for (dir in .env$pub_index) {
    print(paste("Adding publication:", dir))

    meta <- jsonlite::fromJSON(file.path(dir, "metadata.json"))
    pub_fields <- intersect(names(publication_schema), names(meta))
    temp_meta_path <- tempfile(pattern = "temp_pub_meta_", fileext = ".json")
    jsonlite::write_json(
      meta[pub_fields],
      temp_meta_path,
      auto_unbox = TRUE,
      pretty = TRUE
    )
    DBI::dbExecute(
      con,
      sprintf(
        "INSERT INTO publications SELECT * FROM read_json(?, columns = %s)",
        pub_cols
      ),
      params = list(temp_meta_path)
    )
  }
}

#' Initialize the Database
#'
#' Creates and populates a DuckDB database from metadata and data
#' files in the configured data directory. This function reads dataset
#' and publication metadata, integrates gene expression data, and
#' loads everything into DuckDB tables.
#'
#' @details
#' The function performs the following steps:
#' 1. Initializes configuration and paths via `init()`
#' 2. Creates a DuckDB connection
#' 3. Resets existing tables
#' 4. Creates `datasets` and `publications` tables with appropriate schemas
#' 5. Populates tables with metadata from JSON files
#'    and gene data from QS2 files
#'
#' @return Invisibly returns `NULL`. Creates or updates the DuckDB
#'   database file.
#'
#' @importFrom stats setNames
#' @importFrom utils object.size
#'
#' @export
#'
#' @examples
#' \dontrun{
#' init_db()
#' }
init_db <- function() {
  init()

  con <- DBI::dbConnect(
    duckdb::duckdb(),
    dbdir = .env$db_file, read_only = FALSE
  )
  DBI::dbExecute(con, "INSTALL json; LOAD json;")

  # Reset existing tables
  for (t in DBI::dbListTables(con)) DBI::dbRemoveTable(con, t)
  DBI::dbCreateTable(con, "datasets", dataset_schema)
  DBI::dbCreateTable(con, "publications", publication_schema)

  insert_datasets(con)
  insert_publications(con)

  disconnect_db()
}
