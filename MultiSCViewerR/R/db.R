

#' Disconnect from Database
#'
#' Closes the database connection if it exists and is valid.
#'
#' @keywords internal
disconnect_db <- function() {
  if (!is.null(.env$con) && DBI::dbIsValid(.env$con)) {
    DBI::dbDisconnect(.env$con, shutdown = TRUE)
    .env$con <- NULL
  }
}

connect_db <- function() {
  if (is.null(.env$con) || !DBI::dbIsValid(.env$con)) {
    .env$con <- DBI::dbConnect(
      duckdb::duckdb(),
      dbdir = .env$db_file, read_only = TRUE
    )
  }
  .env$con
}
