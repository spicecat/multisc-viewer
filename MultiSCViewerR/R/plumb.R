#' Launch the Plumber2 API
#'
#' Starts the Plumber2 API server for the MultiSC-Viewer daemon.
#'
#' @param port Port to run the API on. Default: 8080
#' @param ... Options passed to \code{plumber2::api_run()}
#' @return A running Plumber API
#' @export
msc_plumb <- function(
  port = plumber2::get_opts("port", 8080L), ...
) {
  options(
    MultiSCViewer.DATA_DIR = Sys.getenv("MULTISC_VIEWER_DATA_DIR", "../data")
  )
  api <- plumber2::api_package("MultiSCViewerR", "daemon")
  api |> plumber2::api_run(port = port, ...)
  api
}
