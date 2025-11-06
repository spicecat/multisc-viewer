#' Runs the MultiSC-Daemon API
#'
#' Launches the Plumber API for the MultiSC-Daemon API.
#'
#' @param port Port to run the API on. Default: 8000
#' @param ... Options passed to \code{plumber::plumb()$run()}
#' @return A running Plumber API
#' @export
msc_plumb <- function(
  port = plumber::get_option_or_env("plumber.port", 8000), ...
) {
  options(
    MultiSCDaemon.DATASETS_DIR =
      normalizePath(
        Sys.getenv("DATASETS_DIR", "../data/datasets"),
        mustWork = TRUE
      ),
    MultiSCDaemon.PUBLICATIONS_DIR =
      normalizePath(
        Sys.getenv("PUBLICATIONS_DIR", "../data/publications"),
        mustWork = TRUE
      ),
    MultiSCDaemon.PLOTS_DIR =
      normalizePath(Sys.getenv("PLOTS_DIR", "../data/plots"),
        mustWork = FALSE
      )
  )
  pr <- plumber::plumb(dir = system.file(
    "plumber", "daemon",
    package = "MultiSCDaemon"
  ))
  pr$run(port = port, ...)
}
