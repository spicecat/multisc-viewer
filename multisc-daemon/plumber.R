#!/usr/bin/env Rscript

plumber::pr("daemon.R") |>
  plumber::pr_run(
    host = "0.0.0.0",
    port = Sys.getenv("PLUMBER_PORT", 8000)
  )
