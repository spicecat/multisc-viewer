#!/usr/bin/env Rscript

plumber::pr("daemon.R") |>
  plumber::pr_run(
    host = "0.0.0.0",
    port = plumber::get_option_or_env("plumber.port", 8000)
  )
