#!/usr/bin/env Rscript

library(plumber)
pr("daemon.R") %>%
  pr_run(host = "0.0.0.0", port = get_option_or_env("plumber.port", 8000))
