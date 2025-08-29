#!/usr/bin/env Rscript

library(plumber)
pr("daemon.R") %>%
  pr_run()
