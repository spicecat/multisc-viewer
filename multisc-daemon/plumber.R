#!/usr/bin/env Rscript

pr <- plumber::plumb("daemon.R")
pr <- plumber::pr_set_api_spec(pr, "openapi.json")
port <- Sys.getenv("PLUMBER_PORT", 8000)
print(sprintf("Starting multisc-daemon on port %s", port))
plumber::pr_run(pr, "0.0.0.0", port)
