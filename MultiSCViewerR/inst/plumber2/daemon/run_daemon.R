port <- as.integer(Sys.getenv("PLUMBER_PORT", "8080"))
MultiSCViewerR::msc_plumb(port = port)
