port <- as.integer(Sys.getenv("PLUMBER_PORT", "8080"))
MultiSCViewerR::init_db()
MultiSCViewerR::msc_plumb(port = port)
