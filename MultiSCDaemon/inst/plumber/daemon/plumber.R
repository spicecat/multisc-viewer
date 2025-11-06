library(MultiSCDaemon)
library(plumber)

MultiSCDaemon:::initialize_config()

#* @apiTitle MultiSC-Daemon API

#* @plumber
function(pr) {
  pr %>%
    pr_set_api_spec("openapi.json")
}

#* Log request information
#* @filter logger
function(req) {
  cat(
    as.character(Sys.time()), "-", req$REQUEST_METHOD, req$PATH_INFO, "-",
    req$HTTP_USER_AGENT, "@", req$REMOTE_ADDR, "\n"
  )
  plumber::forward()
}

#* Check daemon health.
#* @get /
#* @get /health
#* @serializer unboxedJSON
function() list(status = "ok")

#* Get loaded datasets on daemon.
#* @get /loaded
MultiSCDaemon:::loaded

#* Unload datasets on daemon.
#* @post /unload
#* @param ds:[str] Dataset ids. Unset to unload all datasets.
MultiSCDaemon:::unload

#* Load datasets on daemon.
#* @post /load
#* @param ds:[str]* Dataset ids
MultiSCDaemon:::load

#* Get datasets metadata.
#* @get /datasets
#* @serializer unboxedJSON
#* @param ds:[str] Dataset ids. Unset to get all datasets.
MultiSCDaemon:::datasets

#* Get publications metadata.
#* @get /publications
#* @serializer unboxedJSON
#* @param pub:[str] Publication ids. Unset to get all publications.
MultiSCDaemon:::publications

#* Get genes for datasets.
#* @get /genes
#* @serializer unboxedJSON
#* @param ds:[str]* Dataset ids
MultiSCDaemon:::genes

#* Get differentially expressed genes for datasets.
#* @get /degs
#* @serializer unboxedJSON
#* @param ds:[str]* Dataset ids
MultiSCDaemon:::degs

#* Get datasets for genes and differentially expressed genes.
#* @get /genes-rows
#* @param ds:[str]* Dataset ids
MultiSCDaemon:::genes_rows

#* Generate plots for datasets.
#* @post /plots
#* @param ds:[str]* Dataset ids
#* @param gene:[str]* Genes
#* @param pt:[str] Plots to render
#* @param groupBy:str Grouping variable
#* @param splitBy:str Splitting variable
MultiSCDaemon:::plots
