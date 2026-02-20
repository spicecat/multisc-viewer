library(MultiSCViewerR)

MultiSCViewerR:::init()

#* MultiSCViewerR API
#*
#* Backend service for
#* [MultiSC-Viewer](https://github.com/spicecat/multisc-viewer/tree/main)
#* to serve data and plots.
#*
#* @license Apache 2.0 https://www.apache.org/licenses/LICENSE-2.0.html
#* @version 1.0.0
#* @tag state Load and unload datasets on daemon
#* @tag data Fetch metadata for datasets and publications
#* @tag plots Generate plots for datasets
"_API"

#* @plumber
function(api) {
  api |>
    api_on("end", MultiSCViewerR:::disconnect_db)
}

#* Check daemon health
#* @get /
#* @get /health
#* @serializer unboxedJSON
#* @response 200{status:string} OK
#* @tag state
function() list(status = "ok")

#* Get loaded datasets
#* @get /loaded
#* @serializer json
#* @response 200:[string] Loaded dataset ids
#* @tag state
MultiSCViewerR:::loaded

#* Unload datasets
#* @post /unload
#* @serializer json
#* @query ds:[string] Dataset ids. Unset to unload all datasets.
#* @response 200:[string] Loaded dataset ids
#* @tag state
MultiSCViewerR:::unload_ds

#* Load datasets
#* @post /load
#* @serializer json
#* @query ds:[string]* Dataset ids
#* @response 200:[string] Loaded dataset ids
#* @tag state
MultiSCViewerR:::load_ds

#* Get datasets index
#* @get /datasets-index
#* @serializer json
#* @response 200:[{id:string}] List of datasets ids and sizes
#* @tag data
MultiSCViewerR:::get_datasets_index

#* Get datasets metadata
#* @get /datasets
#* @serializer json
#* @query ds:[string] Dataset ids. Unset to get all datasets.
#* @response 200:[{id:string, deg:[{id:string, name:string}]}]
#*   List of datasets metadata
#* @tag data
MultiSCViewerR:::get_datasets

#* Get publications metadata
#* @get /publications
#* @serializer json
#* @query pub:[string] Publication ids. Unset to get all publications.
#* @response 200:[{id:string, datasets:[string]}] List of publications metadata
#* @tag data
MultiSCViewerR:::get_publications

#* Get genes for datasets
#* @get /genes
#* @serializer json
#* @query ds:[string]* Dataset ids
#* @response 200:[{id:string, gene:[string]}] List of datasets genes
#* @tag data
MultiSCViewerR:::get_genes

#* Get differentially expressed genes for datasets
#* @get /degs
#* @serializer json
#* @query ds:[string]* Dataset ids
#* @response 200:[{id:string, deg:[{id:string, name:string}]}]
#*   List of datasets differentially expressed genes
#* @tag data
MultiSCViewerR:::get_degs

#* Query differentially expressed genes across datasets
#* @get /genes-rows
#* @query ds:[string]* Dataset ids
#* @query gene:[string] Genes query filter
#* @query limit:integer Results limit
#* @query offset:integer(0) Results offset
#* @response 200:[{id:string, gene:[string], deg:[string]}] List of genes rows
#* @tag data
MultiSCViewerR:::get_genes_rows

#* Generate plots for datasets
#* @post /plots
#* @query ds:[string]* Dataset ids
#* @query gene:[string]* Genes
#* @query pt:[enum|umap,vln,feat|]* Plots to generate
#* @query groupBy:string Grouping variable
#* @query splitBy:string Splitting variable
#* @response 200:[string] Generated plot ids
#* @tag plots
MultiSCViewerR:::plots
