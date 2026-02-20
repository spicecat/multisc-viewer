<!-- README.md is generated from README.Rmd. Please edit that file -->

# MultiSCDaemon

<!-- badges: start -->

<!-- badges: end -->

Backend service to serve data and plots for
[MultiSC-Viewer](https://github.com/spicecat/multisc-viewer/tree/main).

## Installation

Install the development version from the root of the repository (after
cloning) or directly via `remotes`/`pak`.

```r
# Using pak (recommended)
install.packages("pak")
pak::pkg_install("spicecat/multisc-viewer@master", dependencies = TRUE, upgrade = FALSE)

# Or using remotes
install.packages("remotes")
remotes::install_github("spicecat/multisc-viewer", ref = "master", subdir = "MultiSCDaemon")

# From a local clone
# git clone https://github.com/spicecat/multisc-viewer.git
install.packages("pak")
pak::pkg_install("./multisc-viewer/MultiSCDaemon")
```

This package depends on `plumber2`, `Seurat` and other tidyverse/base
packages; `pak` will resolve these automatically.

## Running the daemon

### From R

```r
# Launch on default port (8080)
MultiSCDaemon::msc_plumb()

# Custom port
MultiSCDaemon::msc_plumb(port = 8123)

# Pass additional plumber2::api_run() args (e.g. host)
MultiSCDaemon::msc_plumb(port = 8080, host = "0.0.0.0")
```

### Environment variables

`msc_plumb()` reads the data root from an environment variable before
starting and sets internal `options()` used by API handlers.

| Variable                  | Default   | Description                                                     |
| ------------------------- | --------- | --------------------------------------------------------------- |
| `PLUMBER_PORT`            | 8080      | HTTP port (alternative to `port` argument)                      |
| `MULTISC_VIEWER_DATA_DIR` | `../data` | Data root containing `datasets/`, `publications/`, and `plots/` |

Set them prior to calling `msc_plumb()`:

```bash
export MULTISC_VIEWER_DATA_DIR="/path/to/data"
export PLUMBER_PORT=8123
R -e 'MultiSCDaemon::msc_plumb()'
```

### Docker

A Dockerfile is provided (see `MultiSCViewerR/Dockerfile`). Build and
run:

```bash
docker build -t MultiSCViewerR MultiSCViewerR
docker run -d -p 8080:8080 MultiSCViewerR
```

Override environment + mount data:

```bash
docker run -d \
  -p 8123:8080 \
  -e MULTISC_VIEWER_DATA_DIR=/data \
  -v $(pwd)/data:/data \
  MultiSCViewerR
```

### Access

Once running: `http://localhost:8080/` serves the API root. Interactive
docs (if enabled) usually at `http://localhost:8080/__docs__/`.

## Function reference

`msc_plumb(port = plumber2::get_opts("port", 8080), ...)`

Arguments:

- `port` – numeric port to bind (defaults to env `PLUMBER_PORT` or 8080)
- `...` – forwarded to `plumber2::api_run()` (e.g. `host`,
  `swagger = TRUE`)

Side-effects: sets options

- `MultiSCViewer.DATA_DIR`

This is derived from `MULTISC_VIEWER_DATA_DIR` (with defaults) and is a
normalized path.

Returns: the running Plumber2 API (invisibly) after calling `api_run()`.

## Example request workflow

```r
library(MultiSCDaemon)
msc_plumb(port = 8080)
# Then in another session or using httr/curl you could query endpoints, e.g.:
## Not run:
# httr::GET("http://localhost:8080/datasets")
```

## Development

Render this README after edits:

```r
devtools::build_readme()
```

Run tests / checks (if added):

```r
devtools::check()
```

## Related

- MultiSC Viewer frontend: `multisc-viewer`
- Data directory layout: `data/`
- Daemon (raw R scripts / Docker): `MultiSCViewerR/`
