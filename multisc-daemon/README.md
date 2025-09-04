# MultiSC-Daemon

Backend service for [MultiSC-Viewer](https://git.jasonxu.dev/JasonXu/plot-viewer) to serve data and plots.

## Run locally

Requires [R](https://www.r-project.org/).

```bash
# Clone the repository
git clone https://git.jasonxu.dev/JasonXu/plot-viewer.git
cd plot-viewer/multisc-daemon

# Install packages
R -e 'install.packages("pak")'
R -e 'pak::pkg_install(c("plumber", "Seurat"))'

# Install optional packages
R -e 'pak::pkg_install(c("url::https://bnprks.r-universe.dev/src/contrib/BPCells_0.3.1.tar.gz", "immunogenomics/presto", "glmGamPoi"))'

# Start the server
Rscript plumber.R
```

View server at <http://localhost:8000/>
View docs at <http://localhost:8000/__docs__/>

## Docker

```bash
docker build -t multisc-daemon .
docker run -d -p 8000:8000 multisc-daemon
```

## Datasets

## Publications
