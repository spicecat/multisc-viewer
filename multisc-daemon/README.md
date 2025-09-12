# MultiSC-Daemon

Backend service for [MultiSC-Viewer](https://git.jasonxu.dev/JasonXu/plot-viewer) to serve data and plots.

## Run locally

Requires [R](https://www.r-project.org/).

```bash
# Clone the repository
git clone https://git.jasonxu.dev/JasonXu/plot-viewer.git
cd plot-viewer/multisc-daemon

# Install system dependencies (Debian/Ubuntu)
sudo apt update
sudo apt-get install libglpk-dev libsodium-dev libx11-dev libxml2-dev pandoc python3

# Install packages
R -e 'install.packages("pak")'
R -e 'pak::pkg_install(c("plumber", "Seurat"))'

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

## Environment variables

The daemon can be configured via environment variables (defaults in parentheses):

- DATASETS_DIR (../data/datasets) — path to datasets directory
- DATASETS_META (meta.json) — datasets metadata filename
- PUBLICATIONS_DIR (../data/publications) — path to publications directory
- PUBLICATIONS_META (meta.json) — publications metadata filename
- PLOT_DIR (../data/plots) — output directory for rendered plots
- PLOT_FILE (plot.png) — filename used per rendered plot
- DATA_FILE (data.rds) — dataset RDS filename within each dataset directory
- GENOTYPE_COLOR_FILE (genotype.colors.rds) — filename for genotype color map
- CLUSTER_COLOR_FILE (cluster.colors.rds) — filename for cluster color map
- PLUMBER_PORT (8000) — HTTP port for the server

Examples

```bash
# Local override
export DATASETS_DIR="/path/to/data/datasets"
export PLOT_DIR="/path/to/data/plots"
export PLUMBER_PORT="8000"
Rscript plumber.R

# Docker override
docker run -d \
	-p 8000:8000 \
	-e DATASETS_DIR=/data/datasets \
	-e PUBLICATIONS_DIR=/data/publications \
	-e PLOT_DIR=/data/plots \
	-v $(pwd)/data:/data \
	multisc-daemon
```

## Datasets

Add datasets to the `data/datasets` directory and `data/datasets/meta.json` file. Each dataset is a directory that includes:

- `data.rds` — a Seurat object
- `genes.json` — an array of genes
- `cluster.colors.rds` — a vector mapping cluster labels to colors
- `genotype.colors.rds` — a vector mapping genotypes/conditions to colors

Generate `genes.json` by running the `genes.R` script from each dataset directory (where `data.rds` is located):

```bash
cd path/to/My_Dataset
Rscript ../../../multisc-daemon/genes.R
```

Example `cluster.colors.rds`

```r
cluster_colors <- c(
	"K_RORB+" = "#00B4F0",
	"K_FAT2+" = "royalblue3",
	"K_TH_SOX6hi" = "#D19300",
	"K_TH_SOX6lo" = "#F17D50"
	# ... add all clusters present in your dataset
)
saveRDS(cluster_colors, "cluster.colors.rds")
```

Example `genotype.colors.rds`

```r
genotype_colors <- c(
	CTRL = "green",
	DLBD = "orange1",
	PD   = "blue",
	PDD  = "purple2",
	HD   = "midnightblue",
	AD   = "red",
	NHD  = "sienna4"
)
saveRDS(genotype_colors, "genotype.colors.rds")
```

- Use `data/datasets/example_meta.json` as a template for `data/datasets/meta.json`.
- Ensure each entry’s `id` matches the dataset directory name (e.g., `data/datasets/<id>`).
- Required fields: `id`, `title`, `year`, `authors`, `PMID`, `region`, `disease`, `cellType`, `species`.
- Optional field: `defaultGenes`.

## Publications

Add publications to the `data/publications/meta.json` file.

- Use `data/publications/example_meta.json` as a template for `data/publications/meta.json`.
- Required fields: `id`, `title`, `year`, `authors`, `PMID`, `journal`, `abstract`, `datasets`.
  - `abstract` is raw HTML.
  - `datasets` is an array of dataset IDs included in the publication.

Example data layout:

```
plot-viewer
├── data
│   ├── datasets
│   │   ├── My_Dataset
│   │   │   ├── cluster.colors.rds
│   │   │   ├── data.rds
│   │   │   ├── genes.json
│   │   │   └── genotype.colors.rds
│   │   ├── Another_Dataset
│   │   │   ├── cluster.colors.rds
│   │   │   ├── data.rds
│   │   │   ├── genes.json
│   │   │   └── genotype.colors.rds
│   │   ├── meta.json
│   └── publications
│       └── meta.json
├── multisc-daemon
│   ├── daemon.R
│   ├── genes.R
│   ├── plumber.R
│   └── ...
├── multisc-viewer
│   └── ...
└── ...
```
