# MultiSC-Daemon

Backend service to serve data and plots for [MultiSC-Viewer](https://git.jasonxu.dev/JasonXu/plot-viewer/src/branch/main).

## Run Locally

Requires [R](https://www.r-project.org/), [Seurat](https://satijalab.org/seurat/), [plumber](https://www.rplumber.io/).

```bash
# Clone the repository
git clone --single-branch -b main https://git.jasonxu.dev/JasonXu/plot-viewer.git
cd plot-viewer/multisc-daemon

# Install system dependencies (Debian/Ubuntu)
sudo apt update
sudo apt-get install libglpk-dev libsodium-dev libx11-dev libxml2-dev pandoc python3

# Install packages
R -e 'install.packages("pak")'
R -e 'pak::pkg_install(c("plumber", "Seurat"))'

# Start the server
Rscript plumber.R

R -e 'MultiSCDaemon::msc_plumb()'
```

View server at <http://localhost:8000/>
View docs at <http://localhost:8000/__docs__/>

## Docker

Requires [Docker](https://www.docker.com/).

```bash
docker build -t multisc-daemon .
docker run -d -p 8000:8000 multisc-daemon
```

## Environment Variables

The daemon can be configured via environment variables (defaults in parentheses):

- PLUMBER_PORT (8000) — HTTP port for the server
- DATASETS_DIR (../data/datasets) — path to datasets directory
- PUBLICATIONS_DIR (../data/publications) — path to publications directory
- PLOTS_DIR (../data/plots) — output directory for rendered plots

Examples

```bash
# Local override
export DATASETS_DIR="/path/to/data/datasets"
export PUBLICATIONS_DIR="/path/to/data/publications"
export PLOTS_DIR="/path/to/data/plots"
export PLUMBER_PORT="8123"
Rscript plumber.R
```

```bash
# Docker override
docker run -d \
	-p 8123:8000 \
	-e DATASETS_DIR=/data/datasets \
	-e PUBLICATIONS_DIR=/data/publications \
	-e PLOTS_DIR=/data/plots \
	-v $(pwd)/data:/data \
	multisc-daemon
```

---

## Related

- MultiSC: [https://git.jasonxu.dev/JasonXu/plot-viewer/src/branch/main/multisc-viewer](https://git.jasonxu.dev/JasonXu/plot-viewer/src/branch/main/README.md)
- Add MultiSC-Data: [https://git.jasonxu.dev/JasonXu/plot-viewer/src/branch/main/data](https://git.jasonxu.dev/JasonXu/plot-viewer/src/branch/main/data/README.md)
- Run MultiSC-Daemon: [https://git.jasonxu.dev/JasonXu/plot-viewer/src/branch/main/multisc-daemon](https://git.jasonxu.dev/JasonXu/plot-viewer/src/branch/main/multisc-daemon/README.md)
- Run MultiSC-Viewer: [https://git.jasonxu.dev/JasonXu/plot-viewer/src/branch/main/multisc-viewer](https://git.jasonxu.dev/JasonXu/plot-viewer/src/branch/main/multisc-viewer/README.md)
