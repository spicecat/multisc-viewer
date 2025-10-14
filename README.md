# MultiSC-Viewer

[MultiSC-Viewer](https://git.jasonxu.dev/JasonXu/plot-viewer/src/branch/main/multisc-viewer) is a web application for visualizing and comparing gene expression in multiple single cell/nucleus datasets across different brain regions, disease conditions, and species.

## Features

- **Multiple datasets visualization and comparison**: currently the only available tool that supports a side-by-side comparative view of more than two datasets simultaneously
- **Interactive plots**: interactive plots controls for real-time exploration of datasets
- **Plot streaming**: renders plots as they are generated, allowing users to view results without waiting for all plots to complete
- **Load balancing and scalability**: manages datasets across multiple plotting daemons for parallelized plot generation.

## Run Locally

Requires [R](https://www.r-project.org/), [Seurat](https://satijalab.org/seurat/), [plumber](https://www.rplumber.io/), [Node.js](https://nodejs.org/), [PM2](https://pm2.keymetrics.io/).

```bash
# Clone the repository
git clone --single-branch -b main https://git.jasonxu.dev/JasonXu/plot-viewer.git
cd plot-viewer

# Install system dependencies (Debian/Ubuntu)
sudo apt update
sudo apt-get install libglpk-dev libsodium-dev libx11-dev libxml2-dev pandoc python3

# Install packages and dependencies
R -e 'install.packages("pak")'
R -e 'pak::pkg_install(c("plumber", "Seurat"))'

npm --prefix multisc-viewer install
npm --prefix multisc-viewer run build
```

```bash
# Start the app
pm2 start ecosystem.config.js

# Scale the number of plotting daemons
pm2 scale multisc-daemon 4
```

```bash
# Stop the app
pm2 stop ecosystem.config.js
```

View app at <http://localhost:3000>.

## Docker

```bash
docker compose up -d
```

---

## Related

- MultiSC: [https://git.jasonxu.dev/JasonXu/plot-viewer/src/branch/main/multisc-viewer](https://git.jasonxu.dev/JasonXu/plot-viewer/src/branch/main/README.md)
- Add MultiSC-Data: [https://git.jasonxu.dev/JasonXu/plot-viewer/src/branch/main/data](https://git.jasonxu.dev/JasonXu/plot-viewer/src/branch/main/data/README.md)
- Run MultiSC-Daemon: [https://git.jasonxu.dev/JasonXu/plot-viewer/src/branch/main/multisc-daemon](https://git.jasonxu.dev/JasonXu/plot-viewer/src/branch/main/multisc-daemon/README.md)
- Run MultiSC-Viewer: [https://git.jasonxu.dev/JasonXu/plot-viewer/src/branch/main/multisc-viewer](https://git.jasonxu.dev/JasonXu/plot-viewer/src/branch/main/multisc-viewer/README.md)
