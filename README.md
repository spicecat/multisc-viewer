# MultiSC-Viewer

[MultiSC-Viewer](https://github.com/spicecat/multisc-viewer/tree/main/multisc-viewer) is a web application for visualizing and comparing gene expression in multiple single cell/nucleus datasets across different brain regions, disease conditions, and species.

## Features

- **Multiple datasets visualization and comparison**: currently the only available tool that supports a side-by-side comparative view of more than two datasets simultaneously
- **Interactive plots**: interactive plots controls for real-time exploration of datasets
- **Plot streaming**: renders plots as they are generated, allowing users to view results without waiting for all plots to complete
- **Load balancing and scalability**: manages datasets across multiple plotting daemons for parallelized plot generation.

## Run Locally

Requires [R](https://www.r-project.org/), [Seurat](https://satijalab.org/seurat/), [plumber2](https://plumber2.posit.co/), [Node.js](https://nodejs.org/), [PM2](https://pm2.keymetrics.io/).

```bash
# Clone the repository
git clone --single-branch https://github.com/spicecat/multisc-viewer.git
cd multisc-viewer

# Install system dependencies (Debian/Ubuntu)
sudo apt update
sudo apt-get install libglpk-dev libsodium-dev libx11-dev libxml2-dev pandoc python3

# Install packages and dependencies
R -e 'install.packages("pak")'
R -e 'pak::pkg_install(c("plumber2", "Seurat"))'

npm --prefix multisc-viewer install
npm --prefix multisc-viewer run build
```

```bash
# Start the app
pm2 start ecosystem.config.js

# Scale the number of plotting daemons
pm2 scale MultiSCViewerR 4
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

- MultiSC-Viewer: [https://github.com/spicecat/multisc-viewer/tree/main/multisc-viewer](https://github.com/spicecat/multisc-viewer/tree/main/README.md)
- Data: [https://github.com/spicecat/multisc-viewer/tree/main/data](https://github.com/spicecat/multisc-viewer/tree/main/data/README.md)
- MultiSC-ViewerR Daemon: [https://github.com/spicecat/multisc-viewer/tree/main/MultiSCViewerR](https://github.com/spicecat/multisc-viewer/tree/main/MultiSCViewerR/README.md)
- MultiSC-Viewer Web App: [https://github.com/spicecat/multisc-viewer/tree/main/multisc-viewer](https://github.com/spicecat/multisc-viewer/tree/main/multisc-viewer/README.md)
