# MultiSC-Viewer

MultiSC-Viewer is a web application to visualize and compare gene expression in multiple single cell/nucleus datasets across different brain regions, disease conditions, and species.
Compare Alzheimer disease (AD), Parkinson disease (PD), and control gene expression and cell clustering.

## Run locally

Requires [Node.js](https://nodejs.org/), [PM2](https://pm2.keymetrics.io/), and [R](https://www.r-project.org/).

```bash
# Clone the repository
git clone https://git.jasonxu.dev/JasonXu/plot-viewer.git
cd plot-viewer

# Install dependencies
npm --prefix multisc-viewer install
npm --prefix multisc-viewer run build
```

```bash
# Start the app
pm2 start ecosystem.json
```

```bash
# Stop the app
pm2 stop ecosystem.json
```

View app at <http://localhost:3000>.

## Docker

```bash
docker compose up -d
```
