## Description

MultiSC-Viewer is a website to compare and visualize gene expression in multiple single cell/nucleus datasets across different brain regions, disease conditions, and species.
Compare Alzheimer disease (AD), Parkinson disease (PD), and control gene expression and cell clustering.

## Installation

```bash
$ git clone https://git.jasonxu.dev/JasonXu/plot-viewer.git
$ cd plot-viewer
$ pnpm -C multisc-viewer i
```

## Running the app

```bash
$ pnpm -C multisc-viewer build
$ pm2 start ecosystem.config.cjs
```

## Stopping the app

```bash
$ pm2 stop ecosystem.config.cjs
```

View app at [localhost:5173](http://localhost:5173)
