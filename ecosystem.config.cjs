const multiscDaemons = {
  name: "multisc-daemon",
  script: "Rscript",
  args: "plumber.R",
  cwd: "./multisc-daemon",
  watch: "daemon.R",
  env: {
    "plumber.port": 8000,
    DATASETS_DIR: "../data/datasets",
    PLOT_DIR: "../data/plots",
    PLOT_FILE: "plot.png",
    DATA_FILE: "data.rds",
    GENOTYPE_COLOR_FILE: "genotype_colors.rds",
    CLUSTER_COLORS_FILE: "cluster_colors.rds",
  },
  instances: 1,
};

const multiscViewer = {
  name: "multisc-viewer",
  script: "build/index.js",
  cwd: "./multisc-viewer",
  watch: false,
  env: {
    NODE_ENV: "production",
    DATASETS_DIR: "../data/datasets",
    DATASETS_META: "meta.json",
    PUBLICATIONS_DIR: "../data/publications",
    PUBLICATIONS_META: "meta.json",
    PLOT_DIR: "../data/plots",
    DAEMON_PORTS: [8000],
  },
};

module.exports = {
  apps: [multiscDaemons, multiscViewer],
};
