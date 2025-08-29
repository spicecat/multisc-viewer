const multiscDaemon = {
  name: "multisc-daemon",
  script: "Rscript",
  args: "plumber.R",
  cwd: "./multisc-daemon",
  watch: "daemon.R",
  instances: 1,
};

const multiscViewer = {
  name: "multisc-viewer",
  script: "pnpm",
  args: "dev",
  cwd: "./multisc-viewer",
  watch: false,
  env: {
    NODE_ENV: "production",
    DAEMON_PORTS: [8000],
  },
};

module.exports = {
  apps: [multiscDaemon, multiscViewer],
};
