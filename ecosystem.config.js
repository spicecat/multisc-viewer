module.exports = {
  apps: [
    {
      name: "MultiSCViewerR",
      cwd: "./MultiSCViewerR",
      script: "inst/plumber2/run_daemon.R",
      interpreter: "Rscript",
      instances: 1,
      exec_mode: "fork",
      watch: true,
      increment_var: "PLUMBER_PORT",
      autorestart: false,
      env: { PLUMBER_PORT: 8080 },
    },
    {
      name: "multisc-viewer",
      cwd: "./multisc-viewer",
      script: "./build/index.js",
      watch: false,
      env: {
        NODE_ENV: "development",
      },
    },
  ],
};
