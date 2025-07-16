const R_DAEMON_COUNT = 3;
const R_PORT_START = 8001;

const rApps = Array.from({ length: R_DAEMON_COUNT }, (_, i) => {
  const port = R_PORT_START + i;
  return {
    name: `multisc-daemon-${port}`,
    script: "Rscript",
    args: `-e "plumber::pr_run(plumber::pr('daemon.r'), port = ${port})"`,
    cwd: "./",
    watch: ["daemon.r"],
    autorestart: true,
    restart_delay: 3000,
  };
});

module.exports = {
  apps: [
    {
      name: "multisc-server",
      script: "npm",
      args: "run dev:server",
      cwd: "./",
      watch: false,
      autorestart: false,
    },
    ...rApps,
  ],
};
