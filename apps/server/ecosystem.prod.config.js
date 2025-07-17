const DAEMON_PORTS = [
  8001, 8002, 8003, 8004, 8005, 8006, 8007, 8008, 8009, 8010,
];

const multiscServer = {
  name: "multisc-server",
  script: "node",
  args: "dist/main.js",
  watch: false,
  env: {
    NODE_ENV: "production",
    PORT: 5000,
    DAEMON_PORTS,
  },
};

const multiscDaemons = DAEMON_PORTS.map((port) => ({
  name: `multisc-daemon-${port}`,
  script: "Rscript",
  args: `-e "plumber::pr_run(plumber::pr('daemon.r'), port = ${port})"`,
  watch: false,
}));

export const apps = [multiscServer, ...multiscDaemons];
