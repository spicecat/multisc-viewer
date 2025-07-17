const DAEMON_PORTS = [8001, 8002, 8003];

const multiscServer = {
  name: "multisc-server",
  script: "nest start",
  watch: ["src"],
  env: {
    NODE_ENV: "development",
    PORT: 5000,
    DAEMON_PORTS: DAEMON_PORTS,
  },
};

const multiscDaemons = DAEMON_PORTS.map((port) => ({
  name: `multisc-daemon-${port}`,
  script: "Rscript",
  args: `-e "plumber::pr_run(plumber::pr('daemon.r'), port = ${port})"`,
  watch: ["daemon.r"],
}));

export const apps = [multiscServer, ...multiscDaemons];