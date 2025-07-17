const multiscServer = {
  name: "multisc-server",
  script: "nest start",
  max_restarts: 3,
  restart_delay: 1000,
  env_production: {
    watch: false,
    NODE_ENV: "production",
    PORT: 5000,
    DAEMON_PORTS: [8001, 8002, 8003], //, 8004, 8005, 8006, 8007, 8008, 8009, 8010],
  },
  env_development: {
    watch: ["src"],
    NODE_ENV: "development",
    PORT: 5000,
    DAEMON_PORTS: [8001, 8002, 8003],
  },
};

const multiscDaemon = (port) => ({
  name: `multisc-daemon-${port}`,
  script: "Rscript",
  args: `--vanilla --quiet -e "plumber::pr_run(plumber::pr('daemon.r'), port = ${port})"`,
  max_restarts: 3,
  restart_delay: 1000,
  watch: false,
  env_production: {
    args: `--version`,
    autorestart: false,
  },
  env_development: {
    args: `--version`,
    autorestart: false,
  },
});

const multiscDaemons = new Map(
  multiscServer.env_production.DAEMON_PORTS.map((port) => [
    port,
    {
      ...multiscDaemon(port),
      env_production: {},
    },
  ])
);

multiscServer.env_development.DAEMON_PORTS.forEach((port) => {
  if (multiscDaemons.has(port))
    multiscDaemons.get(port).env_development = {
      watch: ["daemon.r"],
    };
  else
    multiscDaemons.set(port, {
      ...multiscDaemon(port),
      env_development: {
        watch: ["daemon.r"],
      },
    });
});

module.exports = {
  apps: [multiscServer, ...multiscDaemons.values()],
};
