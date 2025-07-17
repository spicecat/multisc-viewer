export default () => ({
  server: {
    port: process.env.PORT,
  },

  publications: {
    basePath: "publications",
    metaFile: "meta.json",
  },

  datasets: {
    basePath: "datasets",
    metaFile: "meta.json",
    requiredFiles: ["data.rds", "genes.json"],
    genesPath: "genes.json",
  },

  plot: {
    ttl: 30 * 60 * 1000, // 30 minutes
  },

  daemon: {
    ttl: 30 * 60 * 1000, // 30 minutes
    server: "http://localhost",
    ports: (process.env.DAEMON_PORTS ?? "").split(",").map(Number),
  },
});
