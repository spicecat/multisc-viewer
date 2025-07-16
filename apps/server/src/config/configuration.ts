export default () => ({
  server: {
    port: process.env.PORT,
  },

  daemon: {
    maxParallel: parseInt(process.env.DAEMON_MAX_PARALLEL ?? "3", 10),
    ttl: 30 * 60 * 1000, // 30 minutes
    server: "http://localhost",
    ports: [8001, 8002, 8003, 8004, 8005, 8006, 8007, 8008, 8009, 8010],
  },

  // Plot configuration
  plot: {
    ttl: 30 * 60 * 1000, // 30 minutes
  },

  // Publications
  publications: {
    basePath: "publications",
    metaFile: "meta.json",
  },

  // Datasets
  datasets: {
    basePath: "datasets",
    metaFile: "meta.json",
    requiredFiles: ["data.rds", "genes.json"],
    genesPath: "genes.json",
  },
});
