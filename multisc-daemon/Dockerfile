FROM rstudio/plumber:v1.3.0
WORKDIR /app
RUN install2.r --error plumber Seurat rjson
ENV DATASETS_DIR=/data/datasets \
	PLOT_CACHE_DIR=/data/plots/cache \
	PORT=8001

COPY daemon.R ./

EXPOSE 8001
CMD ["R","-e","plumber::pr_run(plumber::pr('daemon.R'),port=as.integer(Sys.getenv('PORT','8001')))" ]
