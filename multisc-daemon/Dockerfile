FROM rocker/r-base:latest AS builder

LABEL org.opencontainers.image.title="multisc-daemon" \
	org.opencontainers.image.description="Plotting daemon for Multisc-Viewer"

RUN apt-get update -qq \
	&& apt-get install -y --no-install-recommends \
	libglpk40 libhdf5-dev libsodium-dev \
	libcurl4-openssl-dev libssl-dev libxml2-dev \
	libgit2-dev \
	&& rm -rf /var/lib/apt/lists/*

RUN install2.r --error --skipinstalled pak

RUN R -q -e 'pak::pkg_install(c("plumber", "Seurat", "rjson"))'
RUN R -q -e 'pak::pkg_install(c("readr", "hexbin", "glmGamPoi", "immunogenomics/presto", "url::https://bnprks.r-universe.dev/src/contrib/BPCells_0.3.1.tar.gz"))'

RUN rm -rf /tmp/downloaded_packages/ /tmp/*.rds

COPY . /api
WORKDIR /api
ENV plumber.port=8000 DATASETS_DIR=/data/datasets PLOT_DIR=/data/plots
EXPOSE 8000
CMD ["Rscript", "plumber.R"]
