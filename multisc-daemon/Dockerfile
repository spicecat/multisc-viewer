FROM rocker/r-base:latest AS builder

LABEL org.opencontainers.image.title="multisc-daemon" \
	org.opencontainers.image.description="Plotting daemon for Multisc-Viewer"

RUN apt-get update -qq \
	&& apt-get install -y --no-install-recommends \
	libhdf5-dev libglpk-dev libsodium-dev libx11-dev libxml2-dev \
	pandoc python3 \
	&& rm -rf /var/lib/apt/lists/*

RUN install2.r --error --skipinstalled pak
RUN R -q -e 'pak::pkg_install(c("plumber", "Seurat"))'
# optional packages
RUN R -q -e 'pak::pkg_install(c("url::https://bnprks.r-universe.dev/src/contrib/BPCells_0.3.1.tar.gz", "immunogenomics/presto", "glmGamPoi"))'

RUN rm -rf /tmp/downloaded_packages/ /tmp/*.rds

COPY . /api
WORKDIR /api
ENV plumber.port=8000 DATASETS_DIR=/data/datasets PLOT_DIR=/data/plots
EXPOSE 8000
CMD ["Rscript", "plumber.R"]
