FROM rocker/r2u:latest

LABEL org.opencontainers.image.title="multisc-daemon" \
    org.opencontainers.image.description="Plotting daemon for Multisc-Viewer (r2u fast-build variant)"

ARG DEBIAN_FRONTEND=noninteractive
RUN apt-get update -qq \
    && apt-get install -y --no-install-recommends \
    libglpk-dev libsodium-dev libxml2-dev \
    pandoc python3 \
    r-cran-plumber r-cran-seurat \
    && rm -rf /var/lib/apt/lists/*

RUN rm -rf /tmp/downloaded_packages/ /tmp/*.rds

COPY . /api
WORKDIR /api
ENV PLUMBER_PORT=8000 DATASETS_DIR=/data/datasets PLOT_DIR=/data/plots
EXPOSE 8000
CMD ["Rscript", "plumber.R"]
