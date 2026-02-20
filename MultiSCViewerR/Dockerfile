FROM rocker/r2u:24.04
LABEL org.opencontainers.image.title="multisc-viewer-daemon" \
    org.opencontainers.image.description="Backend daemon for MultiSC-Viewer"

WORKDIR /app
COPY . /app

RUN R -e "install.packages('pak', repos = 'https://cloud.r-project.org'); pak::pkg_install('.', dependencies = TRUE, upgrade = FALSE)"

ENV MULTISC_VIEWER_DATA_DIR=/data \
    PLUMBER_PORT=8080
EXPOSE 8080

CMD ["R", "--vanilla", "-f", "inst/plumber2/daemon/run_daemon.R"]
