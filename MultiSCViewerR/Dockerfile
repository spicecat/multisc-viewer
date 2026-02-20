FROM rocker/r2u:24.04
LABEL org.opencontainers.image.title="multisc-viewer-daemon" \
    org.opencontainers.image.description="Backend daemon for MultiSC-Viewer"

WORKDIR /app

ENV MULTISC_VIEWER_DATA_DIR=/data \
    PLUMBER_PORT=8080

COPY renv.lock ./
RUN Rscript -e "install.packages('renv'); renv::restore()"

COPY . /app/pkg
RUN Rscript -e "install.packages('remotes'); remotes::install_local('/app/pkg', dependencies = FALSE)"


EXPOSE 8080

CMD ["R", "--vanilla", "-f", "inst/run_daemon.R"]
