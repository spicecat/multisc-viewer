FROM rocker/r-base:latest AS builder

LABEL org.opencontainers.image.title="multisc-daemon" \
	org.opencontainers.image.description="Plotting daemon for Multisc-Viewer"

RUN apt-get update -qq \
	&& apt-get install -y --no-install-recommends \
	libglpk40 libhdf5-dev libsodium-dev \
	libcurl4-openssl-dev libssl-dev libxml2-dev \
	libgit2-dev \
	&& rm -rf /var/lib/apt/lists/*

RUN install2.r \
	--ncpus -1 --error --skipinstalled \
	--repos 'https://packagemanager.posit.co/cran/latest' \
	plumber Seurat rjson

# RUN install2.r \
# 	--ncpus -1 --error --skipinstalled \
# 	readr hexbin glmGamPoi

# RUN install2.r \
# 	--ncpus -1 --error --skipinstalled \
# 	--repos 'https://satijalab.r-universe.dev' \
# 	--repos 'https://bnprks.r-universe.dev' \
# 	presto BPCells \

RUN rm -rf /tmp/downloaded_packages/ /tmp/*.rds

RUN mkdir -p /rlib && cp -a "$(R RHOME)/library" /rlib/site-library

FROM rocker/r-base:latest

RUN apt-get update -qq \
	&& apt-get install -y --no-install-recommends \
	libglpk40 libhdf5-dev \
	libcurl4-openssl-dev libssl-dev libxml2-dev libgit2-dev \
	&& rm -rf /var/lib/apt/lists/*

COPY --from=builder /rlib/site-library /usr/local/lib/R/site-library
ENV R_LIBS_SITE=/usr/local/lib/R/site-library

COPY . /api
WORKDIR /api
EXPOSE 8000
CMD ["Rscript", "plumber.R"]
