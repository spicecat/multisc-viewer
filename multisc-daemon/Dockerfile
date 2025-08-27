FROM rstudio/plumber:latest

LABEL org.opencontainers.image.title="multisc-daemon" \
	org.opencontainers.image.description="Plotting daemon for Multisc-Viewer"

RUN apt-get update -qq \
	&& apt-get install -y --no-install-recommends \
	libeigen3-dev \
	libglpk40 \
	libxml2 \
	&& rm -rf /var/lib/apt/lists/*

RUN install2.r \
	--ncpus -1 --error --skipinstalled \
	--repos 'https://cloud.r-project.org' \
	--repos 'https://satijalab.r-universe.dev' \
	--repos 'https://bnprks.r-universe.dev' \
	Seurat \
	rjson \
	# BPCells \
	# presto \
	# glmGamPoi \
	&& rm -rf /tmp/downloaded_packages/ /tmp/*.rds

COPY . /api
CMD ["/api/plumber.R"]
