FROM ubuntu:22.04

ARG DPLYR_VERSION=1.1.4
ARG TIDYR_VERSION=1.3.1
ARG STRINGI_VERSION=1.8.7
ARG CPP11_VERSION=0.5.3
ARG STRINGR_VERSION=1.6.0
ARG PURRR_VERSION=1.2.1

ENV DEBIAN_FRONTEND=noninteractive
ENV DPLYR_VERSION=${DPLYR_VERSION}
ENV TIDYR_VERSION=${TIDYR_VERSION}
ENV STRINGI_VERSION=${STRINGI_VERSION}
ENV CPP11_VERSION=${CPP11_VERSION}
ENV STRINGR_VERSION=${STRINGR_VERSION}
ENV PURRR_VERSION=${PURRR_VERSION}

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
       bash \
       ca-certificates \
       procps \
       r-base \
       r-base-dev \
       libcurl4-openssl-dev \
       libssl-dev \
       libxml2-dev \
       wget \
    && rm -rf /var/lib/apt/lists/*

RUN Rscript -e 'options(repos=c(CRAN="https://cloud.r-project.org")); install.packages("remotes")' \
    && Rscript -e 'remotes::install_version("dplyr", version=Sys.getenv("DPLYR_VERSION"), repos="https://cloud.r-project.org", upgrade="never")' \
    && Rscript -e 'remotes::install_version("stringi", version=Sys.getenv("STRINGI_VERSION"), repos="https://cloud.r-project.org", upgrade="never"); remotes::install_version("cpp11", version=Sys.getenv("CPP11_VERSION"), repos="https://cloud.r-project.org", upgrade="never"); remotes::install_version("stringr", version=Sys.getenv("STRINGR_VERSION"), repos="https://cloud.r-project.org", upgrade="never"); remotes::install_version("purrr", version=Sys.getenv("PURRR_VERSION"), repos="https://cloud.r-project.org", upgrade="never")' \
    && Rscript -e 'remotes::install_version("tidyr", version=Sys.getenv("TIDYR_VERSION"), repos="https://cloud.r-project.org", dependencies=FALSE, upgrade="never")' \
    && Rscript -e 'stopifnot(as.character(utils::packageVersion("dplyr")) == Sys.getenv("DPLYR_VERSION")); stopifnot(as.character(utils::packageVersion("tidyr")) == Sys.getenv("TIDYR_VERSION")); library(dplyr); library(tidyr); cat("R packages OK\n")'

SHELL ["/bin/bash", "-euo", "pipefail", "-c"]
