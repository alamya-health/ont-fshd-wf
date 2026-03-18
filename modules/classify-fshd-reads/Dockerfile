FROM ubuntu:22.04

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
       bash \
       ca-certificates \
       r-base \
       wget \
    && rm -rf /var/lib/apt/lists/*

RUN /usr/bin/Rscript -e "options(repos='https://cloud.r-project.org'); install.packages('remotes')" \
    && /usr/bin/Rscript -e "remotes::install_version('dplyr', version='1.1.4'); remotes::install_version('tidyr', version='1.3.1')"

SHELL ["/bin/bash", "-euo", "pipefail", "-c"]
