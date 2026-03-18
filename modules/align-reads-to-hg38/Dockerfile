FROM ubuntu:22.04

ARG SAMTOOLS_VERSION=1.22.1
ARG MINIMAP2_TAG=v2.28

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
       bash \
       ca-certificates \
       gcc \
       git \
       libc6-dev \
       libbz2-dev \
       libcurl4-openssl-dev \
       libdeflate-dev \
       liblzma-dev \
       libncurses5-dev \
       libncursesw5-dev \
       make \
       wget \
       zlib1g-dev \
    && rm -rf /var/lib/apt/lists/*

RUN git clone --branch "${MINIMAP2_TAG}" --depth 1 https://github.com/lh3/minimap2 /tmp/minimap2 \
    && make -C /tmp/minimap2 \
    && cp /tmp/minimap2/minimap2 /usr/local/bin/minimap2 \
    && rm -rf /tmp/minimap2 \
    && minimap2 --version

RUN wget -q "https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2" \
    && tar -xjf "samtools-${SAMTOOLS_VERSION}.tar.bz2" \
    && cd "samtools-${SAMTOOLS_VERSION}" \
    && ./configure --without-curses \
    && make -j"$(nproc)" \
    && make install \
    && samtools --version | head -n 1 \
    && cd / \
    && rm -rf "samtools-${SAMTOOLS_VERSION}" "samtools-${SAMTOOLS_VERSION}.tar.bz2"

SHELL ["/bin/bash", "-euo", "pipefail", "-c"]
