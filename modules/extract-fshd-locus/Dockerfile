FROM ubuntu:22.04

ARG SAMTOOLS_VERSION=1.22.1

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
       bash \
       ca-certificates \
       gcc \
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
