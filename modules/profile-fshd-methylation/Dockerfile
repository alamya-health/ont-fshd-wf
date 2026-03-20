FROM ubuntu:22.04

ARG SAMTOOLS_VERSION=1.22.1
ARG MODKIT_VERSION=v0.6.1
ARG MODKIT_ASSET=modkit_${MODKIT_VERSION}_u16_x86_64.tar.gz
ARG MODKIT_URL=https://github.com/nanoporetech/modkit/releases/download/${MODKIT_VERSION}/${MODKIT_ASSET}

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
       bash \
       bzip2 \
       ca-certificates \
       curl \
       gcc \
       libc6-dev \
       libbz2-dev \
       libcurl4-openssl-dev \
       libdeflate-dev \
       liblzma-dev \
       libncurses5-dev \
       libncursesw5-dev \
       make \
       procps \
       python3 \
       python3-pip \
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

RUN set -eu; \
    work=/tmp/modkit_install; mkdir -p "$work"; \
    curl -fL "$MODKIT_URL" -o "$work/modkit.tgz"; \
    tar -xzf "$work/modkit.tgz" -C "$work"; \
    BIN="$(find "$work" -type f -name modkit -perm -u+x | head -n1)"; \
    test -n "$BIN"; \
    install -m 0755 "$BIN" /usr/local/bin/modkit; \
    rm -rf "$work"; \
    modkit --version

RUN pip3 install --no-cache-dir pysam

SHELL ["/bin/bash", "-euo", "pipefail", "-c"]
