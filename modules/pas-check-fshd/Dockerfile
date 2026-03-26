FROM ubuntu:22.04

ARG PYSAM_VERSION=0.23.3

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
       bash \
       ca-certificates \
       procps \
       python3 \
       python3-pip \
    && rm -rf /var/lib/apt/lists/*

RUN pip3 install --no-cache-dir "pysam==${PYSAM_VERSION}"

SHELL ["/bin/bash", "-euo", "pipefail", "-c"]
