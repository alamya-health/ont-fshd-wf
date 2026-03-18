FROM ubuntu:22.04

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
       bash \
       ca-certificates \
       python3 \
       python3-pip \
    && rm -rf /var/lib/apt/lists/*

RUN pip3 install --no-cache-dir pysam

SHELL ["/bin/bash", "-euo", "pipefail", "-c"]
