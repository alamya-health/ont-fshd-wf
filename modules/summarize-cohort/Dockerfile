FROM ubuntu:22.04

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
       bash \
       ca-certificates \
       python3 \
    && rm -rf /var/lib/apt/lists/*

SHELL ["/bin/bash", "-euo", "pipefail", "-c"]
