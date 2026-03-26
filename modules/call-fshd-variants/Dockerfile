FROM hkubal/clair3:v2.0.0

USER root
ENV DEBIAN_FRONTEND=noninteractive

ARG SNIFFLES_VERSION=2.7.0
ARG SNPEFF_VERSION=v4_3t
ARG SNPEFF_URL=https://downloads.sourceforge.net/project/snpeff/snpEff_${SNPEFF_VERSION}_core.zip

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
       ca-certificates \
       curl \
       openjdk-17-jre-headless \
       unzip \
       wget \
    && rm -rf /var/lib/apt/lists/*

RUN python3 -m pip install --no-cache-dir "sniffles==${SNIFFLES_VERSION}"

RUN mkdir -p /opt/snpeff \
    && curl -fL "${SNPEFF_URL}" -o /tmp/snpeff.zip \
    && unzip -q /tmp/snpeff.zip -d /tmp \
    && cp -a /tmp/snpEff/* /opt/snpeff/ \
    && rm -rf /tmp/snpeff.zip /tmp/snpEff \
    && test -f /opt/snpeff/snpEff.jar \
    && test -f /opt/snpeff/SnpSift.jar

SHELL ["/bin/bash", "-euo", "pipefail", "-c"]
