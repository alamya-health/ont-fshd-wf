FROM mambaorg/micromamba:1.5.10

USER root
ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
       bash \
       ca-certificates \
       curl \
       openjdk-17-jre-headless \
       unzip \
       wget \
    && rm -rf /var/lib/apt/lists/*

RUN micromamba install -y -n base -c conda-forge -c bioconda \
      python=3.11 \
      clair3 \
      whatshap \
      sniffles \
      samtools \
      bcftools \
      tabix \
    && micromamba clean --all --yes

RUN mkdir -p /opt/snpeff \
    && curl -fL https://snpeff-public.s3.amazonaws.com/versions/snpEff_latest_core.zip -o /tmp/snpeff.zip \
    && unzip -q /tmp/snpeff.zip -d /tmp \
    && cp -a /tmp/snpEff/* /opt/snpeff/ \
    && rm -rf /tmp/snpeff.zip /tmp/snpEff \
    && test -f /opt/snpeff/snpEff.jar \
    && test -f /opt/snpeff/SnpSift.jar

SHELL ["/bin/bash", "-euo", "pipefail", "-c"]
