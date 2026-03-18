#!/usr/bin/env bash
set -euo pipefail

./scripts/Build_Module_Dockerfiles_ECR_macbook.sh \
  --account-id 160433694210 \
  --region ap-southeast-1 \
  --repo ont-fshd-wf \
  --build-arg CFLAGS="-O2 -pipe -march=x86-64 -mtune=generic" \
  --build-arg CXXFLAGS="-O2 -pipe -march=x86-64 -mtune=generic" \
  --no-role-check
