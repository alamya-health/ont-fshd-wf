#!/usr/bin/env bash
set -euo pipefail

AWS_REGION_DEFAULT="us-west-1"
PIPELINE_REPO_DEFAULT="ont-fshd-wf"

ACCOUNT_ID=""
REGION="${AWS_REGION:-$AWS_REGION_DEFAULT}"
TAG="${TAG:-$(date +%Y%m%d)-$(git rev-parse --short HEAD 2>/dev/null || echo local)}"
PIPELINE_REPO="${PIPELINE_REPO:-$PIPELINE_REPO_DEFAULT}"
MODULE_MAP_FILE="${MODULE_MAP_FILE:-scripts/module_tool_versions.map}"
ECR_LIFECYCLE=1
ROLE_CHECK=1
PLATFORM_AUTO=1
PLATFORM=""
PUSH_LATEST=0
BUILD_ARGS=()

usage() {
  cat <<USAGE
Usage:
  $0 --account-id 123456789012 [options]

Options:
  --account-id <ID>      Required AWS account ID
  --region <REGION>      AWS region (default: ${AWS_REGION_DEFAULT})
  --tag <TAG>            Image tag suffix (default: auto)
  --repo <NAME>          ECR repository name (default: ${PIPELINE_REPO_DEFAULT})
  --module-map <FILE>    Module tag map (default: scripts/module_tool_versions.map)
  --platform <PLAT>      Docker platform (default: auto)
  --latest               Also push <module>-latest
  --build-arg KEY=VAL    Pass-through docker build arg
  --skip-lifecycle       Skip lifecycle policy setup
  --no-role-check        Skip assumed-role validation
  -h | --help            Show help
USAGE
  exit "${1:-0}"
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --account-id) ACCOUNT_ID="${2:-}"; shift 2 ;;
    --region) REGION="${2:-}"; shift 2 ;;
    --tag) TAG="${2:-}"; shift 2 ;;
    --repo) PIPELINE_REPO="${2:-}"; shift 2 ;;
    --module-map) MODULE_MAP_FILE="${2:-}"; shift 2 ;;
    --platform) PLATFORM_AUTO=0; PLATFORM="${2:-}"; shift 2 ;;
    --latest) PUSH_LATEST=1; shift ;;
    --build-arg) BUILD_ARGS+=("$2"); shift 2 ;;
    --skip-lifecycle) ECR_LIFECYCLE=0; shift ;;
    --no-role-check) ROLE_CHECK=0; shift ;;
    -h|--help) usage 0 ;;
    *) echo "Unknown option: $1"; usage 1 ;;
  esac
done

command -v aws >/dev/null || { echo "ERROR: aws CLI not found"; exit 1; }
command -v docker >/dev/null || { echo "ERROR: docker not found"; exit 1; }

if [[ -z "$ACCOUNT_ID" ]]; then
  echo "ERROR: --account-id is required."
  usage 1
fi

if [[ $PLATFORM_AUTO -eq 1 ]]; then
  case "$(uname -m || echo unknown)" in
    arm64|aarch64) PLATFORM="linux/amd64" ;;
    x86_64|amd64) PLATFORM="linux/amd64" ;;
    *) PLATFORM="linux/amd64" ;;
  esac
fi

AWS_ARN="$(aws sts get-caller-identity --query Arn --output text 2>/dev/null || true)"
if [[ -z "$AWS_ARN" ]]; then
  echo "ERROR: Unable to query AWS STS identity."
  exit 1
fi
if [[ $ROLE_CHECK -eq 1 && "$AWS_ARN" != *":assumed-role/"* ]]; then
  echo "ERROR: Builds should use an assumed role."
  echo "Current principal: $AWS_ARN"
  exit 1
fi

ECR_REGISTRY="${ACCOUNT_ID}.dkr.ecr.${REGION}.amazonaws.com"

LIFECYCLE_JSON='{
  "rules": [{
    "rulePriority": 1,
    "description": "Keep last 200 images",
    "selection": {"tagStatus":"any","countType":"imageCountMoreThan","countNumber":200},
    "action": {"type":"expire"}
  }]
}'

aws ecr get-login-password --region "${REGION}" \
| docker login --username AWS --password-stdin "${ECR_REGISTRY}"

if ! aws ecr describe-repositories --repository-names "${PIPELINE_REPO}" --region "${REGION}" >/dev/null 2>&1; then
  aws ecr create-repository \
    --repository-name "${PIPELINE_REPO}" \
    --image-tag-mutability IMMUTABLE \
    --image-scanning-configuration scanOnPush=true \
    --region "${REGION}" >/dev/null
  if [[ "$ECR_LIFECYCLE" -ne 0 ]]; then
    aws ecr put-lifecycle-policy \
      --repository-name "${PIPELINE_REPO}" \
      --lifecycle-policy-text "${LIFECYCLE_JSON}" \
      --region "${REGION}" >/dev/null || true
  fi
fi

declare -A MODULE_TAG_SUFFIX=()
if [[ -f "${MODULE_MAP_FILE}" ]]; then
  while IFS= read -r line; do
    line="${line%%#*}"
    line="$(echo -n "$line" | awk '{$1=$1};1')"
    [[ -z "$line" ]] && continue
    key="${line%%=*}"
    val="${line#*=}"
    MODULE_TAG_SUFFIX["$key"]="$val"
  done < "${MODULE_MAP_FILE}"
fi

mapfile -d '' MODULE_DIRS < <(find modules -mindepth 1 -maxdepth 1 -type d -print0 2>/dev/null || true)
declare -a MODULES=()
for d in "${MODULE_DIRS[@]:-}"; do
  [[ -f "${d}/Dockerfile" ]] || continue
  MODULES+=("$(basename "$d")")
done

if [[ ${#MODULES[@]} -eq 0 ]]; then
  echo "No module Dockerfiles found."
  exit 0
fi

for mod in "${MODULES[@]}"; do
  suffix="${MODULE_TAG_SUFFIX[$mod]:-$TAG}"
  image="${ECR_REGISTRY}/${PIPELINE_REPO}:${mod}-${suffix}"
  echo "Building ${image}"
  build_cmd=(
    docker buildx build
    --platform "${PLATFORM}"
    --push
    -t "${image}"
  )
  for arg in "${BUILD_ARGS[@]:-}"; do
    build_cmd+=(--build-arg "$arg")
  done
  build_cmd+=("modules/${mod}")
  "${build_cmd[@]}"

  if [[ "$PUSH_LATEST" -eq 1 ]]; then
    latest_image="${ECR_REGISTRY}/${PIPELINE_REPO}:${mod}-latest"
    docker buildx imagetools create -t "${latest_image}" "${image}"
  fi
done
