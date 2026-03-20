#!/usr/bin/env bash
set -euo pipefail

# -------------------------------------------------------------
# Single-Repo AWS ECR builder for Nextflow modules
# (Apple Silicon host -> build linux/amd64 -> push to ECR)
# Mirrors the working pattern used in sibling repos.
# -------------------------------------------------------------

AWS_REGION_DEFAULT="ap-southeast-1"
PIPELINE_REPO_DEFAULT="ont-fshd-wf"

ACCOUNT_ID=""
REGION="${AWS_REGION:-$AWS_REGION_DEFAULT}"
TAG="${TAG:-$(date +%Y%m%d)-$(git rev-parse --short HEAD 2>/dev/null || echo local)}"
PIPELINE_REPO="${PIPELINE_REPO:-$PIPELINE_REPO_DEFAULT}"
MODULE_MAP_FILE="${MODULE_MAP_FILE:-module_map.txt}"
ECR_LIFECYCLE=1
ROLE_CHECK=1
PUSH_LATEST=0
PLATFORM_AUTO=1
PLATFORM=""
BUILD_ARGS=()

usage() {
  cat <<USAGE
Single-Repo AWS ECR builder for Nextflow modules (Apple Silicon -> linux/amd64)

Usage:
  $0 --account-id 123456789012 [options]

Options:
  --account-id <ID>      (required) 12-digit AWS account ID
  --region <REGION>      AWS region (default: ${AWS_REGION_DEFAULT}) or env AWS_REGION
  --tag <TAG>            image tag suffix (default: auto: YYYYMMDD-<git short sha>|local)
  --repo <NAME>          ECR repository name (default: ${PIPELINE_REPO_DEFAULT})
  --module-map <FILE>    alias map file (default: module_map.txt)
  --platform <PLAT>      docker platform (default: auto; Apple Silicon -> linux/amd64)
  --latest               also push <module>-latest tags
  --build-arg KEY=VAL    pass through to 'docker buildx build' (repeatable)
  --skip-lifecycle       skip lifecycle policy on the repo
  --no-role-check        do not require assumed-role identity (not recommended)
  -h | --help            show this help
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
    --skip-lifecycle) ECR_LIFECYCLE=0; shift ;;
    --no-role-check) ROLE_CHECK=0; shift ;;
    --latest) PUSH_LATEST=1; shift ;;
    --platform) PLATFORM_AUTO=0; PLATFORM="${2:-}"; shift 2 ;;
    --build-arg) BUILD_ARGS+=("$2"); shift 2 ;;
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
if ! [[ "$ACCOUNT_ID" =~ ^[0-9]{12}$ ]]; then
  echo "ERROR: --account-id must be a 12-digit number."
  exit 1
fi
if [[ -z "$REGION" ]]; then
  echo "ERROR: --region or AWS_REGION must be set."
  exit 1
fi
if [[ -z "$PIPELINE_REPO" ]]; then
  echo "ERROR: --repo must not be empty."
  exit 1
fi

if [[ $PLATFORM_AUTO -eq 1 ]]; then
  host_arch="$(uname -m || echo unknown)"
  case "$host_arch" in
    arm64|aarch64) PLATFORM="linux/amd64" ;;
    x86_64|amd64)  PLATFORM="linux/amd64" ;;
    *)             PLATFORM="linux/amd64" ;;
  esac
fi

AWS_ARN="$(aws sts get-caller-identity --query Arn --output text 2>/dev/null || true)"
if [[ -z "$AWS_ARN" ]]; then
  echo "ERROR: Unable to query AWS STS identity. Check your credentials."
  exit 1
fi
if [[ $ROLE_CHECK -eq 1 && "$AWS_ARN" != *":assumed-role/"* ]]; then
  echo "ERROR: Builds should use an assumed role."
  echo "Current principal: $AWS_ARN"
  echo "(Pass --no-role-check to override, not recommended.)"
  exit 1
fi

ECR_REGISTRY="${ACCOUNT_ID}.dkr.ecr.${REGION}.amazonaws.com"

echo "Region:        ${REGION}"
echo "Account:       ${ACCOUNT_ID}"
echo "Registry:      ${ECR_REGISTRY}"
echo "Repository:    ${PIPELINE_REPO}"
echo "Tag suffix:    ${TAG}"
echo "Platform:      ${PLATFORM}"
echo

LIFECYCLE_JSON='{
  "rules": [{
    "rulePriority": 1,
    "description": "Keep last 200 images; expire older",
    "selection": {"tagStatus":"any","countType":"imageCountMoreThan","countNumber":200},
    "action": {"type":"expire"}
  }]
}'

aws ecr get-login-password --region "${REGION}" \
| docker login --username AWS --password-stdin "${ECR_REGISTRY}"

ensure_repo() {
  local repo="$1"
  if ! aws ecr describe-repositories --repository-names "${repo}" --region "${REGION}" >/dev/null 2>&1; then
    echo "Creating ECR repo: ${repo}"
    aws ecr create-repository \
      --repository-name "${repo}" \
      --image-tag-mutability IMMUTABLE \
      --image-scanning-configuration scanOnPush=true \
      --region "${REGION}" >/dev/null
    if [[ "$ECR_LIFECYCLE" -ne 0 ]]; then
      aws ecr put-lifecycle-policy \
        --repository-name "${repo}" \
        --lifecycle-policy-text "${LIFECYCLE_JSON}" \
        --region "${REGION}" >/dev/null || true
    fi
  fi
}
ensure_repo "${PIPELINE_REPO}"

declare -A MODULE_TO_CONTEXT=()
if [[ -f "${MODULE_MAP_FILE}" ]]; then
  echo "Using module map: ${MODULE_MAP_FILE}"
  while IFS= read -r line; do
    line="${line%%#*}"
    line="$(echo -n "$line" | awk '{$1=$1};1')"
    [[ -z "$line" ]] && continue
    key="${line%%=*}"
    val="${line#*=}"
    key="$(echo -n "$key" | awk '{$1=$1};1')"
    val="$(echo -n "$val" | awk '{$1=$1};1')"
    [[ -z "$key" || -z "$val" ]] && continue
    MODULE_TO_CONTEXT["$key"]="$val"
  done < "${MODULE_MAP_FILE}"
  echo
fi

echo "Scanning modules/ for Dockerfiles..."
mapfile -d '' MODULE_DIRS < <(find modules -mindepth 1 -maxdepth 1 -type d -print0 2>/dev/null || true)
declare -a MODULES=()
for d in "${MODULE_DIRS[@]:-}"; do
  [[ -f "${d}/Dockerfile" ]] || continue
  MODULES+=("$(basename "$d")")
done
if [[ ${#MODULES[@]} -eq 0 ]]; then
  echo "No modules found (no Dockerfiles under modules/*). Nothing to do."
  exit 0
fi
IFS=$'\n' MODULES=($(sort <<<"${MODULES[*]}"))
unset IFS
echo "Found modules:"
printf '  - %s\n' "${MODULES[@]}"
echo

resolve_context() {
  local mod="$1"
  if [[ -n "${MODULE_TO_CONTEXT[$mod]+set}" ]]; then
    echo "${MODULE_TO_CONTEXT[$mod]}"
  else
    echo "modules/${mod}"
  fi
}

for mod in "${MODULES[@]}"; do
  ctx="$(resolve_context "$mod")"
  image="${ECR_REGISTRY}/${PIPELINE_REPO}:${mod}-${TAG}"

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
  build_cmd+=("${ctx}")
  "${build_cmd[@]}"

  if [[ "$PUSH_LATEST" -eq 1 ]]; then
    latest_image="${ECR_REGISTRY}/${PIPELINE_REPO}:${mod}-latest"
    docker buildx imagetools create -t "${latest_image}" "${image}"
  fi
done

echo
echo "Done."
echo "Images are available under: ${ECR_REGISTRY}/${PIPELINE_REPO}:<module>-${TAG}"
