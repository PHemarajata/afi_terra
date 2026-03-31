#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 2 || $# -gt 3 ]]; then
  echo "Usage: $0 <dockerhub-user-or-repo> <tag> [platform]" >&2
  echo "Example: $0 phemarajata614 0.1.1 linux/amd64" >&2
  exit 1
fi

image_repo="$1/afi-terra"
image_tag="$2"
platform="${3:-linux/amd64}"
image_ref="${image_repo}:${image_tag}"

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

echo "Building ${image_ref} for ${platform} from ${repo_root}/Dockerfile"
docker build --platform "${platform}" -t "${image_ref}" "${repo_root}"

echo "Pushing ${image_ref}"
docker push "${image_ref}"

echo
echo "Use this Terra override value in the batch JSON:"
echo "${image_ref}"