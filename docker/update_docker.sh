#!/bin/sh
# Build and publish one release. VERSION is the tag without the leading "v".
set -eu

if [ "$#" -ne 1 ] || [ -z "$1" ]; then
    echo "Usage: sh $0 <version-without-v>" >&2
    exit 1
fi
VERSION=$1
case "$VERSION" in
    *[!a-zA-Z0-9_.-]* | [!0-9]*)
        echo "Error: version must start with a digit and contain only letters, digits, ., _, or -" >&2
        exit 1
        ;;
esac

SCRIPT_DIR=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
IMAGE=naechyun/leviosam2
ENV_IMAGE=naechyun/leviosam2_env:$VERSION

docker build --build-arg "VERSION=$VERSION" --tag "$ENV_IMAGE" "$SCRIPT_DIR/env"
docker build --build-arg "ENV_IMAGE=$ENV_IMAGE" \
    --tag "$IMAGE:$VERSION" --tag "$IMAGE:latest" "$SCRIPT_DIR/production"

docker push "$IMAGE:$VERSION"
docker push "$IMAGE:latest"
