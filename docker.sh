#!/usr/bin/env bash

set -x
set -eo pipefail

docker build --pull --rm -f "Dockerfile" -t bsaintjo/cawlr:full "."
>&2 echo "Image successfully built"

docker push bsaintjo/cawlr:full
>&2 echo "Image successfully pushed"