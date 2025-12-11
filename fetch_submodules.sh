#!/bin/bash
set -e

BASE_DIR="$(cd "$(dirname "$0")" && pwd)"

echo "Fetching submodules..."

# 10x_saturate
if [ ! -f "$BASE_DIR/submodules/10x_saturate/README.md" ]; then
    echo "Downloading 10x_saturate..."
    curl -L -o /tmp/10x_saturate.zip https://github.com/zolotarovgl/10x_saturate/archive/refs/heads/main.zip
    unzip -q /tmp/10x_saturate.zip -d /tmp
    cp -r /tmp/10x_saturate-main/* "$BASE_DIR/submodules/10x_saturate/"
    rm -rf /tmp/10x_saturate-main /tmp/10x_saturate.zip
fi

# GeneExt
if [ ! -f "$BASE_DIR/submodules/GeneExt/environment.yaml" ]; then
    echo "Downloading GeneExt..."
    curl -L -o /tmp/GeneExt.zip https://github.com/zolotarovgl/GeneExt/archive/refs/heads/main.zip
    unzip -q /tmp/GeneExt.zip -d /tmp
    cp -r /tmp/GeneExt-main/* "$BASE_DIR/submodules/GeneExt/"
    rm -rf /tmp/GeneExt-main /tmp/GeneExt.zip
fi

echo "Submodules downloaded successfully!"
