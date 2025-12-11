#!/bin/bash
set -e

# Directory of the script
BASE_DIR="$(cd "$(dirname "$0")" && pwd)"

echo "Fetching submodules..."

# 10x_saturate
if [ ! -d "$BASE_DIR/10x_saturate" ]; then
    echo "Downloading 10x_saturate..."
    curl -L -o /tmp/10x_saturate.zip https://github.com/zolotarovgl/10x_saturate/archive/refs/heads/main.zip
    unzip /tmp/10x_saturate.zip -d /tmp
    mv /tmp/10x_saturate-main "$BASE_DIR/submodules/10x_saturate"
    rm /tmp/10x_saturate.zip
fi

# GeneExt
if [ ! -d "$BASE_DIR/GeneExt" ]; then
    echo "Downloading GeneExt..."
    curl -L -o /tmp/GeneExt.zip https://github.com/zolotarovgl/GeneExt/archive/refs/heads/main.zip
    unzip /tmp/GeneExt.zip -d /tmp
    mv /tmp/GeneExt-main "$BASE_DIR/submodules/GeneExt"
    rm /tmp/GeneExt.zip
fi

echo "Submodules downloaded successfully!"
