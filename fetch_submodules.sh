#!/bin/bash
set -euo pipefail

BASE_DIR="$(cd "$(dirname "$0")" && pwd)"
WORK_DIR="$BASE_DIR/work"
ENVS_DIR="$WORK_DIR/envs"
PAVIAN_ENV_DIR="$ENVS_DIR/pavian"

mkdir -p "$BASE_DIR/submodules/10x_saturate"
mkdir -p "$BASE_DIR/submodules/GeneExt"
mkdir -p "$ENVS_DIR"

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

echo "Setting up Pavian Conda environment..."

if [ ! -d "$PAVIAN_ENV_DIR" ]; then
    echo "Creating Conda env at $PAVIAN_ENV_DIR"
    conda env create -p "$PAVIAN_ENV_DIR" -f "$BASE_DIR/modules/local/tools/pavian/environment.yml"
else
    echo "Conda env already exists at $PAVIAN_ENV_DIR"
fi

echo "Installing sankeyD3 into Pavian env..."
conda run -p "$PAVIAN_ENV_DIR" R --vanilla -e '
if (!requireNamespace("remotes", quietly = TRUE) && !requireNamespace("devtools", quietly = TRUE)) {
  install.packages("remotes", repos = "https://cloud.r-project.org")
}
if (requireNamespace("remotes", quietly = TRUE)) {
  remotes::install_github("fbreitwieser/sankeyD3", upgrade = "never", dependencies = TRUE)
} else {
  devtools::install_github("fbreitwieser/sankeyD3", upgrade = "never", dependencies = TRUE)
}
'

echo "Pavian environment is ready at $PAVIAN_ENV_DIR"
