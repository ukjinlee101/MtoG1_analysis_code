#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 2 ]]; then
  echo "Usage: $0 <input.cool> <output.cool>"
  exit 1
fi

cool_file="$1"
num_cpus="$2"

eval "$(mamba shell hook --shell bash)"
mamba activate hicexplorer_3.7.4

# run the correction

# recompute KR weights into bins.weight
cooler balance --force -p "${num_cpus}" \
    "${cool_file}"