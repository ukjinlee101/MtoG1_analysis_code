#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 2 ]]; then
  echo "Usage: $0 <input.cool> <output.cool>"
  exit 1
fi

input_file="$1"
output_file="$2"

eval "$(mamba shell hook --shell bash)"
mamba activate hicexplorer_3.7.4

# run the correction
hicCorrectMatrix correct \
  -m "$input_file" \
  --correctionMethod KR \
  --outFileName "$output_file"