#!/usr/bin/env bash
# Unattended / nohup helper: M275 + M277 only (see MY_DATA/M275 and MY_DATA/M277).
# M265 is not run from here.
# Usage:
#   export MAPSEQ_USE_DASK_AGGREGATE=1 DASK_TEMP=~/scratch
#   conda run -n mapseq_processing bash scripts/examples/run_MY_DATA_all_pairs.sh
set -eo pipefail
REPO="$(cd "$(dirname "$0")/../.." && pwd)"
cd "$REPO"
export MAPSEQ_USE_DASK_AGGREGATE="${MAPSEQ_USE_DASK_AGGREGATE:-1}"
export DASK_TEMP="${DASK_TEMP:-${HOME}/scratch}"
mkdir -p "$DASK_TEMP"
echo "=== $(date -u) START M275 + M277 only ==="
bash scripts/examples/run_M275_M277_pairs.sh
echo "=== $(date -u) DONE M275 + M277 ==="
