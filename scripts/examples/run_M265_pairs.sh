#!/usr/bin/env bash
# Copy or edit this file — runs both M265 FASTQ pairs using run_mapseq_one_pair.sh
# From repo root:
#   bash scripts/examples/run_M265_pairs.sh
#
set -eo pipefail
REPO="$(cd "$(dirname "$0")/../.." && pwd)"
RUN="$REPO/scripts/run_mapseq_one_pair.sh"

CONFIG="$REPO/MY_DATA/M265.mapseq.conf"
SAMPLEINFO="$REPO/MY_DATA/M265_sampleinfo.xlsx"

run_one() {
  local out="$1" r1="$2" r2="$3"
  if [[ -n "${DASK_TEMP:-}" ]]; then
    "$RUN" "$CONFIG" "$SAMPLEINFO" "$out" "$r1" "$r2" "$DASK_TEMP"
  else
    "$RUN" "$CONFIG" "$SAMPLEINFO" "$out" "$r1" "$r2"
  fi
}

run_one \
  "$REPO/MY_DATA/M265/M265_1/mapseq_pipeline/M265_1_S1_R1_001" \
  "$REPO/MY_DATA/M265/M265_1/M265_1_S1_R1_001.fastq.gz" \
  "$REPO/MY_DATA/M265/M265_1/M265_1_S1_R2_001.fastq.gz"

run_one \
  "$REPO/MY_DATA/M265/M265_2/mapseq_pipeline/M265_2_S1_R1_001" \
  "$REPO/MY_DATA/M265/M265_2/M265_2_S1_R1_001.fastq.gz" \
  "$REPO/MY_DATA/M265/M265_2/M265_2_S1_R2_001.fastq.gz"

echo "Both M265 pairs finished."
