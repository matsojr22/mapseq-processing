#!/usr/bin/env bash
# Full MAPseq chain for MY_DATA M275 and M277 (one FASTQ pair each).
# FASTQs: MY_DATA/M275/M275_EK_S1_R*.fastq.gz, MY_DATA/M277/M277_MJ_S1_R*.fastq.gz
# conda activate mapseq_processing && bash scripts/examples/run_M275_M277_pairs.sh
# Fresh rerun: rm -rf MY_DATA/M275/mapseq_pipeline/* MY_DATA/M277/mapseq_pipeline/*
set -eo pipefail
REPO="$(cd "$(dirname "$0")/../.." && pwd)"
RUN="$REPO/scripts/run_mapseq_one_pair.sh"

run_one() {
  local cfg="$1" si="$2" out="$3" r1="$4" r2="$5"
  if [[ -n "${DASK_TEMP:-}" ]]; then
    "$RUN" "$cfg" "$si" "$out" "$r1" "$r2" "$DASK_TEMP"
  else
    "$RUN" "$cfg" "$si" "$out" "$r1" "$r2"
  fi
}

run_one \
  "$REPO/MY_DATA/M275.mapseq.conf" \
  "$REPO/MY_DATA/M275_sampleinfo.xlsx" \
  "$REPO/MY_DATA/M275/mapseq_pipeline/M275_EK_S1_R1_001" \
  "$REPO/MY_DATA/M275/M275_EK_S1_R1_001.fastq.gz" \
  "$REPO/MY_DATA/M275/M275_EK_S1_R2_001.fastq.gz"

run_one \
  "$REPO/MY_DATA/M277.mapseq.conf" \
  "$REPO/MY_DATA/M277_sampleinfo.xlsx" \
  "$REPO/MY_DATA/M277/mapseq_pipeline/M277_MJ_S1_R1_001" \
  "$REPO/MY_DATA/M277/M277_MJ_S1_R1_001.fastq.gz" \
  "$REPO/MY_DATA/M277/M277_MJ_S1_R2_001.fastq.gz"

echo "M275 and M277 pipelines finished."
