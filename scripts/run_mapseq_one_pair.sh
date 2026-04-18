#!/usr/bin/env bash
#
# Run the full MAPseq CLI chain for ONE R1/R2 pair — plain bash, no Python batching.
# Activate your conda env first, e.g.: conda activate mapseq_processing
#
# Usage:
#   ./scripts/run_mapseq_one_pair.sh CONFIG SAMPLEINFO OUTDIR R1.fastq.gz R2.fastq.gz [DASK_TEMP_DIR]
#
# OUTDIR is created; intermediate files go under OUTDIR/reads.out, aggregated.out, …
# project_id is read from CONFIG ([project] project_id=... must appear as a line project_id=VALUE).
#
# Example:
#   cd ~/git/mapseq-processing
#   ./scripts/run_mapseq_one_pair.sh \
#     MY_DATA/M265.mapseq.conf \
#     MY_DATA/M265_sampleinfo.xlsx \
#     MY_DATA/M265/M265_1/mapseq_pipeline/M265_1_S1_R1_001 \
#     MY_DATA/M265/M265_1/M265_1_S1_R1_001.fastq.gz \
#     MY_DATA/M265/M265_1/M265_1_S1_R2_001.fastq.gz \
#     ~/scratch
#
# Dask (aggregate_reads only — see README):
#   MAPSEQ_USE_DASK_AGGREGATE=1 uses -k on aggregate_reads (lazy TSV read + disk spill).
#   Temp dir: 6th arg, or env DASK_TEMP, or ~/scratch.
#
set -euo pipefail

REPO="$(cd "$(dirname "$0")/.." && pwd)"
# Conda envs usually expose `python`; macOS system `python3` often lacks deps — prefer `python`.
if [[ -n "${PYTHON:-}" ]]; then
  PY="$PYTHON"
elif command -v python >/dev/null 2>&1; then
  PY=python
else
  PY=python3
fi
SCRIPTS="$REPO/scripts"

CONFIG="${1:?usage: $0 CONFIG SAMPLEINFO OUTDIR R1 R2 [DASK_TEMP]}"
SAMPLEINFO="${2:?}"
OUTDIR="${3:?}"
R1="${4:?}"
R2="${5:?}"
# Optional 6th arg: directory for Dask temp during aggregate_reads (-t), and default spill dir when using MAPSEQ_USE_DASK_AGGREGATE=1
DASK_TEMP_ARG="${6:-}"

if [[ ! -f "$CONFIG" ]]; then echo "missing CONFIG: $CONFIG"; exit 1; fi
if [[ ! -f "$SAMPLEINFO" ]]; then echo "missing SAMPLEINFO: $SAMPLEINFO"; exit 1; fi
if [[ ! -f "$R1" || ! -f "$R2" ]]; then echo "missing FASTQ(s)"; exit 1; fi

PID="$(grep -E '^[[:space:]]*project_id[[:space:]]*=' "$CONFIG" | head -1 | sed 's/.*=[[:space:]]*//;s/[[:space:]]*$//')"
if [[ -z "${PID}" ]]; then
  echo "Could not find project_id= in $CONFIG"
  exit 1
fi

mkdir -p "$OUTDIR"

run() {
  echo ""
  echo ">>> $*"
  "$@"
}

run "$PY" "$SCRIPTS/process_fastq_pairs.py" -v -c "$CONFIG" -L "$OUTDIR/reads.log" \
  -o "$OUTDIR/reads.out/${PID}.reads.tsv" -s "$SAMPLEINFO" "$R1" "$R2"

# aggregate_reads: optional Dask (-k) lowers RAM while aggregating the huge reads TSV.
# process_fastq_pairs itself is still all in-memory pandas (no Dask there).
AGG_EXTRA=()
if [[ "${MAPSEQ_USE_DASK_AGGREGATE:-}" == "1" ]]; then
  AGG_EXTRA=( -k )
  DT="${DASK_TEMP_ARG:-${DASK_TEMP:-${HOME}/scratch}}"
  mkdir -p "$DT"
  AGG_EXTRA+=( -t "$DT" )
elif [[ -n "$DASK_TEMP_ARG" ]]; then
  mkdir -p "$DASK_TEMP_ARG"
  AGG_EXTRA+=( -t "$DASK_TEMP_ARG" )
fi

run "$PY" "$SCRIPTS/aggregate_reads.py" -v -c "$CONFIG" -L "$OUTDIR/aggregated.log" \
  -o "$OUTDIR/aggregated.out/${PID}.aggregated.tsv" \
  "${AGG_EXTRA[@]}" \
  "$OUTDIR/reads.out/${PID}.reads.tsv"

run "$PY" "$SCRIPTS/filter_split.py" -v -c "$CONFIG" -L "$OUTDIR/filtered.log" \
  -o "$OUTDIR/filtered.out/${PID}.filtered.tsv" \
  "$OUTDIR/aggregated.out/${PID}.aggregated.tsv"

run "$PY" "$SCRIPTS/make_readtable.py" -v -c "$CONFIG" -L "$OUTDIR/readtable.log" \
  -o "$OUTDIR/readtable.out/${PID}.readtable.tsv" -s "$SAMPLEINFO" \
  "$OUTDIR/filtered.out/${PID}.filtered.tsv"

run "$PY" "$SCRIPTS/align_collapse.py" -v -c "$CONFIG" -L "$OUTDIR/collapsed.log" \
  -o "$OUTDIR/collapsed.out/${PID}.collapsed.tsv" \
  "$OUTDIR/readtable.out/${PID}.readtable.tsv"

run "$PY" "$SCRIPTS/make_vbctable.py" -v -c "$CONFIG" -L "$OUTDIR/vbctable.log" \
  -o "$OUTDIR/vbctable.out/${PID}.vbctable.tsv" -s "$SAMPLEINFO" \
  "$OUTDIR/collapsed.out/${PID}.collapsed.tsv"

run "$PY" "$SCRIPTS/filter_vbctable.py" -v -c "$CONFIG" -L "$OUTDIR/vbcfiltered.log" \
  -o "$OUTDIR/vbcfiltered.out/${PID}.vbcfiltered.tsv" \
  "$OUTDIR/vbctable.out/${PID}.vbctable.tsv"

run "$PY" "$SCRIPTS/make_matrices.py" -v -c "$CONFIG" -L "$OUTDIR/matrices.log" \
  -O "$OUTDIR/matrices.out" \
  "$OUTDIR/vbcfiltered.out/${PID}.vbcfiltered.tsv"

echo ""
echo "Done. Outputs under: $OUTDIR"
echo "read_count PDFs (if any) are under: $OUTDIR/readtable.out/"
