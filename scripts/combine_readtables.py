#!/usr/bin/env python
import argparse
import logging
import os
import shutil
import sys
from typing import Dict, List

import pandas as pd

gitpath = os.path.expanduser("~/git/mapseq-processing")
sys.path.append(gitpath)

from mapseq.core import fix_mapseq_df_types, load_mapseq_df, write_mapseq_df  # noqa: E402
from mapseq.plotting import write_readtable_frequency_plots_bundle  # noqa: E402


GROUP_KEY = ["rtprimer", "label", "type", "vbc_read", "umi", "libtag"]


def parse_args():
    p = argparse.ArgumentParser(
        description="Combine two MAPseq readtables (e.g. resequenced halves) by summing read_count on strict identity key."
    )
    p.add_argument("readtable_a", help="First input readtable (.tsv or .parquet).")
    p.add_argument("readtable_b", help="Second input readtable (.tsv or .parquet).")
    p.add_argument(
        "-o",
        "--outdir",
        required=True,
        help="Output directory (e.g. MY_DATA/readtables_sampleinfos_for_svgs/M265_combined.readtable.out).",
    )
    p.add_argument(
        "--project-id",
        default="M265",
        help="Project id for output filenames (default: M265).",
    )
    p.add_argument(
        "--sampleinfo",
        default=None,
        help="Optional sampleinfo .xlsx/.tsv to copy into outdir for provenance.",
    )
    p.add_argument(
        "--write-tsv",
        action="store_true",
        help="Also write TSV alongside parquet (default writes parquet only).",
    )
    p.add_argument(
        "--make-plots",
        action="store_true",
        help="Also regenerate aggregate read_count frequency PDF+SVG in outdir.",
    )
    p.add_argument("-v", "--verbose", action="store_true", help="Enable INFO logging.")
    p.add_argument("-d", "--debug", action="store_true", help="Enable DEBUG logging.")
    return p.parse_args()


def _required_cols_present(df: pd.DataFrame, required: List[str]) -> None:
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f"Missing required columns: {missing}")


def _build_agg_map(df: pd.DataFrame, group_key: List[str]) -> Dict[str, str]:
    agg: Dict[str, str] = {"read_count": "sum"}
    for c in df.columns:
        if c in group_key or c == "read_count":
            continue
        agg[c] = "first"
    return agg


def main() -> int:
    args = parse_args()
    level = logging.WARN
    if args.verbose:
        level = logging.INFO
    if args.debug:
        level = logging.DEBUG
    logging.basicConfig(
        format="%(asctime)s [%(levelname)s] %(filename)s:%(lineno)d %(funcName)s(): %(message)s",
        level=level,
    )

    outdir = os.path.abspath(args.outdir)
    os.makedirs(outdir, exist_ok=True)

    logging.info("Loading readtable A: %s", args.readtable_a)
    a = load_mapseq_df(args.readtable_a, fformat="readtable", use_dask=False)
    logging.info("Loading readtable B: %s", args.readtable_b)
    b = load_mapseq_df(args.readtable_b, fformat="readtable", use_dask=False)

    _required_cols_present(a, GROUP_KEY + ["read_count"])
    _required_cols_present(b, GROUP_KEY + ["read_count"])

    sum_a = int(pd.to_numeric(a["read_count"], errors="coerce").fillna(0).sum())
    sum_b = int(pd.to_numeric(b["read_count"], errors="coerce").fillna(0).sum())
    logging.warning("read_count sums: A=%s B=%s total_expected=%s", sum_a, sum_b, sum_a + sum_b)

    combined = pd.concat([a, b], ignore_index=True, copy=False)
    combined["read_count"] = pd.to_numeric(combined["read_count"], errors="coerce").fillna(0).astype("int64")

    agg_map = _build_agg_map(combined, GROUP_KEY)
    logging.info("Grouping by %s (%s rows before)", GROUP_KEY, len(combined))
    combined = combined.groupby(GROUP_KEY, observed=True, dropna=False).agg(agg_map).reset_index()
    combined = fix_mapseq_df_types(combined, fformat="readtable")
    logging.warning("Combined rows=%s read_count_sum=%s", len(combined), int(combined["read_count"].sum()))

    base_out = os.path.join(outdir, f"{args.project_id}.readtable")
    outformats = ["parquet"]
    if args.write_tsv:
        outformats = ["tsv", "parquet"]
    write_mapseq_df(combined, base_out, outformats=outformats)

    if args.sampleinfo is not None:
        dst = os.path.join(outdir, os.path.basename(args.sampleinfo))
        shutil.copy2(args.sampleinfo, dst)
        logging.info("Copied sampleinfo to %s", dst)

    if args.make_plots:
        write_readtable_frequency_plots_bundle(
            combined,
            outdir,
            args.project_id,
            column="read_count",
            scale="log10",
        )
        logging.info("Wrote readtable-frequency-plot bundle under %s", outdir)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())

