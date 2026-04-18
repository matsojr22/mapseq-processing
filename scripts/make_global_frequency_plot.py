#!/usr/bin/env python
import argparse
import logging
import math
import os
import sys
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd
import seaborn as sns


def parse_args():
    p = argparse.ArgumentParser(
        description="Create a pooled global uniform-scale logrank readtable frequency plot."
    )
    p.add_argument(
        "-O",
        "--outdir",
        default="MY_DATA/readtables_sampleinfos_for_svgs/GLOBAL.readtable.out",
        help="Output directory for global plot bundle.",
    )
    p.add_argument(
        "--project-id",
        default="GLOBAL",
        help="Project id used in output filenames (default: GLOBAL).",
    )
    p.add_argument(
        "--inputs",
        nargs="+",
        default=[
            "MY_DATA/readtables_sampleinfos_for_svgs/M265_combined.readtable.out/M265.readtable.parquet",
            "MY_DATA/readtables_sampleinfos_for_svgs/M275.readtable.parquet",
            "MY_DATA/readtables_sampleinfos_for_svgs/M277.readtable.parquet",
            "MY_DATA/readtables_sampleinfos_for_svgs/M300.readtable.parquet",
            "MY_DATA/readtables_sampleinfos_for_svgs/M312.readtable.parquet",
        ],
        help="Readtable parquet/tsv inputs to pool (default: m265_combined + m275 + m277 + m300 + m312).",
    )
    p.add_argument(
        "--column",
        default="read_count",
        help="Column to use for the distribution (default: read_count).",
    )
    p.add_argument(
        "--points",
        type=int,
        default=12000,
        help="Number of log-spaced rank points to sample (default: 12000).",
    )
    p.add_argument(
        "--y-cutoff",
        type=float,
        default=None,
        help="Optional horizontal cutoff line to draw at read_count (e.g. 2).",
    )
    p.add_argument(
        "--cutoff-suffix",
        default=None,
        help="If set along with --y-cutoff, write additional *_<suffix> outputs (e.g. cutoff2).",
    )
    p.add_argument("-v", "--verbose", action="store_true", help="Enable INFO logging.")
    p.add_argument("-d", "--debug", action="store_true", help="Enable DEBUG logging.")
    return p.parse_args()


def pooled_value_counts(series_list: List[pd.Series]) -> pd.Series:
    """
    Return pooled value_counts across multiple series.
    Index: value (e.g. read_count), Value: frequency (#rows with that value).
    """
    pooled: Dict[int, int] = {}
    for ser in series_list:
        vc = ser.value_counts(dropna=False)
        for val, freq in vc.items():
            ival = int(val)
            pooled[ival] = pooled.get(ival, 0) + int(freq)
    out = pd.Series(pooled, dtype="int64")
    out.sort_index(ascending=False, inplace=True)
    return out


def sample_rank_curve_from_value_counts(vc_desc: pd.Series, ranks: np.ndarray) -> np.ndarray:
    """
    Given value_counts for y-values sorted descending by y (index is y),
    return y at each requested 1-based rank using cumulative frequencies.
    """
    if vc_desc.empty:
        return np.array([], dtype=float)
    freqs = vc_desc.to_numpy(dtype=np.int64)
    values = vc_desc.index.to_numpy(dtype=np.int64)
    cumfreq = np.cumsum(freqs)
    # ranks are 1..N inclusive
    idx = np.searchsorted(cumfreq, ranks, side="left")
    idx = np.clip(idx, 0, len(values) - 1)
    return values[idx].astype(float)

def infer_run_label(path_str: str) -> str:
    p = Path(path_str)
    name = p.name
    if name.endswith(".readtable.parquet"):
        return name[: -len(".readtable.parquet")]
    if name.endswith(".readtable.tsv"):
        return name[: -len(".readtable.tsv")]
    if name.endswith(".parquet"):
        return name[: -len(".parquet")]
    if name.endswith(".tsv"):
        return name[: -len(".tsv")]
    if p.parent.name.endswith(".readtable.out") and p.name == "M265.readtable.parquet":
        return p.parent.name[: -len(".readtable.out")]
    return p.stem


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

    repo_root = Path(__file__).resolve().parents[1]
    sys.path.insert(0, str(repo_root))
    from mapseq.core import load_mapseq_df  # noqa: E402
    from mapseq.plotting import make_logticks, setup_svg_font_defaults  # noqa: E402

    outdir = Path(args.outdir).resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    inputs = [str(Path(p).resolve()) for p in args.inputs]
    sources_path = outdir / f"{args.project_id}.sources.txt"
    sources_path.write_text("\n".join(inputs) + "\n")

    series_list: List[pd.Series] = []
    per_run: List[Tuple[str, pd.Series]] = []
    total_rows = 0
    total_read_mass = 0
    for p in inputs:
        logging.info("Loading %s", p)
        df = load_mapseq_df(p, fformat="readtable", use_dask=False)
        if args.column not in df.columns:
            raise ValueError(f"Missing column {args.column!r} in {p}")
        ser = pd.to_numeric(df[args.column], errors="coerce").fillna(0).astype("int64")
        series_list.append(ser)
        per_run.append((infer_run_label(p), ser))
        total_rows += int(len(ser))
        total_read_mass += int(ser.sum())

    vc_global = pooled_value_counts(series_list)
    pooled_n = int(vc_global.sum())
    pooled_sum = int((vc_global.index.to_numpy(dtype=np.int64) * vc_global.to_numpy(dtype=np.int64)).sum())
    logging.warning("pooled N=%s (sum rows=%s)", pooled_n, total_rows)
    logging.warning("pooled sum(%s)=%s (sum rows sum=%s)", args.column, pooled_sum, total_read_mass)

    if pooled_n <= 0:
        raise ValueError("No data to plot (pooled N <= 0).")

    points = max(200, int(args.points))
    ranks = np.unique(np.round(np.logspace(0, math.log10(pooled_n), num=points)).astype(np.int64))
    ranks = ranks[(ranks >= 1) & (ranks <= pooled_n)]
    y = sample_rank_curve_from_value_counts(vc_global, ranks)

    # Plot: uniform-scale logrank (x log, y log)
    import matplotlib.pyplot as plt  # noqa: E402
    from matplotlib.backends.backend_pdf import PdfPages  # noqa: E402

    setup_svg_font_defaults()
    pdf_out = outdir / f"{args.project_id}.readtable-frequency-plot_uniform-scale_logrank.pdf"
    svg_out = outdir / f"{args.project_id}.readtable-frequency-plot_uniform-scale_logrank.svg"
    pdf_overlay = outdir / f"{args.project_id}.readtable-frequency-plot_uniform-scale_logrank_overlay.pdf"
    svg_overlay = outdir / f"{args.project_id}.readtable-frequency-plot_uniform-scale_logrank_overlay.svg"
    pdf_out_cut = None
    svg_out_cut = None
    pdf_overlay_cut = None
    svg_overlay_cut = None
    if args.y_cutoff is not None and args.cutoff_suffix is not None:
        pdf_out_cut = outdir / f"{args.project_id}.readtable-frequency-plot_uniform-scale_logrank_{args.cutoff_suffix}.pdf"
        svg_out_cut = outdir / f"{args.project_id}.readtable-frequency-plot_uniform-scale_logrank_{args.cutoff_suffix}.svg"
        pdf_overlay_cut = outdir / f"{args.project_id}.readtable-frequency-plot_uniform-scale_logrank_overlay_{args.cutoff_suffix}.pdf"
        svg_overlay_cut = outdir / f"{args.project_id}.readtable-frequency-plot_uniform-scale_logrank_overlay_{args.cutoff_suffix}.svg"

    with PdfPages(str(pdf_out)) as pdfpages:
        fig, ax = plt.subplots(figsize=(8, 6))
        fig.suptitle(f"{args.project_id} pooled read_count frequency (uniform log rank)")
        ax.plot(ranks.astype(float), y, linewidth=1.25)
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel("Rank (log scale)")
        ax.set_ylabel(f"log10( {args.column} )")
        ax.set_yticks(make_logticks(float(y.max()) if len(y) else 1.0))
        ax.set_xlim(left=max(0.8, float(ranks.min()) * 0.9), right=None)
        ax.set_ylim(bottom=0, top=None)
        if args.y_cutoff is not None:
            ax.axhline(float(args.y_cutoff), color="red", linestyle="--", linewidth=1.0)
        ax.text(
            float(ranks.max()),
            float(y.max()) if len(y) else 1.0,
            s=f"n={pooled_n}\nsum={pooled_sum}",
            fontsize=11,
            horizontalalignment="right",
            verticalalignment="top",
        )
        pdfpages.savefig(fig)
        fig.savefig(str(svg_out), format="svg")
        plt.close(fig)

    if pdf_out_cut is not None and svg_out_cut is not None:
        with PdfPages(str(pdf_out_cut)) as pdfpages:
            fig, ax = plt.subplots(figsize=(8, 6))
            fig.suptitle(
                f"{args.project_id} pooled read_count frequency (uniform log rank; cutoff={args.y_cutoff})"
            )
            ax.plot(ranks.astype(float), y, linewidth=1.25)
            ax.set_xscale("log")
            ax.set_yscale("log")
            ax.set_xlabel("Rank (log scale)")
            ax.set_ylabel(f"log10( {args.column} )")
            ax.set_yticks(make_logticks(float(y.max()) if len(y) else 1.0))
            ax.set_xlim(left=max(0.8, float(ranks.min()) * 0.9), right=None)
            ax.set_ylim(bottom=0, top=None)
            ax.axhline(float(args.y_cutoff), color="red", linestyle="--", linewidth=1.0)
            ax.text(
                float(ranks.max()),
                float(y.max()) if len(y) else 1.0,
                s=f"n={pooled_n}\nsum={pooled_sum}",
                fontsize=11,
                horizontalalignment="right",
                verticalalignment="top",
            )
            pdfpages.savefig(fig)
            fig.savefig(str(svg_out_cut), format="svg")
            plt.close(fig)

    # Overlay per-run curves (same axis conventions)
    try:
        palette = sns.color_palette("Set2", n_colors=len(per_run))
    except Exception:
        palette = None

    with PdfPages(str(pdf_overlay)) as pdfpages:
        fig, ax = plt.subplots(figsize=(8, 6))
        fig.suptitle(f"{args.project_id} per-run read_count frequency (uniform log rank)")
        for i, (label, ser) in enumerate(per_run):
            vc = ser.value_counts(dropna=False)
            vc.index = vc.index.astype(int)
            vc = vc.sort_index(ascending=False)
            n_i = int(vc.sum())
            if n_i <= 0:
                continue
            ranks_i = np.unique(
                np.round(np.logspace(0, math.log10(n_i), num=points)).astype(np.int64)
            )
            ranks_i = ranks_i[(ranks_i >= 1) & (ranks_i <= n_i)]
            y_i = sample_rank_curve_from_value_counts(vc, ranks_i)
            color = palette[i] if palette is not None else None
            ax.plot(ranks_i.astype(float), y_i, linewidth=1.15, label=label, color=color)

        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel("Rank (log scale)")
        ax.set_ylabel(f"log10( {args.column} )")
        ax.set_yticks(make_logticks(float(y.max()) if len(y) else 1.0))
        ax.set_xlim(left=0.8, right=None)
        ax.set_ylim(bottom=0, top=None)
        if args.y_cutoff is not None:
            ax.axhline(float(args.y_cutoff), color="red", linestyle="--", linewidth=1.0)
        ax.legend(loc="best", fontsize=9, frameon=False)
        pdfpages.savefig(fig)
        fig.savefig(str(svg_overlay), format="svg")
        plt.close(fig)

    if pdf_overlay_cut is not None and svg_overlay_cut is not None:
        with PdfPages(str(pdf_overlay_cut)) as pdfpages:
            fig, ax = plt.subplots(figsize=(8, 6))
            fig.suptitle(
                f"{args.project_id} per-run read_count frequency (uniform log rank; cutoff={args.y_cutoff})"
            )
            for i, (label, ser) in enumerate(per_run):
                vc = ser.value_counts(dropna=False)
                vc.index = vc.index.astype(int)
                vc = vc.sort_index(ascending=False)
                n_i = int(vc.sum())
                if n_i <= 0:
                    continue
                ranks_i = np.unique(
                    np.round(np.logspace(0, math.log10(n_i), num=points)).astype(np.int64)
                )
                ranks_i = ranks_i[(ranks_i >= 1) & (ranks_i <= n_i)]
                y_i = sample_rank_curve_from_value_counts(vc, ranks_i)
                color = palette[i] if palette is not None else None
                ax.plot(
                    ranks_i.astype(float),
                    y_i,
                    linewidth=1.15,
                    label=label,
                    color=color,
                )

            ax.set_xscale("log")
            ax.set_yscale("log")
            ax.set_xlabel("Rank (log scale)")
            ax.set_ylabel(f"log10( {args.column} )")
            ax.set_yticks(make_logticks(float(y.max()) if len(y) else 1.0))
            ax.set_xlim(left=0.8, right=None)
            ax.set_ylim(bottom=0, top=None)
            ax.axhline(float(args.y_cutoff), color="red", linestyle="--", linewidth=1.0)
            ax.legend(loc="best", fontsize=9, frameon=False)
            pdfpages.savefig(fig)
            fig.savefig(str(svg_overlay_cut), format="svg")
            plt.close(fig)

    logging.warning("Wrote %s", pdf_out)
    logging.warning("Wrote %s", svg_out)
    logging.warning("Wrote %s", pdf_overlay)
    logging.warning("Wrote %s", svg_overlay)
    if pdf_out_cut is not None and svg_out_cut is not None:
        logging.warning("Wrote %s", pdf_out_cut)
        logging.warning("Wrote %s", svg_out_cut)
    if pdf_overlay_cut is not None and svg_overlay_cut is not None:
        logging.warning("Wrote %s", pdf_overlay_cut)
        logging.warning("Wrote %s", svg_overlay_cut)
    logging.warning("Wrote %s", sources_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

