#!/usr/bin/env python
import argparse
import logging
import math
import os
import re
import sys
from configparser import ConfigParser
from typing import Dict, Optional, Tuple

import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.backends.backend_pdf import PdfPages

DEFAULT_SPECS = [
    (1, 1000000),
    (2, 1000000),
    (2, None),
]


def natural_sort_key(value: str):
    text = str(value)
    parts = re.split(r"(\d+)", text)
    key = []
    for part in parts:
        if part.isdigit():
            key.append(int(part))
        else:
            key.append(part.lower())
    return key


def setup_svg_font_defaults() -> None:
    """
    Keep SVG text editable and request Helvetica as default.
    """
    plt.rcParams["svg.fonttype"] = "none"
    plt.rcParams["font.family"] = "sans-serif"
    plt.rcParams["font.sans-serif"] = ["Helvetica", "Arial", "DejaVu Sans"]


def normalize_rtprimer(value):
    if pd.isna(value):
        return None
    text = str(value).strip()
    if text == "":
        return None

    norm = text.upper().replace(" ", "")
    if norm.startswith("BC"):
        norm = norm[2:]

    try:
        return str(int(float(norm)))
    except ValueError:
        return norm


def sanitize_filename(text: str) -> str:
    clean = re.sub(r"[^A-Za-z0-9._-]+", "-", text.strip())
    clean = re.sub(r"-{2,}", "-", clean).strip("-")
    return clean or "unknown"


def infer_project_id(readtable_path: str, config_path: Optional[str]) -> str:
    if config_path:
        cp = ConfigParser()
        cp.read(config_path)
        if cp.has_option("project", "project_id"):
            return cp.get("project", "project_id")

    base = os.path.basename(readtable_path)
    if ".readtable" in base:
        return base.split(".readtable")[0]
    return base.split(".")[0]


def load_sample_info_compat(file_name: str, sheet_name: str = "Sample information") -> pd.DataFrame:
    """
    Parse sampleinfo in the same spirit as mapseq.core.load_sample_info,
    but without importing pipeline modules.
    """
    sheet_to_sample = {
        "Tube # by user": "usertube",
        "Our Tube #": "ourtube",
        "Sample names provided by user": "samplename",
        "Site information": "siteinfo",
        "Spike-in Ratio": "si_ratio",
        "Read Count Minimum": "min_reads",
        "RT primers for MAPseq": "rtprimer",
        "Brain": "brain",
        "Region": "region",
        "Matrix Column": "matrixcolumn",
    }

    sample_columns = [
        "usertube",
        "ourtube",
        "samplename",
        "siteinfo",
        "si_ratio",
        "rtprimer",
        "brain",
        "region",
        "matrixcolumn",
        "min_reads",
    ]

    if file_name.endswith(".xlsx"):
        edf = pd.read_excel(file_name, sheet_name=sheet_name, header=1, dtype=str)
        sdf = pd.DataFrame()
        for ecol in edf.columns:
            ecol_stp = str(ecol).strip()
            if ecol_stp in sheet_to_sample:
                sdf[sheet_to_sample[ecol_stp]] = edf[ecol]
        sdf = sdf[sdf["rtprimer"].isna() == False]
        sdf.fillna("", inplace=True)
        for scol in sample_columns:
            if scol not in sdf.columns:
                if scol == "samplename":
                    sdf[scol] = sdf.get("ourtube", "")
                elif scol == "region":
                    sdf[scol] = sdf.get("rtprimer", "")
                elif scol == "si_ratio":
                    sdf[scol] = 1.0
                elif scol == "min_reads":
                    sdf[scol] = 1
                else:
                    sdf[scol] = ""
        sdf.loc[sdf["min_reads"] == "", "min_reads"] = 1
        sdf.loc[sdf["si_ratio"] == "", "si_ratio"] = 1.0
    elif file_name.endswith(".tsv"):
        sdf = pd.read_csv(file_name, sep="\t", index_col=0, keep_default_na=False, dtype=str, comment="#")
        for scol in sample_columns:
            if scol not in sdf.columns:
                sdf[scol] = ""
    else:
        raise ValueError("sampleinfo must be .xlsx or .tsv")

    sdf["min_reads"] = pd.to_numeric(sdf["min_reads"], errors="coerce").fillna(1).astype(int)
    sdf["si_ratio"] = pd.to_numeric(sdf["si_ratio"], errors="coerce").fillna(1.0).astype(float)
    return sdf


def load_readtable(infile: str) -> pd.DataFrame:
    if infile.endswith(".parquet"):
        df = pd.read_parquet(infile)
    elif infile.endswith(".tsv"):
        df = pd.read_csv(infile, sep="\t", index_col=0)
    else:
        raise ValueError("readtable must be .tsv or .parquet")

    required = {"label", "rtprimer", "read_count"}
    missing = required.difference(set(df.columns))
    if missing:
        raise ValueError(f"readtable missing required columns: {sorted(missing)}")

    # Handle categorical/object/string columns safely by converting first.
    df["label"] = df["label"].astype("string").fillna("").str.strip()
    df["rtprimer"] = df["rtprimer"].astype("string").fillna("").str.strip()
    df["read_count"] = pd.to_numeric(df["read_count"], errors="coerce").fillna(0).astype(int)
    return df


def build_rt_to_human_map(sampdf: pd.DataFrame) -> Dict[str, str]:
    if "rtprimer" not in sampdf.columns:
        raise ValueError("sampleinfo is missing required column: rtprimer")
    if "samplename" not in sampdf.columns:
        raise ValueError("sampleinfo is missing required column: samplename")
    if "brain" not in sampdf.columns:
        raise ValueError("sampleinfo is missing required column: brain")

    working = sampdf.copy()
    working["_rt_norm"] = working["rtprimer"].map(normalize_rtprimer)
    working["brain"] = working["brain"].fillna("").astype(str).str.strip()
    working["samplename"] = working["samplename"].fillna("").astype(str).str.strip()

    mapping: Dict[str, str] = {}
    for _, row in working.iterrows():
        rt_norm = row["_rt_norm"]
        if not rt_norm:
            continue

        brain = row["brain"]
        area = row["samplename"]
        if brain and area:
            human = f"{brain} - {area}"
        elif brain:
            human = brain
        elif area:
            human = area
        else:
            continue

        if rt_norm in mapping and mapping[rt_norm] != human:
            logging.warning(
                "Conflicting sampleinfo mapping for rtprimer=%s: '%s' vs '%s'. Keeping first.",
                rt_norm,
                mapping[rt_norm],
                human,
            )
            continue
        mapping[rt_norm] = human

    return mapping


def build_label_to_human_map(readtable_df: pd.DataFrame, rt_to_human: Dict[str, str]) -> Dict[str, str]:
    working = readtable_df.copy()
    working["_rt_norm"] = working["rtprimer"].map(normalize_rtprimer)

    mapping: Dict[str, str] = {}
    grouped = working.groupby("label", dropna=False)
    for label, gdf in grouped:
        if not label:
            continue

        rt_counts = gdf["_rt_norm"].dropna().value_counts()
        if rt_counts.empty:
            logging.warning("No rtprimer found for label=%s; leaving label unchanged.", label)
            mapping[label] = label
            continue

        if len(rt_counts) > 1:
            logging.warning(
                "Multiple rtprimers mapped to label=%s (%s). Using most frequent=%s.",
                label,
                ", ".join(rt_counts.index.astype(str).tolist()),
                rt_counts.index[0],
            )
        rt_norm = str(rt_counts.index[0])
        mapping[label] = rt_to_human.get(rt_norm, label)
        if mapping[label] == label:
            logging.warning(
                "No sampleinfo match for label=%s (rtprimer=%s); leaving BC label.",
                label,
                rt_norm,
            )

    return mapping


def filter_for_plot(df: pd.DataFrame, groupby: str, column: str, min_count: int) -> pd.DataFrame:
    fdf = df.copy()
    if min_count > 1:
        fdf = fdf[fdf[column] >= int(min_count)]
    fdf = fdf[fdf[groupby] != ""]
    fdf = fdf.reset_index(drop=True)
    return fdf


def make_logticks(max_value: int):
    ticklist = []
    i = 1
    while i < max_value:
        ticklist.append(i)
        i *= 10
    ticklist.append(i)
    return ticklist


def calc_freq_threshold(df: pd.DataFrame, fraction: float = 0.9, column: str = "read_count") -> int:
    ser = df[column].copy()
    ser.sort_values(ascending=False, inplace=True)
    ser.reset_index(drop=True, inplace=True)
    idx = int(len(ser) * fraction)
    return int(ser.iloc[idx] + 1)


def counts_axis_plot(
    ax,
    df: pd.DataFrame,
    column: str,
    title: str,
    nranks: Optional[int] = None,
    proportion: float = 0.85,
) -> None:
    pdf = df.sort_values(by=[column], ascending=False)
    h = calc_freq_threshold(pdf, fraction=proportion, column=column)
    nranks_initial = len(pdf)
    if nranks is not None:
        pdf = pdf.iloc[:nranks]
    pdf = pdf.reset_index(drop=True).reset_index()

    s = int(pdf[column].sum())
    n = len(df)
    t = int(pdf[column].max())

    ax.plot(pdf["index"], pdf[column], linewidth=1.25)
    lx = int(pdf["index"].max())
    ly = int(pdf[column].max())
    ax.set_xlabel("Rank")
    ax.set_yscale("log")
    ax.set_yticks(make_logticks(ly))
    ax.set_ylabel(f"log10( {column} )")
    plot_title = f"{title} log10()"

    ax.set_xlim(left=0, right=None)
    ax.set_ylim(bottom=0, top=None)
    ax.text(
        lx,
        ly,
        s=f"n={n}\ntop={t}\nsum={s}\nest_{proportion}_threshold={h}",
        fontsize=11,
        horizontalalignment="right",
        verticalalignment="top",
    )
    if nranks is not None:
        plot_title = f"{plot_title} nranks=[{nranks}/{nranks_initial}]"
    ax.set_title(plot_title, fontsize=10)


def write_overall_pdf(df: pd.DataFrame, outfile: str, column: str) -> None:
    with PdfPages(outfile) as pdfpages:
        fig, ax = plt.subplots(figsize=(8, 6))
        fig.suptitle("Overall read count frequency (log10)")
        counts_axis_plot(ax, df, column=column, title=f"{column} freqplot", nranks=None)
        pdfpages.savefig(fig)
        plt.close(fig)


def write_combined_pdf_and_individual_svgs(
    df: pd.DataFrame,
    outdir: str,
    svg_dir: str,
    project_id: str,
    groupby: str,
    column: str,
    min_count: int,
    nranks: Optional[int],
    label_display: Dict[str, str],
) -> Tuple[str, int]:
    ranks = "all" if nranks is None else str(nranks)
    outfile = os.path.join(
        outdir,
        f"{project_id}.all.{column}.by{groupby}.c{min_count}.r{ranks}.pdf",
    )
    title = f"{project_id}:all {column} frequency (log10)"
    groups = sorted(list(df[groupby].unique()), key=natural_sort_key)
    page_dims = (11.7, 8.27)
    plots_per_page = 9
    num_figs = max(1, int(math.ceil(float(len(groups)) / float(plots_per_page))))

    figlist = []
    axlist = []
    for _ in range(num_figs):
        fig, axes = plt.subplots(nrows=3, ncols=3, figsize=page_dims, layout="constrained")
        fig.suptitle(title)
        figlist.append(fig)
        for ax in axes.flat:
            axlist.append(ax)

    svg_count = 0
    for i, group in enumerate(groups):
        gdf = df[df[groupby] == group]
        panel_label = label_display.get(group, group)
        panel_title = f"{panel_label} {column}"
        counts_axis_plot(
            axlist[i],
            gdf,
            column=column,
            title=panel_title,
            nranks=nranks,
        )

        ind_fig, ind_ax = plt.subplots(figsize=(8, 6), layout="constrained")
        ind_fig.suptitle(title)
        counts_axis_plot(
            ind_ax,
            gdf,
            column=column,
            title=panel_title,
            nranks=nranks,
        )
        svg_name = (
            f"{sanitize_filename(panel_label)}__{sanitize_filename(str(group))}"
            f".{column}.c{min_count}.r{ranks}.svg"
        )
        ind_fig.savefig(os.path.join(svg_dir, svg_name), format="svg")
        plt.close(ind_fig)
        svg_count += 1

    # Hide unused axes for a cleaner PDF page.
    for j in range(len(groups), len(axlist)):
        axlist[j].axis("off")

    with PdfPages(outfile) as pdfpages:
        for fig in figlist:
            pdfpages.savefig(fig)
            plt.close(fig)

    return outfile, svg_count


def parse_args():
    parser = argparse.ArgumentParser(
        description="Regenerate readtable rank plots with human-readable titles and per-panel SVGs."
    )
    parser.add_argument("readtable", type=str, help="Input readtable (.tsv or .parquet).")
    parser.add_argument("-s", "--sampleinfo", required=True, type=str, help="Sample info (.xlsx or .tsv).")
    parser.add_argument(
        "-S",
        "--samplesheet",
        default="Sample information",
        type=str,
        help="Sheet tab name for sampleinfo xlsx.",
    )
    parser.add_argument("-O", "--outdir", default=None, type=str, help="Output directory. Defaults to readtable dir.")
    parser.add_argument("-c", "--config", default=None, type=str, help="Optional config file for project_id.")
    parser.add_argument("--groupby", default="label", type=str, help="Grouping column in readtable.")
    parser.add_argument("--column", default="read_count", type=str, help="Count column to plot.")
    parser.add_argument("--svg-subdir", default="svgs", type=str, help="Subdirectory under outdir for individual SVG files.")
    parser.add_argument("-d", "--debug", action="store_true", help="Enable debug logging.")
    parser.add_argument("-v", "--verbose", action="store_true", help="Enable info logging.")
    return parser.parse_args()


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

    setup_svg_font_defaults()

    outdir = os.path.abspath(args.outdir) if args.outdir else os.path.dirname(os.path.abspath(args.readtable))
    svg_dir = os.path.join(outdir, args.svg_subdir)
    os.makedirs(outdir, exist_ok=True)
    os.makedirs(svg_dir, exist_ok=True)

    project_id = infer_project_id(args.readtable, args.config)
    logging.info("project_id=%s", project_id)

    readtable_df = load_readtable(args.readtable)
    sample_df = load_sample_info_compat(args.sampleinfo, args.samplesheet)

    rt_to_human = build_rt_to_human_map(sample_df)
    label_display = build_label_to_human_map(readtable_df, rt_to_human)

    if args.groupby not in readtable_df.columns:
        raise ValueError(f"groupby column not found in readtable: {args.groupby}")
    if args.column not in readtable_df.columns:
        raise ValueError(f"count column not found in readtable: {args.column}")

    total_svgs = 0
    for min_count, nranks in DEFAULT_SPECS:
        plot_df = filter_for_plot(readtable_df, args.groupby, args.column, min_count=min_count)
        if len(plot_df) == 0:
            logging.warning("No rows remain for min_count=%s; skipping.", min_count)
            continue
        pdf_file, svg_count = write_combined_pdf_and_individual_svgs(
            plot_df,
            outdir=outdir,
            svg_dir=svg_dir,
            project_id=project_id,
            groupby=args.groupby,
            column=args.column,
            min_count=min_count,
            nranks=nranks,
            label_display=label_display,
        )
        total_svgs += svg_count
        logging.info("Wrote PDF: %s", pdf_file)

    repo_root = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
    if repo_root not in sys.path:
        sys.path.insert(0, repo_root)
    from mapseq.plotting import write_readtable_frequency_plots_bundle

    write_readtable_frequency_plots_bundle(
        readtable_df, outdir, project_id, column=args.column, scale="log10"
    )
    logging.info("Wrote readtable-frequency-plot bundle under %s", outdir)
    logging.info("Wrote %s individual SVGs under %s", total_svgs, svg_dir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
