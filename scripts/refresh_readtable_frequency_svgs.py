#!/usr/bin/env python
import argparse
import logging
import sys
from pathlib import Path


def infer_project_id(readtable_path: Path) -> str:
    name = readtable_path.name
    if name.endswith(".readtable.parquet"):
        return name[: -len(".readtable.parquet")]
    if ".readtable." in name:
        return name.split(".readtable.")[0]
    return readtable_path.stem


def parse_args():
    p = argparse.ArgumentParser(
        description="Refresh readtable frequency plot bundle (PDF+SVG) for readtables in a directory."
    )
    p.add_argument(
        "--root",
        default="MY_DATA/readtables_sampleinfos_for_svgs",
        help="Directory containing *.readtable.parquet and corresponding *.readtable.out/ dirs.",
    )
    p.add_argument(
        "--pattern",
        default="*.readtable.parquet",
        help="Glob pattern under root to find readtables (default: *.readtable.parquet).",
    )
    p.add_argument(
        "-v", "--verbose", action="store_true", help="Enable INFO logging."
    )
    p.add_argument(
        "-d", "--debug", action="store_true", help="Enable DEBUG logging."
    )
    return p.parse_args()


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
    from mapseq.plotting import write_readtable_frequency_plots_bundle  # noqa: E402

    root = Path(args.root).expanduser().resolve()
    if not root.exists():
        raise FileNotFoundError(str(root))

    readtables = sorted(root.glob(args.pattern))
    if not readtables:
        logging.error("No readtables found under %s with pattern %s", root, args.pattern)
        return 2

    wrote = 0
    for rt_path in readtables:
        project_id = infer_project_id(rt_path)
        outdir = root / f"{project_id}.readtable.out"
        if not outdir.exists():
            logging.warning("Missing outdir for %s: %s (skipping)", project_id, outdir)
            continue

        logging.info("Loading %s", rt_path)
        df = load_mapseq_df(str(rt_path), fformat="readtable", use_dask=False)
        logging.info("Writing readtable-frequency-plot bundle under %s", outdir)
        write_readtable_frequency_plots_bundle(
            df, str(outdir), project_id, column="read_count", scale="log10"
        )
        wrote += 1

    logging.warning("Refreshed aggregate plots for %s experiment(s).", wrote)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
