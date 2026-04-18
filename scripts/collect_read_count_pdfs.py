#!/usr/bin/env python
"""
Discover MAPseq experiment folders under MY_DATA/, optionally run the MAPseq
pipeline **per R1/R2 pair** by invoking each stage script in order (same argv
pattern as scripts/process_all.py, but **without** calling process_all.py).
Each stage is run with **-v** (verbose), like the README examples, not **-d** (debug).

Stages: process_fastq_pairs → aggregate_reads → filter_split → make_readtable →
align_collapse → make_vbctable → filter_vbctable → make_matrices.

Each pair writes under ``<version>/<pipeline_subdir>/<pair_slug>/`` so outputs
never clobber each other. *read_count*.pdf files (from make_readtable) are copied
to ``<my-data-root>/pdfs/<experiment>/<version>/<pair_slug>/``.

Default is dry-run: no subprocesses are started unless you pass --execute.

**Run with your Conda env** (dependencies are not on the system Python), e.g.::

  conda activate mapseq_processing
  python scripts/collect_read_count_pdfs.py --my-data-root MY_DATA --execute

Or: ``conda run -n mapseq_processing python scripts/collect_read_count_pdfs.py ...``

Requires: sampleinfo .xlsx where needed (reads / readtable / vbctable stages).
Missing FASTQ pairs are skipped (e.g. incomplete downloads).

Usage:
  conda activate mapseq_processing
  python scripts/collect_read_count_pdfs.py --my-data-root MY_DATA
  python scripts/collect_read_count_pdfs.py --my-data-root MY_DATA --execute
  python scripts/collect_read_count_pdfs.py --my-data-root MY_DATA --halt readtable --execute
"""
from __future__ import annotations

import argparse
import csv
import logging
import os
import re
import shutil
import subprocess
import sys
from configparser import ConfigParser
from pathlib import Path
from typing import Iterable, List, Optional, Sequence, Tuple

# Repo root: .../mapseq-processing
_SCRIPT_DIR = Path(__file__).resolve().parent
_REPO_ROOT = _SCRIPT_DIR.parent
_DEFAULT_CONFIG = _REPO_ROOT / "etc" / "mapseq.conf"
_PROCESS_FASTQ = _SCRIPT_DIR / "process_fastq_pairs.py"

# Mirror scripts/process_all.py (do not import process_all — it pulls heavy deps).
STEPLIST: List[str] = [
    "reads",
    "aggregated",
    "filtered",
    "readtable",
    "collapsed",
    "vbctable",
    "vbcfiltered",
    "matrices",
]
STEPMAP = {
    "reads": "process_fastq_pairs",
    "aggregated": "aggregate_reads",
    "filtered": "filter_split",
    "readtable": "make_readtable",
    "collapsed": "align_collapse",
    "vbctable": "make_vbctable",
    "vbcfiltered": "filter_vbctable",
    "matrices": "make_matrices",
}
DIRMAP = {
    "process_fastq_pairs": "reads",
    "aggregate_reads": "aggregated",
    "filter_split": "filtered",
    "make_readtable": "readtable",
    "align_collapse": "collapsed",
    "make_vbctable": "vbctable",
    "filter_vbctable": "vbcfiltered",
    "make_matrices": "matrices",
}


def _setup_logging(verbose: bool) -> None:
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format="%(asctime)s [%(levelname)s] %(message)s",
    )


def _read_project_id(config_path: Path) -> str:
    cp = ConfigParser()
    read = cp.read(config_path)
    if not read:
        raise ValueError(f"Could not read config: {config_path}")
    return cp.get("project", "project_id").strip()


def _steps_for_halt(halt: Optional[str]) -> List[str]:
    if halt is None:
        return list(STEPLIST)
    if halt not in STEPLIST:
        raise ValueError(f"--halt must be one of {STEPLIST}, got {halt!r}")
    i = STEPLIST.index(halt)
    return STEPLIST[: i + 1]


def _pair_slug(r1: Path) -> str:
    """Filesystem-safe label from R1 filename (one directory per pair)."""
    name = r1.name
    for suf in (".fastq.gz", ".fq.gz", ".fastq", ".fq", ".FASTQ.GZ", ".FQ.GZ"):
        if name.endswith(suf):
            name = name[: -len(suf)]
            break
    safe = re.sub(r"[^\w.\-]+", "_", name)
    safe = safe.strip("._")
    return safe or "pair"


def _pipeline_command(
    step: str,
    *,
    step_index: int,
    steps: Sequence[str],
    outdir: Path,
    config: Path,
    sampleinfo: Path,
    project_id: str,
    r1: Path,
    r2: Path,
    dask_temp: Optional[Path],
) -> List[str]:
    sprog = STEPMAP[step]
    sname = DIRMAP[sprog]
    script = _SCRIPT_DIR / f"{sprog}.py"
    outdir_res = outdir.resolve()
    log_file = outdir_res / f"{step}.log"
    cmd: List[str] = [
        sys.executable,
        str(script),
        "-v",
        "-c",
        str(config.resolve()),
        "-L",
        str(log_file),
    ]

    if sname == "matrices":
        soutdir = outdir_res / f"{sname}.out"
        prev_sname = DIRMAP[STEPMAP[steps[step_index - 1]]]
        infile = outdir_res / f"{prev_sname}.out" / f"{project_id}.{prev_sname}.tsv"
        cmd.extend(["-O", str(soutdir)])
        cmd.append(str(infile))
        return cmd

    soutfile = outdir_res / f"{sname}.out" / f"{project_id}.{sname}.tsv"
    cmd.extend(["-o", str(soutfile)])

    if sname == "reads":
        cmd.extend(
            [
                "-s",
                str(sampleinfo.resolve()),
                str(r1.resolve()),
                str(r2.resolve()),
            ]
        )
    elif sname in ("readtable", "vbctable"):
        prev_sname = DIRMAP[STEPMAP[steps[step_index - 1]]]
        infile = outdir_res / f"{prev_sname}.out" / f"{project_id}.{prev_sname}.tsv"
        cmd.extend(["-s", str(sampleinfo.resolve()), str(infile)])
    else:
        prev_sname = DIRMAP[STEPMAP[steps[step_index - 1]]]
        infile = outdir_res / f"{prev_sname}.out" / f"{project_id}.{prev_sname}.tsv"
        cmd.append(str(infile))
        if sprog == "aggregate_reads" and dask_temp is not None:
            inf = cmd.pop()
            cmd.extend(["-k", "-t", str(dask_temp.resolve()), inf])

    return cmd


def _run_pipeline_for_pair(
    *,
    config: Path,
    sampleinfo: Path,
    pair_outdir: Path,
    r1: Path,
    r2: Path,
    halt: Optional[str],
    dask_temp: Optional[Path],
    execute: bool,
) -> int:
    """Run (or print) the pipeline for a single R1/R2 pair. Returns last exit code (0 = ok)."""
    try:
        project_id = _read_project_id(config)
    except Exception as exc:
        logging.error("Could not read [project] project_id from %s: %s", config, exc)
        return 1

    steps = _steps_for_halt(halt)
    if execute:
        pair_outdir.mkdir(parents=True, exist_ok=True)

    for step_idx, step in enumerate(steps):
        cmd = _pipeline_command(
            step,
            step_index=step_idx,
            steps=steps,
            outdir=pair_outdir,
            config=config,
            sampleinfo=sampleinfo,
            project_id=project_id,
            r1=r1,
            r2=r2,
            dask_temp=dask_temp,
        )
        logging.info("Would run:\n  %s", " \\\n  ".join(cmd))

        if not execute:
            continue

        if step_idx > 0:
            infile = Path(cmd[-1])
            if not infile.is_file():
                logging.error("Missing expected input for step %s: %s", step, infile)
                return 1

        proc = subprocess.run(cmd, cwd=str(_REPO_ROOT))
        if proc.returncode != 0:
            logging.error("Step %s failed with code %s", step, proc.returncode)
            return int(proc.returncode)

    return 0


def _find_fastq_r1_files(directory: Path) -> List[Path]:
    """Non-recursive R1 FASTQ paths."""
    out: List[Path] = []
    for pat in ("*R1*.fastq.gz", "*R1*.fastq", "*_R1_*.fastq.gz", "*_R1_*.fastq"):
        out.extend(directory.glob(pat))
    # de-dupe, stable sort
    return sorted(set(out), key=lambda p: p.name.lower())


def _r1_to_r2_path(r1: Path) -> Optional[Path]:
    name = r1.name
    if "_R1_" in name:
        r2_name = name.replace("_R1_", "_R2_", 1)
    elif ".R1." in name:
        r2_name = name.replace(".R1.", ".R2.", 1)
    else:
        # e.g. ...R1_001... without underscore pair
        m = re.search(r"_R1([._])", name)
        if m:
            r2_name = name.replace("_R1" + m.group(1), "_R2" + m.group(1), 1)
        else:
            return None
    return r1.with_name(r2_name)


def _pair_fastqs(r1_files: Sequence[Path]) -> Tuple[List[Tuple[Path, Path]], List[Path]]:
    pairs: List[Tuple[Path, Path]] = []
    skipped: List[Path] = []
    for r1 in r1_files:
        r2 = _r1_to_r2_path(r1)
        if r2 is None or not r2.is_file():
            skipped.append(r1)
            continue
        pairs.append((r1, r2))
    # stable order by first file of each pair
    pairs.sort(key=lambda t: t[0].name.lower())
    return pairs, skipped


def _find_sampleinfo(
    version_dir: Path,
    experiment_dir: Path,
    my_data_root: Optional[Path] = None,
    experiment_name: Optional[str] = None,
) -> Optional[Path]:
    """Prefer EXP_sampleinfo.xlsx, then *sampleinfo*.xlsx in version then experiment.

    Also checks ``<MY_DATA>/<Experiment>_sampleinfo.xlsx`` (e.g. MY_DATA/M265_sampleinfo.xlsx
    next to folder MY_DATA/M265/).
    """
    candidates: List[Path] = []

    def collect(base: Path) -> None:
        if not base.is_dir():
            return
        for p in base.iterdir():
            if not p.is_file():
                continue
            low = p.name.lower()
            if not low.endswith((".xlsx", ".xls")):
                continue
            if "sampleinfo" in low or low.startswith("exp_sampleinfo"):
                candidates.append(p)

    collect(version_dir)
    if not candidates:
        collect(experiment_dir)

    if not candidates and my_data_root is not None and experiment_name:
        root = my_data_root.resolve()
        if root.is_dir():
            sibling = root / f"{experiment_name}_sampleinfo.xlsx"
            if sibling.is_file():
                return sibling
            for pat in (
                f"{experiment_name}_sampleinfo.xls",
                f"{experiment_name}*sampleinfo*.xlsx",
                f"{experiment_name}*sampleinfo*.xls",
            ):
                for p in sorted(root.glob(pat)):
                    if p.is_file() and "sampleinfo" in p.name.lower():
                        return p

    if not candidates:
        return None

    def key(p: Path) -> Tuple[int, str]:
        n = p.name.lower()
        if n == "exp_sampleinfo.xlsx":
            return (0, n)
        if n.startswith("exp_") and "sampleinfo" in n:
            return (1, n)
        if "sampleinfo" in n:
            return (2, n)
        return (3, n)

    candidates.sort(key=key)
    return candidates[0]


def _find_config(
    experiment_name: str,
    version_dir: Path,
    experiment_dir: Path,
    fallback: Path,
    my_data_root: Optional[Path] = None,
) -> Path:
    for base in (version_dir, experiment_dir):
        exact = base / f"{experiment_name}.mapseq.conf"
        if exact.is_file():
            return exact
        matches = sorted(base.glob("*.mapseq.conf"))
        if matches:
            return matches[0]
    if my_data_root is not None and my_data_root.is_dir():
        sib = my_data_root / f"{experiment_name}.mapseq.conf"
        if sib.is_file():
            return sib
        for p in sorted(my_data_root.glob(f"{experiment_name}*.mapseq.conf")):
            if p.is_file():
                return p
    if fallback.is_file():
        return fallback
    return fallback


def _discover_versions(experiment_dir: Path) -> List[Tuple[str, Path]]:
    """
    Return list of (version_label, version_path).
    - If R1 FASTQs live directly under experiment_dir, one version '__root__'.
    - Each immediate subdirectory that contains R1 FASTQs is its own version.
    """
    exp = experiment_dir.resolve()
    name = exp.name
    versions: List[Tuple[str, Path]] = []

    root_r1 = _find_fastq_r1_files(exp)
    subs = sorted([p for p in exp.iterdir() if p.is_dir()], key=lambda p: p.name.lower())
    sub_with_fastq = [s for s in subs if _find_fastq_r1_files(s)]

    if root_r1:
        versions.append(("__root__", exp))
    for s in sub_with_fastq:
        # Avoid duplicating exp if __root__ already used same path
        if s.resolve() == exp.resolve():
            continue
        versions.append((s.name, s))

    if not versions and not root_r1 and not sub_with_fastq:
        return []

    if not versions and sub_with_fastq:
        for s in sub_with_fastq:
            versions.append((s.name, s))

    # Deduplicate by path
    seen = set()
    uniq: List[Tuple[str, Path]] = []
    for label, p in versions:
        rp = p.resolve()
        if rp in seen:
            continue
        seen.add(rp)
        uniq.append((label, p))
    return uniq


def _find_read_count_pdfs(root: Path) -> List[Path]:
    out: List[Path] = []
    if not root.is_dir():
        return out
    for dirpath, _, filenames in os.walk(root):
        for f in filenames:
            if not f.endswith(".pdf"):
                continue
            if "read_count" not in f:
                continue
            out.append(Path(dirpath) / f)
    return sorted(out)


def _dest_pdf_dir(
    pdf_dest_root: Path,
    experiment_name: str,
    version_label: str,
    pair_slug: str,
) -> Path:
    safe_ver = version_label.replace(os.sep, "_").replace("..", "_")
    safe_pair = pair_slug.replace(os.sep, "_").replace("..", "_")
    return pdf_dest_root / experiment_name / safe_ver / safe_pair


def _append_manifest(
    manifest_path: Path,
    rows: Iterable[dict],
    write_header: bool,
) -> None:
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "experiment",
        "version",
        "pair",
        "action",
        "detail",
        "source_pdf",
        "dest_pdf",
    ]
    mode = "w" if write_header else "a"
    with manifest_path.open(mode, newline="", encoding="utf-8") as fp:
        w = csv.DictWriter(fp, fieldnames=fieldnames, delimiter="\t")
        if write_header:
            w.writeheader()
        for row in rows:
            w.writerow({k: row.get(k, "") for k in fieldnames})


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Run MAPseq pipeline per FASTQ pair (explicit stage scripts) and collect read_count PDFs."
    )
    parser.add_argument(
        "--my-data-root",
        type=Path,
        default=_REPO_ROOT / "MY_DATA",
        help="Root folder containing experiment directories (default: <repo>/MY_DATA)",
    )
    parser.add_argument(
        "--pdf-dest-root",
        type=Path,
        default=None,
        help="Under this path, PDFs go to {pdf-dest-root}/{experiment}/{version}/{pair}/ "
        "(default: <my-data-root>/pdfs)",
    )
    parser.add_argument(
        "--pipeline-out-subdir",
        type=str,
        default="mapseq_pipeline",
        help="Under each version folder, pair outputs go in <this>/<pair_slug>/",
    )
    parser.add_argument(
        "--config",
        type=Path,
        default=_DEFAULT_CONFIG,
        help="Fallback mapseq.conf if none found next to sampleinfo",
    )
    parser.add_argument(
        "--execute",
        action="store_true",
        help="Actually run pipeline stage scripts (default: dry-run only)",
    )
    parser.add_argument(
        "--halt",
        metavar="STAGE",
        choices=STEPLIST,
        default=None,
        help="Stop after this stage (inclusive), e.g. readtable for PDF-only. Default: full pipeline through matrices.",
    )
    parser.add_argument(
        "--dask-temp",
        type=Path,
        default=None,
        help="Dask temp dir for aggregate_reads; also passes -k (Dask) to lower RAM on that step.",
    )
    parser.add_argument(
        "--no-skip-if-dest-has-pdf",
        action="store_false",
        dest="skip_if_dest_has_pdf",
        default=True,
        help="By default, skip pipeline if destination already has a *read_count*.pdf; "
        "pass this to always run when executing",
    )
    parser.add_argument(
        "--max-failures",
        type=int,
        default=0,
        help="Stop after this many subprocess failures (0 = no limit)",
    )
    parser.add_argument(
        "--only-experiments",
        nargs="+",
        metavar="NAME",
        default=None,
        help="Only process these experiment folder names (case-insensitive), e.g. M265 M275",
    )
    parser.add_argument("-v", "--verbose", action="store_true")
    args = parser.parse_args()

    _setup_logging(args.verbose)

    my_data: Path = args.my_data_root.resolve()
    pdf_dest_root = (
        args.pdf_dest_root.resolve()
        if args.pdf_dest_root
        else (my_data / "pdfs").resolve()
    )
    manifest_path = pdf_dest_root / "pdf_manifest.tsv"

    if not my_data.is_dir():
        logging.error("MY_DATA root does not exist: %s", my_data)
        return 2

    if not _PROCESS_FASTQ.is_file():
        logging.error("process_fastq_pairs.py not found: %s", _PROCESS_FASTQ)
        return 2

    experiments = sorted(
        [p for p in my_data.iterdir() if p.is_dir() and p.name != "pdfs"],
        key=lambda p: p.name.lower(),
    )
    only_set = None
    if args.only_experiments:
        only_set = {n.lower() for n in args.only_experiments}

    failures = 0
    manifest_needs_header = not manifest_path.exists()

    for exp_dir in experiments:
        exp_name = exp_dir.name
        if only_set is not None and exp_name.lower() not in only_set:
            continue
        versions = _discover_versions(exp_dir)
        if not versions:
            logging.warning("[%s] No FASTQ versions found; skipping", exp_name)
            _append_manifest(
                manifest_path,
                [
                    {
                        "experiment": exp_name,
                        "version": "",
                        "pair": "",
                        "action": "skip",
                        "detail": "no_fastq_versions",
                        "source_pdf": "",
                        "dest_pdf": "",
                    }
                ],
                write_header=manifest_needs_header,
            )
            manifest_needs_header = False
            continue

        for version_label, version_path in versions:
            sampleinfo = _find_sampleinfo(
                version_path, exp_dir, my_data_root=my_data, experiment_name=exp_name
            )
            config = _find_config(
                exp_name,
                version_path,
                exp_dir,
                args.config.resolve(),
                my_data_root=my_data,
            )

            r1_list = _find_fastq_r1_files(version_path)
            pairs, missing_r2 = _pair_fastqs(r1_list)

            if not pairs:
                logging.warning(
                    "[%s / %s] No complete R1/R2 pairs (%d R1 missing R2)",
                    exp_name,
                    version_label,
                    len(missing_r2),
                )
                _append_manifest(
                    manifest_path,
                    [
                        {
                            "experiment": exp_name,
                            "version": version_label,
                            "pair": "",
                            "action": "skip",
                            "detail": "no_complete_pairs",
                            "source_pdf": "",
                            "dest_pdf": "",
                        }
                    ],
                    write_header=manifest_needs_header,
                )
                manifest_needs_header = False
                continue

            if sampleinfo is None:
                logging.warning(
                    "[%s / %s] No sampleinfo .xlsx; cannot run pipeline. "
                    "Copying any existing PDFs from per-pair output dirs only.",
                    exp_name,
                    version_label,
                )

            pair_base = version_path / args.pipeline_out_subdir
            dask_temp = args.dask_temp.resolve() if args.dask_temp else None

            for r1, r2 in pairs:
                pair_slug = _pair_slug(r1)
                pair_outdir = pair_base / pair_slug
                dest_dir = _dest_pdf_dir(
                    pdf_dest_root, exp_name, version_label, pair_slug
                )
                existing_dest = _find_read_count_pdfs(dest_dir)
                rows_for_manifest: List[dict] = []

                do_print_commands = sampleinfo is not None and not args.execute
                do_run_subprocess = (
                    args.execute
                    and sampleinfo is not None
                    and not (
                        args.skip_if_dest_has_pdf
                        and bool(existing_dest)
                    )
                )

                rc = 0
                if sampleinfo is not None and (do_print_commands or do_run_subprocess):
                    rc = _run_pipeline_for_pair(
                        config=config,
                        sampleinfo=sampleinfo,
                        pair_outdir=pair_outdir,
                        r1=r1,
                        r2=r2,
                        halt=args.halt,
                        dask_temp=dask_temp,
                        execute=do_run_subprocess,
                    )
                elif (
                    args.execute
                    and sampleinfo is not None
                    and args.skip_if_dest_has_pdf
                    and existing_dest
                ):
                    logging.info(
                        "[%s / %s / %s] Destination already has read_count PDFs; skipping run",
                        exp_name,
                        version_label,
                        pair_slug,
                    )

                if rc != 0:
                    failures += 1
                    logging.error(
                        "[%s / %s / %s] pipeline failed with code %s",
                        exp_name,
                        version_label,
                        pair_slug,
                        rc,
                    )
                    rows_for_manifest.append(
                        {
                            "experiment": exp_name,
                            "version": version_label,
                            "pair": pair_slug,
                            "action": "pipeline_failed",
                            "detail": str(rc),
                            "source_pdf": "",
                            "dest_pdf": "",
                        }
                    )
                    _append_manifest(
                        manifest_path,
                        rows_for_manifest,
                        write_header=manifest_needs_header,
                    )
                    manifest_needs_header = False
                    if args.max_failures and failures >= args.max_failures:
                        logging.error(
                            "Max failures (%s) reached; stopping.", args.max_failures
                        )
                        return 1
                    continue

                pdfs = _find_read_count_pdfs(pair_outdir)
                dest_dir.mkdir(parents=True, exist_ok=True)
                for pdf in pdfs:
                    target = dest_dir / pdf.name
                    if args.execute:
                        shutil.copy2(pdf, target)
                        rows_for_manifest.append(
                            {
                                "experiment": exp_name,
                                "version": version_label,
                                "pair": pair_slug,
                                "action": "copied",
                                "detail": "",
                                "source_pdf": str(pdf),
                                "dest_pdf": str(target),
                            }
                        )
                        logging.info("Copied %s -> %s", pdf, target)
                    else:
                        rows_for_manifest.append(
                            {
                                "experiment": exp_name,
                                "version": version_label,
                                "pair": pair_slug,
                                "action": "dry_run_copy",
                                "detail": "",
                                "source_pdf": str(pdf),
                                "dest_pdf": str(target),
                            }
                        )
                        logging.info(
                            "[dry-run] would copy %s -> %s", pdf, target
                        )

                if not pdfs:
                    rows_for_manifest.append(
                        {
                            "experiment": exp_name,
                            "version": version_label,
                            "pair": pair_slug,
                            "action": "no_pdfs",
                            "detail": "no read_count PDFs under %s" % pair_outdir,
                            "source_pdf": "",
                            "dest_pdf": str(dest_dir),
                        }
                    )

                _append_manifest(
                    manifest_path,
                    rows_for_manifest,
                    write_header=manifest_needs_header,
                )
                manifest_needs_header = False

                if args.max_failures and failures >= args.max_failures:
                    logging.error(
                        "Max failures (%s) reached; stopping.", args.max_failures
                    )
                    return 1

    logging.info("Done. Manifest: %s", manifest_path)
    return 0 if failures == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
