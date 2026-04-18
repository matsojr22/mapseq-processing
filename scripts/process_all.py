#!/usr/bin/env python
#
# Top level processing script.
# Runs all MAPseq processing steps from FASTQ pairs.
#
import argparse
import datetime
import logging
import os
import subprocess
import sys

from configparser import ConfigParser

# Ensure repo root (one level above scripts/) is importable.
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.dirname(SCRIPT_DIR)
if REPO_ROOT not in sys.path:
    sys.path.append(REPO_ROOT)


STEPLIST = [
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


def _validate_inputs(config_file, sampleinfo_file, infiles):
    """Validate top-level files and paired FASTQ argv."""
    if not os.path.isfile(config_file):
        raise FileNotFoundError(f"config file does not exist: {config_file}")

    if not os.path.isfile(sampleinfo_file):
        raise FileNotFoundError(f"sampleinfo file does not exist: {sampleinfo_file}")

    if not infiles or len(infiles) < 2:
        raise ValueError("at least one R1/R2 pair is required (2+ FASTQ files).")

    if len(infiles) % 2 != 0:
        raise ValueError("FASTQ inputs must be even: R1 R2 [R1b R2b ...].")

    for fn in infiles:
        if not os.path.isfile(fn):
            raise FileNotFoundError(f"FASTQ input does not exist: {fn}")


def _steps_until(halt):
    """Return the effective step list without mutating globals."""
    if halt is None:
        return list(STEPLIST)
    if halt not in STEPLIST:
        raise ValueError(f"invalid halt step {halt!r}. Must be one of {STEPLIST}")
    return STEPLIST[: STEPLIST.index(halt) + 1]


def _step_output_exists(step_name, soutfile, soutdir):
    """
    Determine if a step output already exists.
    For matrices, treat non-empty output directory as done.
    """
    if step_name == "matrices":
        if soutdir is None or not os.path.isdir(soutdir):
            return False
        return len(os.listdir(soutdir)) > 0
    return (soutfile is not None) and os.path.exists(soutfile)


def process_mapseq_all(
    config_file,
    sampleinfo_file,
    infiles,
    outdir=None,
    force=False,
    halt=None,
    aggregate_use_dask=False,
    aggregate_dask_temp=None,
    stage_loglevel_flag="-v",
):
    """
    Performs end-to-end default processing.
    Executes each pipeline script in an external process.
    """
    logging.info(
        "%s input files. config=%s sampleinfo=%s outdir=%s force=%s halt=%s "
        "aggregate_use_dask=%s",
        len(infiles),
        config_file,
        sampleinfo_file,
        outdir,
        force,
        halt,
        aggregate_use_dask,
    )

    _validate_inputs(
        config_file=config_file,
        sampleinfo_file=sampleinfo_file,
        infiles=infiles,
    )

    cp = ConfigParser()
    loaded = cp.read(config_file)
    if len(loaded) < 1:
        raise ValueError(f"No valid configuration loaded from {config_file}")
    if not cp.has_option("project", "project_id"):
        raise ValueError(f"Config missing required [project] project_id: {config_file}")
    project_id = cp.get("project", "project_id")

    if outdir is None:
        outdir = os.path.abspath("./")
    else:
        outdir = os.path.abspath(outdir)
    os.makedirs(outdir, exist_ok=True)

    steps = _steps_until(halt)
    logging.debug("effective STEPLIST=%s", steps)

    for idx, step in enumerate(steps):
        runstep = True
        soutdir = None
        soutfile = None
        sprog = STEPMAP[step]
        sname = DIRMAP[sprog]
        logging.debug("handling step=%s sprog=%s sname=%s", step, sprog, sname)

        if sname == "matrices":
            soutdir = os.path.join(outdir, f"{sname}.out")
        else:
            soutfile = os.path.join(outdir, f"{sname}.out/{project_id}.{sname}.tsv")

        infile = None
        if sname != "reads":
            instep = steps[idx - 1]
            inprog = STEPMAP[instep]
            insname = DIRMAP[inprog]
            infile = os.path.join(outdir, f"{insname}.out/{project_id}.{insname}.tsv")

        log_file = os.path.join(outdir, f"{step}.log")
        cmd = [
            sys.executable,
            os.path.join(SCRIPT_DIR, f"{sprog}.py"),
        ]
        if stage_loglevel_flag in ("-d", "-v"):
            cmd.append(stage_loglevel_flag)
        cmd.extend(
            [
                "-c",
                os.path.abspath(config_file),
                "-L",
                log_file,
            ]
        )

        if soutfile is not None:
            cmd.extend(["-o", soutfile])

        if soutdir is not None:
            cmd.extend(["-O", soutdir])

        if sname in ("readtable", "vbctable"):
            cmd.extend(["-s", os.path.abspath(sampleinfo_file)])

        if sname == "reads":
            cmd.extend(["-s", os.path.abspath(sampleinfo_file)])
            for fn in infiles:
                cmd.append(os.path.abspath(fn))
        else:
            cmd.append(infile)
            if not os.path.exists(infile):
                raise FileNotFoundError(f"required infile missing for step {step}: {infile}")

        if sname == "aggregated" and aggregate_use_dask:
            cmd.append("-k")
            if aggregate_dask_temp is not None:
                cmd.extend(["-t", os.path.abspath(aggregate_dask_temp)])

        logging.debug("made command=%s", " ".join(cmd))

        if (not force) and _step_output_exists(sname, soutfile, soutdir):
            runstep = False
            logging.info("output exists for step=%s; skipping", step)

        if runstep:
            logging.info("running step=%s (log=%s)", step, log_file)
            start = datetime.datetime.now()
            cp = subprocess.run(cmd, cwd=REPO_ROOT)
            end = datetime.datetime.now()
            elapsed = end - start
            logging.info(
                "finished step=%s rc=%s in %s seconds",
                step,
                cp.returncode,
                elapsed.seconds,
            )
            if cp.returncode != 0:
                raise RuntimeError(
                    f"step '{step}' failed with rc={cp.returncode}: {' '.join(cmd)}"
                )


if __name__ == "__main__":
    FORMAT = (
        "%(asctime)s (UTC) [ %(levelname)s ] %(filename)s:%(lineno)d "
        "%(name)s.%(funcName)s(): %(message)s"
    )
    logging.basicConfig(format=FORMAT)
    logging.getLogger().setLevel(logging.WARN)

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "-d",
        "--debug",
        action="store_true",
        dest="debug",
        help="debug logging",
    )

    parser.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        dest="verbose",
        help="verbose logging",
    )

    parser.add_argument(
        "-c",
        "--config",
        metavar="config",
        required=True,
        default=os.path.expanduser("~/git/mapseq-processing/etc/mapseq.conf"),
        type=str,
        help="config file.",
    )

    parser.add_argument(
        "-s",
        "--sampleinfo",
        metavar="sampleinfo",
        required=True,
        default=None,
        type=str,
        help="XLS sampleinfo file.",
    )

    parser.add_argument(
        "-O",
        "--outdir",
        metavar="outdir",
        required=False,
        default=None,
        type=str,
        help="outdir. output file base dir if not given.",
    )

    parser.add_argument(
        "-f",
        "--force",
        action="store_true",
        default=False,
        help="Recalculate even if output exists.",
    )

    parser.add_argument(
        "-H",
        "--halt",
        metavar="halt",
        required=False,
        default=None,
        choices=STEPLIST,
        type=str,
        help=f"Stage name to stop after: {STEPLIST}",
    )

    parser.add_argument(
        "-k",
        "--aggregate_use_dask",
        action="store_true",
        default=False,
        help="Pass -k to aggregate_reads.py (Dask for aggregation step).",
    )

    parser.add_argument(
        "-t",
        "--aggregate_dask_temp",
        metavar="aggregate_dask_temp",
        required=False,
        default=None,
        type=str,
        help="Optional temp path passed as -t to aggregate_reads.py when -k is set.",
    )

    parser.add_argument(
        "-L",
        "--logfile",
        metavar="logfile",
        required=False,
        default=None,
        type=str,
        help="Logfile for this wrapper process.",
    )

    parser.add_argument(
        "infiles",
        metavar="infiles",
        type=str,
        nargs="*",
        default=None,
        help="Initial FASTQ input to process (R1 R2 [R1b R2b ...]).",
    )

    args = parser.parse_args()

    if args.debug:
        logging.getLogger().setLevel(logging.DEBUG)
    elif args.verbose:
        logging.getLogger().setLevel(logging.INFO)

    logging.debug("indirs=%s", args.infiles)

    cp = ConfigParser()
    cp.read(args.config)
    cdict = {
        section: dict(cp.items(section))
        for section in cp.sections()
    }
    logging.debug("Running with config. %s: %s", args.config, cdict)
    logging.debug("infiles=%s", args.infiles)

    if args.logfile is not None:
        log = logging.getLogger()
        formatter = logging.Formatter(FORMAT)
        log_stream = logging.FileHandler(filename=args.logfile)
        log_stream.setFormatter(formatter)
        log.addHandler(log_stream)

    outdir = os.path.abspath("./")
    if args.outdir is not None:
        outdir = os.path.abspath(args.outdir)
    os.makedirs(outdir, exist_ok=True)

    if args.aggregate_dask_temp is not None and not args.aggregate_use_dask:
        logging.warning(
            "--aggregate_dask_temp provided without --aggregate_use_dask; ignoring temp path."
        )

    stage_loglevel_flag = "-v"
    if args.debug:
        stage_loglevel_flag = "-d"
    elif args.verbose:
        stage_loglevel_flag = "-v"

    try:
        process_mapseq_all(
            config_file=args.config,
            sampleinfo_file=args.sampleinfo,
            infiles=args.infiles,
            outdir=outdir,
            force=args.force,
            halt=args.halt,
            aggregate_use_dask=args.aggregate_use_dask,
            aggregate_dask_temp=args.aggregate_dask_temp,
            stage_loglevel_flag=stage_loglevel_flag,
        )
    except Exception as exc:
        logging.error("process_all failed: %s", exc)
        sys.exit(1)

    logging.info("Done process_all.")
   