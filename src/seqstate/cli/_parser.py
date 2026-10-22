from __future__ import annotations

import argparse
import os
from pathlib import Path

from seqstate import __package__ as _module
from seqstate import __version__


def configure_general_arguments(
    parser: argparse.ArgumentParser,
) -> argparse._ArgumentGroup:

    group = parser.add_argument_group()

    group.add_argument(
        "-V",
        "--version",
        action="version",
        version=f"{_module} {__version__}",
        help="Show the program version number and exit.",
    )
    group.add_argument(
        "-j",
        "--jobs",
        default=max(0, len(os.sched_getaffinity(0))),  # works with slurms cpus-per-task
        type=int,
        help="The number of jobs to run in parallel sections.",
    )
    group.add_argument(
        "--seed",
        default=42,
        type=int,
        help="The seed to use for initializing pseudo-random number generators.",
    )

    return group


def configure_build_vocab_input(
    parser: argparse.ArgumentParser,
) -> argparse._ArgumentGroup:

    group = parser.add_argument_group(
        "Build Vocab",
        "Arguments to build a Vocab of Pfams from a gene and a feature table",
    )
    group.add_argument(
        "-g",
        "--genes",
        required=False,
        help="If provided, only a Vocab for the features of the supplied genes will be generated",
        type=Path,
    )
    group.add_argument(
        "-f",
        "--features",
        required=True,
        help="Path to the Pfam features table as .tsv or .parquet file",
        type=Path,
    )
    group.add_argument(
        "--domain-number",
        required=False,
        type=int,
        help="Number of feature domains per gene to generate Vocab from. Default all",
    )
    group.add_argument(
        "--unique-domains",
        action="store_true",
        help="Consider only unique domains per Gene for building the Vocab. (Affects db creation too!)",
    )
    group.add_argument(
        "--min-domain-occurence",
        default=0,
        type=int,
        help="Minimum occurence of a pfam domain to be included in Vocab. Default no minimum.",
    )
    group.add_argument(
        "-o",
        "--output",
        required=True,
        help="Path to the output vocab to write",
        type=Path,
    )
    group.add_argument(
        "--genome-database",
        required=False,
        help="Optional. Also requires genes. If passed will write a GenomeDatabase too!",
        type=Path,
    )

    return group
