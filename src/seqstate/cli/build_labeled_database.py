from __future__ import annotations

import argparse
from pathlib import Path

import polars
import rich
from rich.console import Console

from seqstate.cli._parser import configure_general_arguments
from seqstate.cli.utils import create_versioned_path
from seqstate.database import LabeledGenomeDatabase
from seqstate.states import FeatureTable, GeneTable
from seqstate.utils import is_parquet_file
from seqstate.vocab import Labels, Vocab

"""
export PYTHONPATH='/home/ray/Documents/seqstate/src/'
GENES='/home/ray/Documents/seqstate/test/test_data/genes_1000.train.genes.tsv'
DOMAINS='/home/ray/Documents/seqstate/test/test_data/pfams_matching.train.features.tsv'
FEATURES="/home/ray/Documents/seqstate/test/test_data/kofams_matching.train.features.tsv"
VOCAB='/home/ray/Documents/seqstate/test/test_data/genes_1000.train.vocab.txt'
LABELS='kofam_labels.txt'

GENES='/home/ray/Documents/seqstate/test/test_data/genes_1000.test.genes.tsv'
DOMAINS='/home/ray/Documents/seqstate/test/test_data/pfams_matching.test.features.tsv'
FEATURES='/home/ray/Documents/seqstate/test/test_data/kofams_matching.test.features.tsv'

python -m seqstate.cli build-labeled-database -g $GENES --pfams $DOMAINS --labels $FEATURES --pfam-vocab $VOCAB --label-vocab $LABELS -o database_test.sqlite3
"""


def configure_build_labeled_database_input(
    parser: argparse.ArgumentParser,
) -> argparse._ArgumentGroup:

    group = parser.add_argument_group(
        "Build Database",
        "Arguments to build a labeled genome SQL database",
    )
    group.add_argument(
        "-g",
        "--genes",
        required=True,
        help="Path to the genes table as .tsv or .parquet file",
        type=Path,
    )
    group.add_argument(
        "--pfams",
        required=True,
        help="Path to the Pfam features table as .tsv or .parquet file",
        type=Path,
    )
    group.add_argument(
        "--labels",
        required=True,
        help="Path to the label feature table as .tsv or .parquet file",
        type=Path,
    )
    group.add_argument(
        "--label-col",
        required=True,
        help="Column name of the labels in the label feature table",
        type=str,
    )
    group.add_argument(
        "-o",
        "--output",
        required=True,
        help="Path to the output database to write",
        type=Path,
    )
    group.add_argument(
        "--pfam-vocab",
        required=True,
        help="Path to the input vocab to encode domains with",
        type=Path,
    )
    group.add_argument(
        "--label-vocab",
        required=True,
        help="Path to the label vocab to encode labels with",
        type=Path,
    )
    group.add_argument(
        "--n-samples",
        required=False,
        help="How many labeled annotated genes to sample.",
        type=int,
    )

    return group


def configure_parser(parser: argparse.ArgumentParser):
    """Configure the main parser with argument groups for different namespaces"""

    configure_general_arguments(parser)

    configure_build_labeled_database_input(parser)

    parser.set_defaults(run=run)


def run(args: argparse.Namespace, console: Console) -> int:
    """Run function for the CLI."""

    # abort early if output database file alread exits
    if args.output.exists():
        rich.print(f"[bold red]Error[/]: File {str(args.output)!r} already exists")
        return 1

    if is_parquet_file(args.labels):
        labels_lf = polars.scan_parquet(args.labels)
    else:
        labels_lf = polars.scan_csv(args.labels, separator="\t")

    # TODO this should have something like a n features per gene thing?
    labeled_genes = GeneTable.load(args.genes).label(
        labels_lf, label_name=args.label_col
    )

    pfam_features = FeatureTable.load(args.pfams)

    rich.print("Loading Vocab and Labels")
    with args.pfam_vocab.open() as src:
        pfam_vocab = Vocab.load(src)

    with args.label_vocab.open() as src:
        label_vocab = Labels.load(src)

    oo = create_versioned_path(args.output)

    labeled_db = LabeledGenomeDatabase.build_labeled_database(
        labeled_genes=labeled_genes,
        pfam_features=pfam_features,
        pfam_vocab=pfam_vocab,
        label_vocab=label_vocab,
        output=oo,
        console=console,
        seed=args.seed,
    )

    labeled_db.close()

    console.print(
        f"[bold green]{'Finished':>12}[/] building database to {str(oo.resolve())}."
    )

    return 0
