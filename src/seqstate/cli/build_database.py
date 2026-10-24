from __future__ import annotations

import argparse
from pathlib import Path

from rich.console import Console

from seqstate.cli._parser import configure_general_arguments
from seqstate.cli.utils import create_versioned_path
from seqstate.database import GenomeDatabase
from seqstate.tables import FeatureTable, GeneTable
from seqstate.vocab import Vocab

"""
export PYTHONPATH='/home/ray/Documents/seqstate/src/'
GENES='/home/ray/Documents/seqstate/test/test_data/genes_1000.train.genes.tsv'
FEATURES='/home/ray/Documents/seqstate/test/test_data/genes_1000.train.pfam.features.tsv'
VOCAB='/home/ray/Documents/seqstate/test/test_data/genes_1000.vocab.txt'

GENES='/home/ray/Documents/seqstate/test/test_data/genes_1000.test.genes.tsv'
FEATURES='/home/ray/Documents/seqstate/test/test_data/genes_1000.test.pfam.features.tsv'

python -m seqstate.cli build-database -g $GENES -f $FEATURES --vocab $VOCAB -o database_test.sqlite3 --unique-domains
"""


def configure_build_database_input(
    parser: argparse.ArgumentParser,
) -> argparse._ArgumentGroup:

    group = parser.add_argument_group(
        "Build Database",
        "Arguments to build a genome SQL database",
    )
    group.add_argument(
        "-f",
        "--features",
        required=True,
        help="Path to the Pfam features table as .tsv or .parquet file",
        type=Path,
    )
    group.add_argument(
        "-g",
        "--genes",
        required=True,
        help="Path to the genes table as .tsv or .parquet file",
        type=Path,
    )
    group.add_argument(
        "-o",
        "--output",
        required=True,
        help="Path to the output database to write",
        type=Path,
    )
    group.add_argument(
        "--vocab",
        required=True,
        help="Path to the input vocab to encode domains with",
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

    return group


def configure_parser(parser: argparse.ArgumentParser):
    """Configure the main parser with argument groups for different namespaces"""

    configure_general_arguments(parser)

    configure_build_database_input(parser)

    parser.set_defaults(run=run)


def run(args: argparse.Namespace, console: Console) -> int:
    """Run function for the CLI."""

    genes = GeneTable.load(args.genes)

    features = FeatureTable.load(args.features)

    console.print("Loading Vocab")
    with args.vocab.open() as src:
        vocab = Vocab.load(src)

    oo = create_versioned_path(args.output)

    # generate the database file
    db = GenomeDatabase.build_database(
        genes=genes,
        features=features,
        vocab=vocab,
        output=oo,
        console=console,
        domain_number=args.domain_number,
        only_unique_domains_per_gene=args.unique_domains,
        index_db=True,
    )
    db.close()

    console.print(
        f"[bold green]{'Finished':>12}[/] building database to {str(oo.resolve())}."
    )

    return 0
