"""
Build a Vocab from a ``features.tsv`` table.
If a ``genes.tsv`` table is supplied,
only features which appear in the gene table are saved to the Vocab.
"""

from __future__ import annotations

from argparse import ArgumentParser, Namespace
from pathlib import Path

from rich.console import Console

from seqstate.cli._parser import (configure_build_vocab_input,
                                  configure_general_arguments)
from seqstate.cli.utils import create_versioned_path
from seqstate.database import GenomeDatabase
from seqstate.states import FeatureTable, GeneTable

"""
export PYTHONPATH='/home/ray/Documents/seqstate/src/'
GENES='/home/ray/Documents/seqstate/test/test_data/genes_1000.train.genes.tsv'
FEATURES='/home/ray/Documents/seqstate/test/test_data/genes_1000.train.pfam.features.tsv'

python -m seqstate.cli build-vocab -g $GENES -f $FEATURES  -o vocab_test.txt --unique-domains --genome-database vocab_db.sqlite3
"""


def configure_parser(parser: ArgumentParser):
    """Configure the main parser with argument groups for different namespaces"""

    configure_general_arguments(parser)

    configure_build_vocab_input(parser)

    parser.set_defaults(run=run)


def run(args: Namespace, console: Console) -> int:
    """Run function for the CLI."""

    # abort early if output database file alread exits
    if args.output.exists():
        console.print(f"[bold red]Error[/]: File {str(args.output)!r} already exists")
        return 1
    if args.genes is None:
        raise ValueError(
            "Not passing genes is disabled for safety reasons. Please conciously pass the training genes!"
        )

    console.print("Loading Features and annotating Genes")
    genes = GeneTable.load(args.genes)
    features = FeatureTable.load(args.features)

    # generate the vocab
    vocab, annotated_features = features.build_vocab(
        genes=genes,
        domain_number=args.domain_number,
        min_domain_occurence=args.min_domain_occurence,
        only_unique_domains_per_gene=args.unique_domains,
    )

    oo = create_versioned_path(args.output)
    console.print(
        f"[bold green]{'Finished':>12}[/] building vocab and writing to {args.output}"
    )
    with Path(oo).open("w") as dst:
        vocab.dump(dst)

    db_oo = create_versioned_path(args.genome_database)

    if args.genome_database is not None and genes is not None:
        db = GenomeDatabase._build_sql_db(
            genes_df=genes.table.collect(),
            features_df=annotated_features,
            vocab=vocab,
            console=console,
            output=db_oo,
            index_db=True,
        )
        db.close()
        console.print(
            f"[bold green]{'Finished':>12}[/] building database to {str(db_oo.resolve())}."
        )

    return 0
