"""
Build a Vocab from a ``features.tsv`` table.
If a ``genes.tsv`` table is supplied,
only features which appear in the gene table are saved to the Vocab.
"""

from __future__ import annotations

from argparse import ArgumentParser, Namespace
from pathlib import Path

import polars
from rich.console import Console

from seqstate.cli._parser import (configure_build_vocab_input,
                                  configure_general_arguments)
from seqstate.cli.utils import create_versioned_path
from seqstate.states import GeneTable
from seqstate.utils import is_parquet_file
from seqstate.vocab import Labels

"""
export PYTHONPATH='/home/ray/Documents/seqstate/src/'
GENES='/home/ray/Documents/seqstate/test/test_data/genes_1000.train.genes.tsv'
FEATURES='/home/ray/Documents/seqstate/test/test_data/kofams_matching.train.features.tsv'

python -m seqstate.cli build-labels -g $GENES -f $FEATURES  -o labels_test.txt --domain-number 1 --min-domain-occurence 1
"""


def configure_parser(parser: ArgumentParser):
    """Configure the main parser with argument groups for different namespaces"""

    configure_general_arguments(parser)

    configure_build_vocab_input(parser)

    parser.add_argument(
        "--label-col",
        required=True,
        type=str,
        help="Label column name of the label in the feature table.",
    )

    parser.set_defaults(run=run)


def run(args: Namespace, console: Console) -> int:
    """Run function for the CLI."""

    # NOTE
    # this fn is currently specific for building KOfam labels
    # (which require additonal filtering!)
    # not all labels will have i_e values for example ...

    # TODO fix this to work with stuff that doesnt have i_evalues....

    # abort early if output database file alread exits
    if args.output.exists():
        console.print(f"[bold red]Error[/]: File {str(args.output)!r} already exists")
        return 1

    genes = GeneTable.load(args.genes)

    if is_parquet_file(args.features):
        features = polars.scan_parquet(args.features).select(
            ["protein_id", args.label_col, "i_evalue"]
        )
    else:
        features = polars.scan_csv(args.features, separator="\t").select(
            ["protein_id", args.label_col, "i_evalue"]
        )

    console.print("Loading Features and annotating Genes")
    console.print("Filtering out features below an evalue of 1e-9.")
    console.print(
        f"Filtering out features occuring less than {args.min_domain_occurence} times"
    )

    # generate the labels
    labels = Labels.build_labels(
        features=features,
        genes=genes.table,
        min_label_occurence=args.min_domain_occurence,
        label_col=args.label_col,
    )

    oo = create_versioned_path(args.output)
    with Path(oo).open("w") as dst:
        labels.dump(dst)

    console.print(
        f"[bold green]{'Finished':>12}[/] building labels and writing to {str(oo.resolve())}"
    )

    return 0
