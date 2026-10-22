import argparse
from pathlib import Path

from rich.console import Console

from seqstate.cli._parser import configure_general_arguments
from seqstate.cli.utils import create_versioned_path
from seqstate.database import GenomeDatabase

"""
export PYTHONPATH='/home/ray/Documents/seqstate/src/'
DB1=/home/ray/Documents/seqstate/test/test_data/genes_1000.test.sqlite3
DB2=/home/ray/Documents/seqstate/test/test_data/genes_1000.train.sqlite3
python -m seqstate.cli merge-dbs -i $DB1 $DB2 -o out
"""


def configure_merging_databases(
    parser: argparse.ArgumentParser,
) -> argparse._ArgumentGroup:

    group = parser.add_argument_group(
        title="Merge Database Arguments",
        description="Merge a list of genome SQL databases into a new larger one.",
    )

    group.add_argument(
        "-i",
        "--input",
        nargs="+",
        type=Path,
        required=True,
        help="Input sql databases with genome schema",
    )

    group.add_argument(
        "-o",
        "--output",
        type=Path,
        required=True,
        help="Path to the output database to write",
    )

    return group


def configure_parser(parser: argparse.ArgumentParser):
    """Configure the main parser with argument groups for different namespaces"""

    configure_general_arguments(parser)

    configure_merging_databases(parser)

    parser.set_defaults(run=run)


def run(args: argparse.Namespace, console: Console) -> int:
    """Entry point for the CLI to merge a list of genome databases"""

    oo = create_versioned_path(args.output)

    db = GenomeDatabase.merge_databases(
        output=oo,
        databases=args.input,
        console=console,
    )

    db.conn.close()
    console.print(
        f"[bold green]{'Finished':>12}[/] Merged databases to {str(oo.resolve())}."
    )

    return 0
