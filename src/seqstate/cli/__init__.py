# PYTHON_ARGCOMPLETE_OK
import argparse
from pathlib import Path

import rich
from rich.console import Console
from rich_argparse import ArgumentDefaultsRichHelpFormatter

try:
    import argcomplete
except ImportError:
    argcomplete = None  # type: ignore

from seqstate import __package__ as _module
from seqstate.cli import (build_database, build_labeled_database, build_labels,
                          build_vocab, call_genes, merge_databases, preprocess)


def validate_path_arguments(
    args: argparse.Namespace, parser: argparse.ArgumentParser
) -> int:
    """
    Validate arguments Path type arguments.
    """
    for action in parser._actions:
        # Check if the argument's type is Path
        if action.type is Path:
            # skip output
            if action.dest == "output":
                continue
            arg_value = getattr(args, action.dest)  # Get argument value
            # if required, it must be a file
            if action.required and not arg_value.is_file():
                rich.print(
                    f"[bold red]Error[/]: Required input file '{arg_value}' does not exist."
                )
                return 1
            # if not required but set, it must be a file
            elif (
                not action.required
                and arg_value is not None
                and not arg_value.is_file()
            ):
                rich.print(
                    f"[bold red]Error[/]: Optional input file '{arg_value}' does not exist."
                )
                return 1

    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog=_module,
        formatter_class=ArgumentDefaultsRichHelpFormatter,
        add_help=False,
    )

    commands = parser.add_subparsers(title="commands", required=True, metavar="COMMAND")
    call_genes.configure_parser(
        commands.add_parser(
            "call-genes",
            formatter_class=ArgumentDefaultsRichHelpFormatter,
            help="Call genes in genome fna or genbank files with Pyrodigal.",
        )
    )
    preprocess.configure_parser(
        commands.add_parser(
            "preprocess",
            formatter_class=ArgumentDefaultsRichHelpFormatter,
            help="Preprocess genome fasta or genbank files to Genome sqlite3 dbs.",
        )
    )
    merge_databases.configure_parser(
        commands.add_parser(
            "merge-dbs",
            formatter_class=ArgumentDefaultsRichHelpFormatter,
            help="Merge a list of sql databases in genome schema.",
        )
    )
    build_vocab.configure_parser(
        commands.add_parser(
            "build-vocab",
            formatter_class=ArgumentDefaultsRichHelpFormatter,
            help="Build a vocab file of Pfams to build a Genome database",
        )
    )
    build_labels.configure_parser(
        commands.add_parser(
            "build-labels",
            formatter_class=ArgumentDefaultsRichHelpFormatter,
            help="Build a label file for gene feature labels",
        )
    )
    build_database.configure_parser(
        commands.add_parser(
            "build-database",
            formatter_class=ArgumentDefaultsRichHelpFormatter,
            help="Build an SQL Genome Database",
        )
    )
    build_labeled_database.configure_parser(
        commands.add_parser(
            "build-labeled-database",
            formatter_class=ArgumentDefaultsRichHelpFormatter,
            help="Build an SQL Genome Database labeled with some label information",
        )
    )

    return parser


def run(argv: list[str] | None = None, console: Console | None = None) -> int:
    if console is None:
        console = Console()

    parser = build_parser()
    if argcomplete is not None:
        argcomplete.autocomplete(parser)
    args = parser.parse_args(argv)

    validation_code = validate_path_arguments(args, parser)
    if validation_code != 0:
        return validation_code

    try:
        return args.run(args, console)
    except Exception as err:
        console.print_exception()
        return getattr(err, "code", 1)
    else:
        return 0
