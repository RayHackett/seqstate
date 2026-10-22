import argparse
from contextlib import ExitStack
from pathlib import Path

from rich.console import Console

from seqstate.cli._parser import configure_general_arguments
from seqstate.cli.utils import create_versioned_path
from seqstate.states import SequencePathCollection

"""
export PYTHONPATH='/home/ray/Documents/seqstate/src/'

python -m seqstate.cli call-genes -i ~/Downloads/GCA_000013265.1/genome.gbff -o out --type faa+gff
"""


def configure_call_genes_input(
    parser: argparse.ArgumentParser,
) -> argparse._ArgumentGroup:

    group = parser.add_argument_group(
        "Build Database",
        "Call genes with Pyrodigal",
    )

    group.add_argument(
        "-i",
        "--input",
        nargs="+",
        required=True,
        type=Path,
        default="Input genomic sequence file in gbf, fna, or embl format",
    )
    group.add_argument(
        "-o",
        "--output-dir",
        required=True,
        type=Path,
        default="Input genomic sequence file in gbf, fna, or embl format",
    )
    group.add_argument(
        "--output-name",
        type=str,
        required=False,
        default="genome",
        help="The stem to give to the output file(s)",
    )
    group.add_argument(
        "--type",
        required=True,
        nargs="+",
        type=str,
        default="What type of output to write the genes to. any of 'table', 'gbk', 'faa+gff', 'contig'",
    )

    return group


def configure_parser(parser: argparse.ArgumentParser):
    """Configure the main parser with argument groups for different namespaces"""

    configure_general_arguments(parser)

    configure_call_genes_input(parser)

    parser.set_defaults(run=run)


def run(args: argparse.Namespace, console: Console) -> int:
    """Entry point for the CLI to call_genes"""

    accepted = set(["table", "gbk", "faa+gff", "contig"])
    for i in args.type:
        if i not in accepted:
            raise ValueError(f"Unexpected output type {i}. Accepting only {accepted}")

    outdir = args.output_dir
    outdir.mkdir(parents=True, exist_ok=True)

    with ExitStack() as ctx:

        if "contig" in args.type:
            contig_file = ctx.enter_context(
                open(Path(outdir, f"{args.output_name}.contigs.csv"), "w")
            )
        else:
            contig_file = None

        gene_collection = (
            SequencePathCollection(
                args.input,
                console=console,
            )
            .load_files(
                streaming=True,
            )
            .extract_genes(
                call_genes=True,
                jobs=args.jobs,
                streaming=len(args.type)
                == 1,  # stream when there is only one output type
                callback_file=contig_file,
            )
        )

    if "faa+gff" in args.type:
        # rn this only checks for the version on gff...
        # it might overwrite the faa file...
        name = Path(outdir, f"{args.output_name}.genes.gff")
        oo = create_versioned_path(name)
        gene_collection.dump(oo)

    if "gbk" in args.type:
        name = Path(outdir, f"{args.output_name}.genes.gbk")
        oo = create_versioned_path(name)
        featured_seqs = gene_collection.to_featured_sequence_collection()
        featured_seqs.split_and_dump(
            oo,
            compressed=False,
        )

    if "table" in args.type:
        name = Path(outdir, f"{args.output_name}.genes.parquet")
        oo = create_versioned_path(name)
        gene_table = gene_collection.create_gene_table()
        gene_table.dump(oo)

    console.print(
        f"[bold green]{'Finished':>12}[/] Wrote called genes to {str(oo.resolve())}."
    )

    # db.close()
    # exit successfully
    return 0
