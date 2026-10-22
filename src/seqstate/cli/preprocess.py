import argparse
from pathlib import Path

from rich.console import Console

from seqstate.cli._parser import configure_general_arguments
from seqstate.cli.utils import create_versioned_path
from seqstate.database import GenomeDatabase
from seqstate.hmm_collection import HMMSet
from seqstate.states import SequencePathCollection
from seqstate.vocab import Vocab

"""
This relies on hmms of version:
https://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam37.1/Pfam-A.hmm.gz

# downloaded at:
PFAM=/home/ray/Documents/seqstate/data/hmms/Pfam-A.Pfam37.1.h3m

export PYTHONPATH='/home/ray/Documents/seqstate/src/'
VOCAB=/home/ray/Documents/seqstate/test/test_data/genes_1000.vocab.txt

python -m seqstate.cli preprocess -i ~/Downloads/GCA_000013265.1/genome.gbff -v $VOCAB --hmms $PFAM -o out
"""


def configure_preprocessing(
    parser: argparse.ArgumentParser,
) -> argparse._ArgumentGroup:

    group = parser.add_argument_group(
        title="Preprocessing Arguments",
        description="Going from genomic input, call genes and annotate with HMMs",
    )
    group.add_argument(
        "-i",
        "--input",
        nargs="+",
        type=Path,
        required=True,
        help="Input fna fasta or genbank files.",
    )
    group.add_argument(
        "-v",
        "--vocab",
        type=Path,
        required=False,
        help="Vocab file for Pfams. This is required to build a database.",
    )
    group.add_argument(
        "--hmms",
        type=Path,
        required=True,
        help="Path to hmms to use for domain annotation.",
    )
    group.add_argument(
        "--hmm-type",
        choices=("pfam", "cazy", "other"),
        required=False,
        default="pfam",
        help="The type of hmms passed via --hmms flag.",
    )
    group.add_argument(
        "-o", "--output-dir", type=Path, required=True, help="Output Directory"
    )
    group.add_argument(
        "--output-name",
        type=str,
        required=False,
        default="genome",
        help="The stem to give to the output file(s)",
    )
    group.add_argument(
        "-f",
        "--format",
        choices=("fasta", "genbank", "faa+gff"),
        required=False,
        default=None,
        help="The file format of the input. Optional. Will detect fasta, genbank and embl files.",
    )
    group.add_argument(
        "--min-record-length",
        type=int,
        required=False,
        default=0,
        help="Minimum record length to pass loading.",
    )
    group.add_argument(
        "--min-gene-number",
        type=int,
        required=False,
        default=0,
        help="Minimum number of genes a sequence must have during gene calling. Does not apply to CDS features!",
    )
    group.add_argument(
        "--cds-feature",
        type=str,
        default="CDS",
        required=False,
        help="Extract genes from annotated records using a feature rather than running de-novo gene-calling.",
    )
    group.add_argument(
        "--call-genes",
        action="store_true",
        help="Wether to force re-calling the genes with pyrodigal.",
    )
    group.add_argument(
        "--locus-tag",
        type=str,
        required=False,
        default="locus_tag",
        help=" The name of the feature qualifier to use for naming extracted genes when using the ``--cds-feature`` flag.",
    )
    group.add_argument(
        "-p",
        "--p-filter",
        type=float,
        required=False,
        default=1e-09,
        help="The p-value cutoff for protein domains to be included.",
    )
    group.add_argument(
        "--index-db",
        action="store_true",
        help="Wether to index the returned sql database. False by default to index after merge.",
    )
    group.add_argument(
        "--write-tables",
        action="store_true",
        help="Wether to write gene and feature tables as parquet files too",
    )
    # --disentangle creates issues wrt. HMM annotation because it may lead to
    # different domains being selected and applied to a domain depending on
    # whether we annotate with the full HMM or domains from a whitelist (e.g. a vocab)
    # so we disable and hide it by default.
    group.add_argument(
        "--disentangle",
        action="store_true",
        help=argparse.SUPPRESS,
        # help=(
        #     "Disentangle overlapping domains in each gene by keeping only "
        #     "the domains with the lowest p-value over a given position."
        # ),
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

    configure_preprocessing(parser)

    parser.set_defaults(run=run)


def run(args: argparse.Namespace, console: Console) -> int:
    """Entry point for the CLI to preprocessing"""

    if args.vocab is None and not args.write_tables:
        raise ValueError(
            "Either pass vocab to generate a genome DB or set the write-tables flag to generate some output."
        )

    if args.vocab is not None:
        console.print("Loading Vocab")
        with args.vocab.open() as src:
            vocab = Vocab.load(src)
    else:
        vocab = None

    # NOTE currently only using one set of HMMs!
    if args.hmm_type == "pfam":
        hmms = [
            # anything starting with PF + number gets cut before a .
            HMMSet(
                id="Pfam",
                path=str(args.hmms),
                version="37.1",
                size=23794,
                # md5="768b7d9e42fb08af15cc578f15a89cf9", # for hmm.gz
                md5="3a2921464b634d314568069ddffc6cf0",  # for h3m
                url="ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam37.1/Pfam-A.hmm.gz",
                relabel_with=r"s/(PF\d+)\.\d+/\1/",
            ),
        ]
    elif args.hmm_type == "cazy":
        hmms = [
            HMMSet(
                id="CAZy",
                path=str(args.hmms),
                version="14.0",
                size=875,
                md5="6649757f0bf464d18dc1ad3ed98c956c",  # h3m file
                # md5="8d6d98e77bbca00fbe9cf6744d5b9a4f", # for gzip -n of the txt
                url="http://dbcan-hcc.unl.edu/download/dbCAN-HMMdb-V14.txt",
                relabel_with=r"s/\.hmm$//",
            ),
        ]
    else:
        raise NotImplementedError(
            "Currently you can only annotate with one hmm file."
            "Select one of either 'pfam' or 'cazy'"
        )
        # hmms = [
        #     HMM(
        #         id="user_provided",
        #         path=args.hmm,
        #         version="",
        #         size=875,
        #         url="",
        #     ),
        # ]

    gene_collection = (
        SequencePathCollection(
            args.input,
            format=args.format,
            console=console,
        )
        .load_files(
            streaming=False,
            min_record_length=args.min_record_length,
        )
        .extract_genes(
            min_gene_number=args.min_gene_number,
            cds_feature=args.cds_feature,
            locus_tag=args.locus_tag,
            call_genes=args.call_genes,
            jobs=args.jobs,
            streaming=False,  # materialize the genes since annotation will do so anyway
        )
    )

    # Note that I currently allow individual records without genes.
    if len(gene_collection) == 0:
        raise ValueError("No Genes were found!")

    gene_table = gene_collection.create_gene_table()

    annotated_genes = gene_collection.annotate_domains(
        hmms=hmms,
        whitelist=vocab.classes_ if vocab is not None else None,
        pvalue_filter=args.p_filter,
    )

    if args.disentangle:
        if vocab is None:
            # disentangling with a vocab vs without one produces different results
            # Therefore the user must explicitly ask for disentangling with a vocab
            # to make sure they are aware of this behaviour
            raise ValueError("Disentangling after without a vocab is not encouraged!")
        feature_table = annotated_genes.disentangle_domains().create_feature_table()
    else:
        # by default we do not disentangle
        # this is safer and for a multilabel model, its fine too
        feature_table = annotated_genes.create_feature_table()

    outdir = args.output_dir
    outdir.mkdir(parents=True, exist_ok=True)

    if args.write_tables:
        gt_oo = create_versioned_path(Path(outdir, f"{args.output_name}.genes.parquet"))
        ft_oo = create_versioned_path(
            Path(outdir, f"{args.output_name}.features.parquet")
        )
        gene_table.dump(gt_oo)
        console.print(
            f"[bold green]{'Finished':>12}[/] Wrote GeneTable to {str(gt_oo.resolve())}."
        )
        feature_table.dump(ft_oo)
        console.print(
            f"[bold green]{'Finished':>12}[/] Wrote FeatureTable to {str(gt_oo.resolve())}."
        )

    # NOTE only with a vocab can we build a database!
    if vocab is not None:
        db_oo = create_versioned_path(Path(outdir, f"{args.output_name}.sqlite3"))
        db = GenomeDatabase.build_database(
            genes=gene_table,
            features=feature_table,
            vocab=vocab,
            output=db_oo,
            console=console,
            domain_number=args.domain_number,
            only_unique_domains_per_gene=args.unique_domains,
            index_db=args.index_db,  # false by default to disable indexing
        )

    db.close()
    console.print(
        f"[bold green]{'Finished':>12}[/] Built databases to {str(db_oo.resolve())}."
    )
    # exit successfully
    return 0
