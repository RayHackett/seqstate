from __future__ import annotations

import collections
import contextlib
import copy
import gzip
import warnings
from collections.abc import Callable
from functools import partial
from itertools import chain, groupby
from pathlib import Path
from typing import (IO, Collection, Iterable, Iterator, Literal, Sequence,
                    Sized, TextIO, Union, overload)

# batched was introduced in python 3.12
try:
    from itertools import batched  # ty:ignore[unresolved-import]
except ImportError:
    from itertools import islice

    def batched(iterable, n):
        it = iter(iterable)
        while batch := tuple(islice(it, n)):
            yield batch


import gb_io
import polars
import pyhmmer
import pyrodigal
import rich.progress
from pyhmmer.hmmer import hmmsearch
from rich.console import Console

from seqstate import Domain, Gene, Strand
from seqstate.featureizer import Featureizer
from seqstate.gene_callers import CDSFinder, PyrodigalFinder
from seqstate.hmm_collection import HMMSet
from seqstate.utils import is_parquet_file, zopen
from seqstate.vocab import Labels, Vocab


class SequencePathCollection(Sequence[Path]):
    """A collection of file paths to genomic sequence files"""

    def __init__(
        self,
        paths: Iterable[Path],
        format: str | None = None,
        console: Console | None = None,
    ):

        self.paths = list(paths)
        self.console = console if console is not None else Console()
        self.format = format

        for path in self.paths:
            if not isinstance(path, Path):
                try:
                    path = Path(path)
                except TypeError as e:
                    raise TypeError(
                        f"Expected a str or path-like object, got {type(path).__name__}"
                    ) from e
            if not path.is_file():
                raise FileNotFoundError(f"Path {str(path.resolve())} is missing!")

    def __len__(self) -> int:
        return len(self.paths)

    @overload
    def __getitem__(self, idx: int) -> Path: ...

    @overload
    def __getitem__(self, idx: slice) -> list[Path]: ...

    def __getitem__(self, idx: slice | int):  # ty:ignore[invalid-method-override]
        return self.paths[idx]

    @staticmethod
    def guess_sequences_format(path: Union[str, Path]) -> str | None:
        with zopen(path, mode="rb") as fh:
            for head in iter(lambda: fh.read(256), b""):
                head = head.strip()
                if head:
                    break
            else:
                return None

        if head.startswith(b">"):
            return "fasta"
        if head.startswith(b"LOCUS"):
            return "genbank"
        if head.startswith(b"ID "):
            return "embl"
        if head.startswith(b"##gff-version 3"):
            return "gff"
        return None

    @classmethod
    def _load_records(
        cls,
        path: Path,
        format: str | None = None,
        console: Console | None = None,
        streaming: bool = True,
    ) -> Iterable[gb_io.Record]:

        if console is None:
            console = Console()

        # guess format or use the one given in CLI
        if format is not None:
            format = format.lower()
            console.print("Using user-provided sequence format", repr(format))
        else:
            console.print("Detecting sequence format from file contents")
            format = cls.guess_sequences_format(path)
            if format is None:
                raise RuntimeError(f"Failed to detect format of {path!r}")
            console.print("Detected format of input as", repr(format))

        # load sequences from given file
        try:
            with contextlib.ExitStack() as ctx:
                n = 0
                console.print("Loading sequences from genomic file", repr(path))

                if not streaming:
                    # We can only dislay the progress bar if we are not streaming.
                    # Since we may only have one live progress
                    progress = ctx.enter_context(
                        rich.progress.Progress(
                            *rich.progress.Progress.get_default_columns(),
                            rich.progress.DownloadColumn(),
                            rich.progress.TimeElapsedColumn(),
                            console=console,
                        )
                    )
                    file = ctx.enter_context(
                        progress.open(path, "rb", description="Loading input")
                    )
                else:
                    file = ctx.enter_context(open(path, "rb"))
                decompressed = ctx.enter_context(zopen(file, mode="rb"))

                if format == "genbank":
                    for record in gb_io.iter(decompressed):
                        if record.accession:
                            # NOTE we want only the primary accession number!
                            # TODO this is a bit dirty. probably should work this way
                            record.accession = record.accession.split()[0]
                            yield record
                            n += 1
                        else:
                            raise ValueError(
                                f"Encountered a record without valid accession in {str(path)}"
                            )
                elif format in ["fasta", "embl", "ncbi"]:
                    alphabet = pyhmmer.easel.Alphabet.dna()
                    for seq in ctx.enter_context(
                        pyhmmer.easel.SequenceFile(
                            decompressed, digital=True, alphabet=alphabet
                        )
                    ):
                        yield gb_io.Record(
                            # TODO same as above this is a bit dirty!
                            accession=seq.name.split()[0],  # get primary
                            sequence=seq.alphabet.decode(seq.sequence).encode("ascii"),
                            length=len(seq),
                            molecule_type="dna",
                        )
                        n += 1
                elif format == "faa+gff":
                    alphabet = pyhmmer.easel.Alphabet.amino()
                    for seq in ctx.enter_context(
                        pyhmmer.easel.SequenceFile(
                            decompressed,
                            digital=True,
                            alphabet=alphabet,
                            format="fasta",
                        )
                    ):
                        yield gb_io.Record(
                            # TODO same as above. dirty seq name handling
                            accession=seq.name.split()[0],  # get primary
                            sequence=seq.textize().sequence.lower().encode("ascii"),
                            molecule_type="protein",
                        )
                else:
                    raise NotImplementedError(f"Unknown format {format}")

        except ValueError as err:
            raise RuntimeError(f"Failed to load sequences: {err}") from err
        else:
            console.print(f"Loaded {n} sequences from {str(path)!r}")

    def _iter_sequences(
        self,
        paths: list[Path],
        streaming: bool,
        min_record_length: int = 0,
    ) -> Iterable[gb_io.Record]:
        unique_contigs: set[str] = set()

        discarded_seqs = 0
        for p in paths:
            for record in self._load_records(
                path=p,
                format=self.format,
                console=self.console,
                streaming=streaming,
            ):
                if record.accession is None or not record.accession:
                    raise ValueError(f"A sequence from {str(p)} has no accession!")
                if record.accession in unique_contigs:
                    raise ValueError(f"Duplicate sequence accession {record.accession}")
                unique_contigs.add(record.accession)

                if record.length is not None and record.length >= min_record_length:
                    yield record
                else:
                    discarded_seqs += 1

        if min_record_length > 0:
            self.console.print(
                f"Discarded {discarded_seqs} sequences below {min_record_length} bps"
            )

    def load_files(
        self,
        min_record_length: int = 0,
        streaming: bool = False,
    ) -> SequenceCollection:

        if min_record_length > 0:
            self.console.print(
                f"Skipping sequences with fewer than {min_record_length} bases"
            )

        if self.format == "faa+gff":
            if streaming:
                raise ValueError("No streaming support with faa+gff format!")

            if len(self.paths) != 2:
                raise ValueError("Expecting one fasta faa and one gff file!")

            if min_record_length > 0:
                raise ValueError(
                    "A minimum record length is not supported for faa+gff input"
                )

            for path in self.paths:
                guess = self.guess_sequences_format(path)
                if path.suffix == ".faa" or guess == "fasta":
                    prot_records = list(
                        self._load_records(
                            format="faa+gff",
                            path=path,
                            streaming=False,
                            console=self.console,
                        )
                    )
                    prot_dict = {prot.accession: prot.sequence for prot in prot_records}
                elif path.suffix == ".gff" or guess == "gff":
                    # parse a tab delinated gff file with polars!.
                    annotations = (
                        polars.read_csv(
                            path,
                            has_header=False,
                            separator="\t",
                            comment_prefix="#",
                            null_values=".",
                            schema={
                                "seqid": polars.String,
                                "source": polars.String,
                                "type": polars.String,
                                "start": polars.UInt32,  # 1-based
                                "end": polars.UInt32,  # 1-based
                                "score": polars.Float32,
                                "strand": polars.String,
                                "phase": polars.String,
                                "attributes": polars.String,
                            },
                        )
                        # NOTE we will correct for start and end
                        # at the extract Genes step
                        .filter(polars.col("type") == "CDS")
                        .with_columns(
                            first=polars.col("attributes").str.split(";").list.first()
                        )
                        .filter(polars.col("first").str.starts_with("ID="))
                        .with_columns(
                            protein_id=polars.col("first").str.strip_prefix("ID=")
                        )
                        .select(["seqid", "start", "end", "strand", "protein_id"])
                    )
                    # take only those of type CDS and which have an ID tag
                    # ID tag is the protein accession
                else:
                    raise ValueError(f"Expecting a fasta or gff file from {str(path)}.")

            # we return a sequence collection which can be "called" by the CDS finder.
            # it has feature kind: CDS with qualifiers id and translation
            seqs: dict[str, gb_io.Record] = {}  # dict of seq accession gb_io.Record
            for row in annotations.iter_rows():
                if row[0] not in seqs:
                    seqs[row[0]] = gb_io.Record(sequence=b"", accession=row[0])

                # start -1 to make it 0-based, end exclusive
                location = gb_io.Between(start=row[1] - 1, end=row[2])
                if Strand.from_sign(row[3]) == Strand.Reverse:
                    location = gb_io.Complement(location)
                seqs[row[0]].features.append(
                    gb_io.Feature(
                        kind="CDS",
                        location=location,
                        qualifiers=[
                            # instead I could let CDSFinder assign new names?
                            gb_io.Qualifier(key="locus_tag", value=row[4]),
                            gb_io.Qualifier(
                                key="translation",
                                value=prot_dict.pop(row[4]).decode("ascii"),
                            ),
                            # gb_io.Qualifier(
                            #     key="transl_table",
                            #     value=str(11), # assuming a default here?
                            # ),
                        ],
                    )
                )

            if len(prot_dict) > 0:
                raise ValueError(
                    f"There remain {len(prot_dict)} protein sequences",
                    "not mapped to any contig in the GFF file!",
                )
            return SequenceCollection(
                sequences=list(seqs.values()),
                console=self.console,
            )

        else:
            if streaming:
                return SequenceCollection(
                    sequences=self._iter_sequences(
                        paths=self.paths,
                        streaming=streaming,
                        min_record_length=min_record_length,
                    ),
                    console=self.console,
                )
            else:
                return SequenceCollection(
                    sequences=list(
                        self._iter_sequences(
                            paths=self.paths,
                            streaming=streaming,
                            min_record_length=min_record_length,
                        )
                    ),
                    console=self.console,
                )


class SequenceCollection(Iterable[gb_io.Record]):
    """A collection sequence records in gb_io.Record format"""

    def __init__(
        self,
        sequences: Iterable[gb_io.Record],
        console: Console | None = None,
    ):
        self.sequences = sequences
        self.console = console if console is not None else Console()

    def __iter__(self) -> Iterator[gb_io.Record]:
        return iter(self.sequences)

    def split_and_dump(
        self,
        output: Path,
        compressed: bool = True,
        batch_size: int | None = None,
    ) -> None:
        # lazily split sequences into batches and dump them to a file
        if batch_size is None:
            if not compressed:
                with open(output, "wb") as f:
                    gb_io.dump(self.sequences, fh=f)
            else:
                with gzip.GzipFile(
                    filename=str(output),
                    mode="wb",
                    mtime=0,  # equivalent to gzip -n
                ) as f:
                    gb_io.dump(self.sequences, fh=f)  # ty:ignore[invalid-argument-type]
        else:
            output.mkdir(exist_ok=True, parents=True)
            for i, batch in enumerate(batched(self.sequences, batch_size)):
                batch_filename = output / f"batch{i}.gbf"
                if compressed:
                    batch_filename = batch_filename.with_suffix(".gbf.gz")
                    with gzip.GzipFile(
                        filename=str(batch_filename),
                        mode="wb",
                        mtime=0,
                    ) as f:
                        gb_io.dump(batch, fh=f)  # ty:ignore[invalid-argument-type]
                else:
                    with open(batch_filename, "wb") as f:
                        gb_io.dump(batch, fh=f)
                self.console.print(f"Writing to {str(batch_filename)}")

            self.console.print(
                f"Split the sequences into {i+1} batches of size {batch_size}"
            )

    def extract_genes(
        self,
        min_gene_number: int = 0,
        cds_feature: str = "CDS",
        locus_tag: str = "locus_tag",
        jobs: int = 0,
        call_genes: bool = False,
        streaming: bool = True,
        callback_file: TextIO | None = None,
    ) -> GeneCollection:

        self.console.print("Extracting genes from input sequences")

        if min_gene_number > 0:
            # TODO with the exception of annotation genes!
            self.console.print(
                f"Skipping sequences with fewer than {min_gene_number} genes"
            )

        pyrodial_finder = PyrodigalFinder(metagenome=True, mask=True, n_jobs=jobs)

        if not call_genes:
            self.console.print(f"Using record features named {cds_feature!r}")
            cds_finder = CDSFinder(feature=cds_feature, locus_tag=locus_tag)

        with contextlib.ExitStack() as ctx:

            progress = ctx.enter_context(
                rich.progress.Progress(
                    *rich.progress.Progress.get_default_columns(),
                    rich.progress.MofNCompleteColumn(),
                    rich.progress.TimeElapsedColumn(),
                    console=self.console,
                )
            )

            task = progress.add_task(
                description="Finding ORFs",
                total=(
                    len(self.sequences) if isinstance(self.sequences, Sized) else None
                ),
                unit="contig(s)",
                precision="",
            )

            def callback(
                record: gb_io.Record,
                found: int,
                contig_file: TextIO | None = None,
            ) -> None:
                if contig_file is not None:
                    contig_file.write(f"{record.accession},{found}\n")
                else:
                    self.console.print(
                        f"Found {found} genes in record {repr(record.accession)}"
                    )
                progress.update(task, advance=1)

            if callback_file is not None:
                callback_to_file = partial(callback, contig_file=callback_file)

            # only annotate those sequences lacking features with pyrodigal;
            # parse the rest with the CDS finder
            if not call_genes:
                seq_dict = collections.defaultdict(list)

                for seq in self.sequences:
                    if seq.features is not None and seq.features:
                        seq_dict["annotated"].append(seq)
                    else:
                        seq_dict["unannotated"].append(seq)

                self.console.print(
                    f"Found record features for {len(seq_dict['annotated'])} sequences"
                )

                if seq_dict["unannotated"]:
                    self.console.print(
                        f"Gene-calling unanntotated sequences with pyrodigal {pyrodigal.__version__}"
                    )
                    pyrodigal_genes = pyrodial_finder.find_genes(
                        records=seq_dict["unannotated"],
                        min_gene_number=min_gene_number,
                        progress=callback_to_file if callback_file else callback,
                    )
                else:
                    pyrodigal_genes = []

                # TODO apply min gene number to cds finder too!
                cds_genes = cds_finder.find_genes(
                    records=seq_dict["annotated"], progress=callback
                )

                genes = chain(pyrodigal_genes, cds_genes)

            # use pyrodigal for everything
            else:
                self.console.print(
                    f"Gene-calling all sequences with pyrodigal {pyrodigal.__version__}"
                )

                genes = pyrodial_finder.find_genes(
                    records=self.sequences,
                    min_gene_number=min_gene_number,
                    progress=callback_to_file if callback_file else callback,
                )

            if not streaming:
                genes = list(genes)
                # collect genes because we will materialize them anyway later
                self.console.print(f"Found a total of {len(genes)} genes!")

            return GeneCollection(genes=genes, jobs=jobs, console=self.console)


class GeneCollection(Iterable[Gene]):
    """A collection Gene items"""

    def __init__(
        self,
        genes: Iterable[Gene],
        jobs: int = 0,
        console: Console | None = None,
    ):
        self.genes = genes  # list(genes)
        self.console = console if console is not None else Console()
        self.jobs = jobs

    def __len__(self) -> int:
        if isinstance(self.genes, Sized):
            return len(self.genes)
        else:
            return 0

    # @overload
    # def __getitem__(self, idx: int) -> Gene: ...

    # @overload
    # def __getitem__(self, idx: slice) -> list[Gene]: ...

    # def __getitem__(self, idx: int | slice):
    #     return self.genes[idx]

    def __iter__(self) -> Iterator[Gene]:
        return iter(self.genes)

    def to_featured_sequence_collection(
        self, is_ordered: bool = True
    ) -> SequenceCollection:
        """Use this to include the called genes in the SequeneCollection"""

        # static method to stream featured sequence records from source grouped genes
        def _featurize_records(genes) -> Iterator[gb_io.Record]:
            seen_sequences = set()
            # streaming groupby
            for source, source_genes in groupby(genes, key=lambda gene: gene.source):
                assert source not in seen_sequences, "Genes are not grouped by sequence"
                # clear the existing CDS features because we are going to re-add them
                source.features = [i for i in source.features if i.kind != "CDS"]
                seen_sequences.add(source)
                # sort the source genes
                source_genes = sorted(source_genes, key=lambda gene: gene.start)
                # copy the source; happens once per sequence
                source = source.copy()
                for gene in source_genes:
                    source.features.append(gene.to_feature())

                yield source

        # Streaming option which only works if the genes are already grouped by sequence
        if is_ordered:
            sequences = _featurize_records(self.genes)

        else:
            # Option 2: if the genes are not in order:
            # e.g. if [gene from sequence A, gene from sequence B, gene form sequence A]
            sources = collections.defaultdict(
                copy.deepcopy
            )  # ty:ignore[no-matching-overload]
            for gene in self.genes:
                source_copy = sources[gene.source]
                source_copy.features.append(gene.to_feature())

            for seq in sources.values():
                seq.features.sort(key=lambda feature: feature.location.start)

            sequences = sources.values()

        return SequenceCollection(
            sequences=sequences,
            console=self.console,
        )

    def annotate_domains(
        self,
        hmms: list[HMMSet],
        whitelist: Collection[str] | None = None,
        bit_cutoffs: str | None = None,
        pvalue_filter: float | None = None,
        evalue_filter: float | None = None,
    ) -> AnnotatedGeneCollection:

        self.console.print(
            f"Running HMMER domain annotation with pyhmmer {pyhmmer.__version__}"
        )
        if bit_cutoffs is not None:
            self.console.print(f"Running HMMER with bit cutoffs {bit_cutoffs}")
        if pvalue_filter is not None:
            self.console.print(f"Running HMMER with pvalue_filter {pvalue_filter}")
        if evalue_filter is not None:
            self.console.print(f"Running HMMER with evalue_filter {evalue_filter}")
        if whitelist is not None:
            self.console.print(
                f"Running HMMER with a whitelist of {len(whitelist)} hmms"
            )

        # Run all HMMs over ORFs to annotate with protein domains

        with rich.progress.Progress(
            *rich.progress.Progress.get_default_columns(),
            rich.progress.TimeElapsedColumn(),
            rich.progress.MofNCompleteColumn(),
            console=self.console,
        ) as progress:
            task_hmms = progress.add_task(
                description="Annotating HMMs",
                unit="file" if len(hmms) == 1 else "files",
                total=len(hmms),
                precision="",
            )
            task_domains = progress.add_task(description="", unit="HMMs", precision="")

            def callback(h, t):
                progress.update(task_domains, advance=1)

            for hmm in progress.track(hmms, task_id=task_hmms, total=len(hmms)):
                progress.update(
                    task_domains,
                    description=f" {hmm.id} v{hmm.version}",
                    total=len(whitelist) if whitelist is not None else hmm.size,
                )
                self.console.print(
                    f"Starting annotation with [bold blue]{hmm.id} v{hmm.version}[/]"
                )
                genes = PyHMMER(hmm, self.jobs, whitelist).run(
                    self.genes,
                    progress=callback,
                    bit_cutoffs=bit_cutoffs,
                    pvalue_filter=pvalue_filter,
                    evalue_filter=evalue_filter,
                )
                self.console.print(
                    f"Finished annotation with [bold blue]{hmm.id} v{hmm.version}[/]"
                )
            progress.update(task_id=task_domains, visible=False)

        # Count number of annotated domains
        count = sum(len(gene.protein.domains) for gene in genes)
        self.console.print(f"Found {count} domains across all proteins")

        # NOTE no disentangling;
        # annotating with different sets of HMMs will produce different results
        # after disentangling. therefore I discourage disentangling in general.
        # If nonetheless required, disentangle on the feature table with the Featureizer
        # or (preferred) on an AnnotatedGeneCollection object.

        return AnnotatedGeneCollection(genes=genes, console=self.console)

    def create_gene_table(self) -> GeneTable:
        tbl_data = collections.defaultdict(list)
        for gene in self.genes:
            tbl_data["sequence_id"].append(gene.source.accession)
            tbl_data["protein_id"].append(gene.protein.id)
            tbl_data["start"].append(gene.start)
            tbl_data["end"].append(gene.end)
            tbl_data["strand"].append(gene.strand.sign)

        gene_df = polars.from_dict(
            data=tbl_data,
            schema_overrides={
                "sequence_id": polars.Utf8,
                "protein_id": polars.Utf8,
                "start": polars.UInt32,
                "end": polars.UInt32,
                "strand": polars.Utf8,
            },
        ).sort(by=["sequence_id", "start", "end"])
        return GeneTable(gene_df.lazy())

    def dump(self, output: Path):
        """Dump the GeneCollection to a FASTA (.faa) and GFF3 (.gff) file."""
        output = Path(output)
        output.parent.mkdir(parents=True, exist_ok=True)
        gff_path = output.with_suffix(".gff")
        faa_path = output.with_suffix(".faa")

        with open(faa_path, "w") as faa_f, open(gff_path, "w") as gff_f:
            gff_f.write("##gff-version 3\n")
            for gene in self.genes:
                faa_str, gff_str = gene.to_faa_and_gff()
                faa_f.write(faa_str)
                if not faa_str.endswith("\n"):
                    faa_f.write("\n")
                gff_f.write(gff_str)
                if not gff_str.endswith("\n"):
                    gff_f.write("\n")


class AnnotatedGeneCollection(Sequence[Gene]):
    """A collection Gene items annotated with domain information"""

    def __init__(
        self,
        genes: list[Gene],
        console: Console | None = None,
    ):
        self.genes = genes
        self.console = console if console is not None else Console()

    def __len__(self) -> int:
        return len(self.genes)

    @overload
    def __getitem__(self, idx: int) -> Gene: ...

    @overload
    def __getitem__(self, idx: slice) -> list[Gene]: ...

    def __getitem__(self, idx: int | slice):  # ty:ignore[invalid-method-override]
        return self.genes[idx]

    def disentangle_domains(self) -> AnnotatedGeneCollection:
        self.console.print("Disentangling overlapping domains in each gene")
        return AnnotatedGeneCollection(
            genes=[gene.disentangle() for gene in self.genes],
            console=self.console,
        )

    def create_score_tbl(self) -> polars.DataFrame:
        tbl_data = collections.defaultdict(list)
        for gene in self.genes:
            for domain in gene.protein.domains:
                tbl_data["sequence_id"].append(gene.source.accession)
                tbl_data["protein_id"].append(gene.protein.id)
                tbl_data["start"].append(gene.start)
                tbl_data["end"].append(gene.end)
                tbl_data["strand"].append(gene.strand.sign)
                tbl_data["domain"].append(domain.name)
                tbl_data["hmm"].append(domain.hmm)
                tbl_data["i_evalue"].append(domain.i_evalue)
                tbl_data["pvalue"].append(domain.pvalue)
                tbl_data["domain_start"].append(domain.start)
                tbl_data["domain_end"].append(domain.end)
                tbl_data["score"].append(domain.score)

        schema = {
            "sequence_id": polars.Utf8,
            "protein_id": polars.Utf8,
            "start": polars.UInt32,
            "end": polars.UInt32,
            "strand": polars.Utf8,
            "domain": polars.Utf8,
            "hmm": polars.Utf8,
            "i_evalue": polars.Float64,
            "pvalue": polars.Float64,
            "domain_start": polars.UInt32,
            "domain_end": polars.UInt32,
            "score": polars.Float32,
        }

        return polars.from_dict(
            data=tbl_data,
            schema_overrides=schema,
        ).sort(by=["sequence_id", "start", "end", "domain_start", "domain_end"])

    def create_gene_table(self) -> GeneTable:
        return GeneCollection(self.genes).create_gene_table()

    def create_feature_table(self) -> FeatureTable:
        tbl_data = collections.defaultdict(list)
        for gene in self.genes:
            for domain in gene.protein.domains:
                tbl_data["sequence_id"].append(gene.source.accession)
                tbl_data["protein_id"].append(gene.protein.id)
                tbl_data["start"].append(gene.start)
                tbl_data["end"].append(gene.end)
                tbl_data["strand"].append(gene.strand.sign)
                tbl_data["domain"].append(domain.name)
                tbl_data["hmm"].append(domain.hmm)
                tbl_data["i_evalue"].append(domain.i_evalue)
                tbl_data["pvalue"].append(domain.pvalue)
                tbl_data["domain_start"].append(domain.start)
                tbl_data["domain_end"].append(domain.end)

        schema = {
            "sequence_id": polars.Utf8,
            "protein_id": polars.Utf8,
            "start": polars.UInt32,
            "end": polars.UInt32,
            "strand": polars.Utf8,
            "domain": polars.Utf8,
            "hmm": polars.Utf8,
            "i_evalue": polars.Float64,
            "pvalue": polars.Float64,
            "domain_start": polars.UInt32,
            "domain_end": polars.UInt32,
        }

        return FeatureTable(
            polars.from_dict(
                data=tbl_data,
                schema_overrides=schema,
            )
            .sort(by=["sequence_id", "start", "end", "domain_start", "domain_end"])
            .lazy()
        )


class GeneTable:
    """A Gene Table with all Gene information wrapping around a polars Frame"""

    schema = polars.Schema(
        {
            "sequence_id": polars.Utf8,
            "protein_id": polars.Utf8,
            "start": polars.UInt32,
            "end": polars.UInt32,
            "strand": polars.Utf8,
        }
    )

    def __init__(self, table: polars.LazyFrame):
        self.table = table.select(
            [
                "sequence_id",
                "protein_id",
                "start",
                "end",
                "strand",
            ]
        )
        self._validate()

    def _validate(self):
        if isinstance(self.table, polars.DataFrame):
            if self.table.schema != self.schema:
                raise ValueError(f"Got wrong schema. Expected {self.schema}")

        elif isinstance(self.table, polars.LazyFrame):
            if self.table.collect_schema() != self.schema:
                raise ValueError(
                    f"Got wrong schema {self.table.collect_schema()}. Expected {self.schema}"
                )
        else:
            raise TypeError(f"Expected polars Frame but got {type(self.table)}")

    @classmethod
    def load(cls, file: str | Path) -> GeneTable:
        if is_parquet_file(Path(file)):
            genes = polars.scan_parquet(file)
        else:
            genes = polars.scan_csv(file, separator="\t", schema_overrides=cls.schema)

        genes = genes
        return GeneTable(genes)

    def dump(
        self,
        file: str | Path | IO[bytes],
        format: Literal["parquet", "tsv"] = "parquet",
    ):
        if isinstance(self.table, polars.LazyFrame):
            table = self.table.collect()
        else:
            table = self.table

        if format == "parquet":
            table.write_parquet(file)
        elif format == "tsv":
            table.write_csv(file, separator="\t")
        else:
            raise NotImplementedError

    def label(self, labels: polars.LazyFrame, label_name: str) -> LabeledGeneTable:

        if labels.collect_schema()[label_name] == polars.Boolean:
            mapping = "1:1"
        else:
            mapping = "1:m"

        return LabeledGeneTable(
            table=(
                self.table.join(
                    other=labels,
                    on="protein_id",
                    how="left",
                    validate=mapping,
                ).filter(polars.col(label_name).is_not_null())
            ),
            label_column=label_name,
        )


class LabeledGeneTable:
    """A GeneTable with an additional column for gene-level labels"""

    schema = polars.Schema(
        {
            "sequence_id": polars.Utf8,
            "protein_id": polars.Utf8,
            "start": polars.UInt32,
            "end": polars.UInt32,
            "strand": polars.Utf8,
            "label": polars.String,
        }
    )

    def __init__(
        self,
        table: polars.LazyFrame,
        label_column: str,
    ):
        self.label_column = label_column
        # If the labels are nested (possilby for multilabel - we explode)
        if table.collect_schema()[self.label_column].is_nested():
            table = table.explode(self.label_column)

        self.table = (
            table.select(
                [
                    "sequence_id",
                    "protein_id",
                    "start",
                    "end",
                    "strand",
                    label_column,
                ]
            )
            .with_columns(
                polars.col(label_column).cast(polars.String)  # Cast inner to a string!
            )
            .rename({label_column: "label"})
        )
        self._validate()

    def _validate(self):
        if isinstance(self.table, polars.DataFrame):
            if self.table.schema != self.schema:
                raise ValueError(f"Got wrong schema. Expected {self.schema}")

        elif isinstance(self.table, polars.LazyFrame):
            if self.table.collect_schema() != self.schema:
                raise ValueError(
                    f"Got wrong schema {self.table.collect_schema()}. Expected {self.schema}"
                )
        else:
            raise TypeError(f"Expected polars Frame but got {type(self.table)}")

    def build_label_vocab(
        self,
        min_occurence: int | float | None = None,
    ) -> Labels:
        tbl = self.table.collect()
        tbl = tbl.group_by("label").agg(
            polars.len(),
            fraction=polars.len().truediv(tbl.height),
        )
        if min_occurence is not None:
            if 1 < min_occurence > 0:
                tbl = tbl.filter(polars.col("fraction") >= min_occurence)
            else:
                tbl = tbl.filter(polars.col("len") >= min_occurence)

        return Labels.from_iterator(tbl["label"].sort().to_list())

    @classmethod
    def load(
        cls,
        file: str | Path,
        label_column: str,
    ) -> LabeledGeneTable:

        if is_parquet_file(Path(file)):
            genes = polars.scan_parquet(file)
        else:
            genes = polars.scan_csv(file, separator="\t")

        return LabeledGeneTable(genes, label_column)

    def dump(
        self,
        file: str | Path | IO[bytes],
        format: Literal["parquet", "tsv"] = "parquet",
    ):
        if isinstance(self.table, polars.LazyFrame):
            table = self.table.collect()
        else:
            table = self.table

        table = table.rename({"label": self.label_column})
        if format == "parquet":
            table.write_parquet(file)
        elif format == "tsv":
            table.write_csv(file, separator="\t")
        else:
            raise NotImplementedError


class FeatureTable:
    """A Feature Table with protein domain information wrapping around a polars Frame"""

    schema = polars.Schema(
        {
            "sequence_id": polars.Utf8,
            "protein_id": polars.Utf8,
            "start": polars.UInt32,
            "end": polars.UInt32,
            "strand": polars.Utf8,
            "domain": polars.Utf8,
            "hmm": polars.Utf8,
            "i_evalue": polars.Float64,
            "pvalue": polars.Float64,
            "domain_start": polars.UInt32,
            "domain_end": polars.UInt32,
        }
    )

    def __init__(self, table: polars.LazyFrame):
        self.table = table.select(
            [
                "sequence_id",
                "protein_id",
                "start",
                "end",
                "strand",
                "domain",
                "hmm",
                "i_evalue",
                "pvalue",
                "domain_start",
                "domain_end",
            ]
        )
        self._validate()

    def _validate(self):
        if isinstance(self.table, polars.DataFrame):
            if self.table.schema != self.schema:
                raise ValueError(f"Got wrong schema. Expected {self.schema}")

        elif isinstance(self.table, polars.LazyFrame):
            if self.table.collect_schema() != self.schema:
                raise ValueError(
                    f"Got wrong schema {self.table.collect_schema()}. Expected {self.schema}"
                )
        else:
            raise TypeError(f"Expected polars Frame but got {type(self.table)}")

    @classmethod
    def load(cls, file: str | Path) -> FeatureTable:
        if is_parquet_file(Path(file)):
            features = polars.scan_parquet(file)
        else:
            features = polars.scan_csv(
                file, separator="\t", schema_overrides=cls.schema
            )
        return FeatureTable(features)

    def dump(
        self,
        file: str | Path | IO[bytes],
        format: Literal["parquet", "tsv"] = "parquet",
    ):
        if isinstance(self.table, polars.LazyFrame):
            table = self.table.collect()
        else:
            table = self.table

        if format == "parquet":
            table.write_parquet(file)
        elif format == "tsv":
            table.write_csv(file, separator="\t")
        else:
            raise NotImplementedError

    def build_vocab(
        self,
        genes: GeneTable | None = None,
        domain_number: int | None = None,
        min_domain_occurence: int | None = None,
        only_unique_domains_per_gene: bool = True,
        console: Console | None = None,
        out_database: Path | None = None,
    ) -> tuple[Vocab, polars.DataFrame]:

        if console is None:
            console = Console()

        # I thought quite a bit about the order:
        # -> occurence filter comes last for vocab
        # Thus:

        # NOTE We annotate with genes
        # because currently the data split only splits the gene table
        # the feature table still contains ALL the genes.
        # We do feature selection only on the train set.
        # Otherwise we could skip the annotation step

        if out_database is not None and genes is None:
            raise ValueError(
                "If you want to output a GenomeDatabase, you need to pass genes too!"
            )

        console.print("Filtering out features below a p-value of 1e-9.")
        if domain_number is not None and only_unique_domains_per_gene:
            warnings.warn(
                "Setting 'only_unique_domains_per_gene' implies multilabel while setting 'domain_number' implies multiclass!"
            )
        elif domain_number is not None:
            console.print(f"Selecting the top {domain_number} features per gene.")

        if only_unique_domains_per_gene:
            # multilabel setting!
            # NOTE:
            # This only has an effect on the vocab if min_domain_occurence is set
            # and domains wouldnt land in the vocab
            # unless you count their per-gene duplicates
            featureizer = Featureizer(features=self.table).unique_features_per_gene()
        else:
            featureizer = Featureizer(features=self.table)

        featureizer = featureizer.apply_threshold_filtering(
            threshold_column="pvalue"
        ).select_n_features(n_features=domain_number, threshold_column="pvalue")

        if genes is None:
            console.print(
                "No genes were supplied."
                "Assuming the Feature Table contains only training genes!"
            )
        else:
            genes_df = genes.table.sort(
                by=["sequence_id", "start", "end"]
            ).with_row_index(
                "index", offset=1
            )  # sql is offset by 1

            featureizer = featureizer.annotate(genes=genes_df)

        if min_domain_occurence is not None:
            console.print(
                f"Filtering out features occuring >= {min_domain_occurence} times"
            )
            # as a last step
            # filter out domains which occur less than min_domain_occurence from Vocab
            # In build_database these will be labled as Symbol.UNKNOWN
            domains = (
                featureizer.features.group_by("domain")
                .agg(polars.len())
                .filter(polars.col("len") >= min_domain_occurence)
                .collect()["domain"]
                .sort()
                .to_list()
            )
        else:
            domains = featureizer.features.collect()["domain"].sort().to_list()

        return Vocab.from_iterator(domains), featureizer.features.collect()


class UniversalContainer(object):
    """A container that contains everything."""

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}()"

    def __contains__(self, item: object) -> bool:
        return True


class PyHMMER:
    """A domain annotator that uses `pyhmmer`."""

    def __init__(
        self,
        hmm: HMMSet,
        cpus: int | None = None,
        whitelist: Collection[str] | None = None,
    ) -> None:
        """Prepare a new PyHMMER annotation handler with the given ``hmms``.

        Arguments:
            hmm: Instance of HMM;
            cpus (int, optional): The number of CPUs to allocate for the
                ``hmmsearch`` command. Give ``None`` to use the default.
            whitelist (container of str): If given, a container containing
                the accessions of the individual HMMs to annotate with. If
                `None` is given, annotate with the entire file.

        """
        self.hmm = hmm
        self.cpus = cpus
        self.whitelist = UniversalContainer() if whitelist is None else whitelist

    # TODO should have a check to see if any of the HMMs in the whitelist
    # ARE NOT in the self.hmm file....

    def run(
        self,
        genes: Iterable[Gene],
        progress: Callable[[pyhmmer.plan7.HMM, int], None] | None = None,
        bit_cutoffs: str | None = None,
        evalue_filter: float | None = None,
        pvalue_filter: float | None = None,
    ) -> list[Gene]:
        # collect genes and keep them in original order
        gene_index = list(genes)

        if self.hmm.size is None:
            warnings.warn(
                f"HMM object with name {self.hmm.id} does not have a size. "
                "This might make the i_evalues irreproducible. Do not use them."
            )

        # convert proteins to Easel sequences, namend them after
        # their location in the original input to ignore any duplicate
        # protein identifiers
        esl_sqs = pyhmmer.easel.DigitalSequenceBlock(
            alphabet=pyhmmer.easel.AA(),
        )
        for i, gene in enumerate(gene_index):
            seq = gene.protein.seq
            seq.name = str(i)
            esl_sqs.append(seq)

        # # The hmm collection internally does something like:
        # Only retain the HMMs which are in the whitelist
        # file = ctx.enter_context(zopen(self.hmm.path, mode="rb"))
        # hmm_file = ctx.enter_context(pyhmmer.plan7.HMMFile(file))
        # profiles = (
        #     hmm
        #     for hmm in hmm_file
        #     if self.hmm.relabel(hmm.accession) in self.whitelist
        # )

        profiles = self.hmm.to_HMMCollection().filter_accessions(
            self.whitelist
        )  # ty:ignore[invalid-argument-type]
        # Run search pipeline using the filtered HMMs
        cpus = 0 if self.cpus is None else self.cpus
        hmms_hits = hmmsearch(
            profiles,  # only the hmms in the whitelist
            esl_sqs,
            cpus=cpus,
            callback=progress,
            # We always keep this fixed to the hmm db size
            Z=self.hmm.size,
            bit_cutoffs=bit_cutoffs,
            # domE is the reporting theshold -
            # setting it means we are stricter than the inclusion threshold
            domE=10 if evalue_filter is None else evalue_filter,
        )

        # NOTE normally one would set Z to the total number of sequences
        # This is difficult to account for with batching etc.
        # Instead I missuse Z and fix it to the total number of hmms
        # This means that atleast you get the same results as with hmmscan
        # which is what we actually want.
        # for speed we use hmmsearch and set the Z to hmm.size
        # E-value = P-value * Z, so lower Z means lower E-value

        # NOTE we are taking the reported hits
        # which include all the "included hits"
        # by setting thresholds on reported hits higher than on the included ones
        # we are enforcing that reported == included

        # NOTE we are not loading interpro annotations here but we could
        # Transcribe HMMER hits to GECCO model

        found_domains = collections.defaultdict(list)
        for hits in hmms_hits:
            for hit in hits.reported:
                idx = int(hit.name)
                for domain in hit.domains.reported:
                    # extract name and get InterPro metadata about hit
                    acc = domain.alignment.hmm_accession or domain.alignment.hmm_name
                    accession = self.hmm.relabel(acc)

                    # extract coordinates
                    # pyhmmer seems to follow a 1-based indexing
                    start = domain.alignment.target_from
                    end = domain.alignment.target_to

                    # add the domain to the protein domains of the right gene
                    assert domain.env_from < domain.env_to
                    assert domain.i_evalue >= 0
                    assert domain.pvalue >= 0

                    # applying pvalue filter here - prefer this over e_value
                    if pvalue_filter is not None and domain.pvalue > pvalue_filter:
                        continue

                    found_domains[idx].append(
                        Domain(
                            name=accession,
                            start=start - 1,  # make 0-based
                            end=end,
                            hmm=self.hmm.id,
                            i_evalue=domain.i_evalue,
                            pvalue=domain.pvalue,
                            score=domain.score,
                        )
                    )

        # Create a new Gene with the new protein
        for idx, domains in found_domains.items():
            gene_index[idx] = gene_index[idx].with_protein_domains(domains)

        # return the updated list of genes that was given in argument
        return gene_index
