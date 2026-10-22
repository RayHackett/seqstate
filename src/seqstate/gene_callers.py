import os
import warnings
from collections.abc import Callable
from multiprocessing.pool import Pool, ThreadPool
from typing import Iterable, Iterator, Type, Union

import gb_io
import pyhmmer
import pyrodigal

from seqstate import Gene, Protein, Strand
from seqstate.utils import DummyPool


class CDSFinder:
    """An ORF Finder that simply extracts CDS annotations from records."""

    def __init__(
        self,
        feature: str = "CDS",
        translation_table: int = 11,
        locus_tag: str = "locus_tag",
    ):
        self.feature = feature
        self.translation_table = translation_table
        self.locus_tag = locus_tag

    def find_genes(
        self,
        records: Iterable[gb_io.Record],
        progress: Callable[[gb_io.Record, int], None] | None = None,
    ) -> Iterator[Gene]:
        """Find all genes contained in a sequence of DNA records."""
        ids = set()
        _progress = (lambda x, y: None) if progress is None else progress

        for record in records:
            features = filter(lambda feat: feat.kind == self.feature, record.features)
            i = 1
            for f in features:

                # skip pseudogenes
                if any(q.key == "pseudo" for q in f.qualifiers):
                    continue

                # get the gene translation
                id = f"{record.accession}_{i}"
                tt, prot_seq = None, None
                for q in f.qualifiers:
                    if q.key == "transl_table" and q.value is not None:
                        tt = int(q.value)
                    if q.key == "translation":
                        prot_seq = pyhmmer.easel.TextSequence(
                            sequence=q.value
                        ).digitize(alphabet=pyhmmer.easel.AA())
                    # # get the gene name
                    # if q.key == self.locus_tag and q.value is not None:
                    #     id = q.value
                if tt is None:
                    tt = self.translation_table
                if prot_seq is None:
                    # slice the sequence between start and end and then translate
                    # gb_io uses 0-based, end exclusive intervals just like we do
                    start = min(f.location.start, f.location.end)
                    end = max(f.location.start, f.location.end)
                    try:
                        prot_seq = pyhmmer.easel.DigitalSequence(
                            alphabet=pyhmmer.easel.DNA(),
                            sequence=record.sequence[start:end].upper(),
                        ).translate(genetic_code=pyhmmer.easel.GeneticCode(tt))
                    except ValueError:
                        warnings.warn(
                            f"Invalid sequence in {record.accession} from {start, end}",
                            UserWarning,
                        )
                        # raise ValueError (f"Got sequence {record.sequence[start : end]}") from e
                        prot_seq = prot_seq = pyhmmer.easel.DigitalSequence(
                            alphabet=pyhmmer.easel.DNA(),
                            sequence=b"",
                        ).translate(genetic_code=pyhmmer.easel.GeneticCode(tt))

                protein = Protein(id=id, seq=prot_seq, tt=tt)
                # check IDs are unique
                if protein.id in ids:
                    raise ValueError(
                        f"Duplicate gene identifier found in {record.accession!r}: {protein.id!r}"
                    )
                ids.add(protein.id)
                # wrap the gene into a Gene
                # some f.location might not have the strand
                # if they are nested joins or range types!
                # causes are either unannotated pseudogenes or ribosomal slippage
                if isinstance(f.location, gb_io.Join):
                    strands = [loc.strand for loc in f.location.locations]
                    if strands.count(strands[0]) == len(strands):
                        strand = Strand.from_sign(strands[0])
                    else:
                        raise ValueError(
                            f"Dont know how to handle strand for feature {f}"
                        )
                elif isinstance(f.location, gb_io.Complement) and isinstance(
                    f.location.location, gb_io.Join
                ):
                    strands = [loc.strand for loc in f.location.locations]
                    if strands.count(strands[0]) == len(strands):
                        strand = Strand.from_sign("-")
                    else:
                        raise ValueError(
                            f"Dont know how to handle strand for feature {f}"
                        )
                else:
                    strand = Strand.from_sign(f.location.strand)

                i += 1
                yield Gene(
                    source=record,
                    start=min(f.location.start, f.location.end),
                    end=max(f.location.start, f.location.end),
                    strand=strand,
                    protein=protein,
                )
            _progress(record, i - 1)


class PyrodigalFinder:
    """An ORF Finder that uses the Pyrodigal bindings to Prodigal.

    Prodigal is a fast and reliable protein-coding gene prediction
    for prokaryotic genomes, with support for draft genomes and metagenomes.

    See Also:
        `Doug Hyatt, Gwo-Liang Chen, Philip F. LoCascio, Miriam L. Land,
        Frank W. Larimer and Loren J. Hauser.
        "Prodigal: Prokaryotic Gene Recognition and Translation Initiation
        Site Identification", BMC Bioinformatics 11 (8 March 2010), p119
        <https://doi.org/10.1186/1471-2105-11-119>`_

    """

    def __init__(
        self,
        metagenome: bool = True,
        mask: bool = True,
        n_jobs: int = 1,
    ) -> None:
        """Create a new `PyrodigalFinder` instance.

        Arguments:
            metagenome (bool): Whether or not to run Prodigal in metagenome
                mode, defaults to `True`.
            mask (bool): Whether or not to mask genes running across regions
                containing unknown nucleotides, defaults to `True`.
            n_jobs (int): The number of threads to use to run Pyrodigal in
                parallel. -1 for all, -2 for all but one etc. 0 errors.

        """
        super().__init__()
        self.metagenome = metagenome
        self.mask = mask
        self.n_jobs = n_jobs
        # NOTE by default closed is False
        # which means that partial genes (min length 60) may be called at contig edges
        self.gene_finder = pyrodigal.GeneFinder(meta=metagenome, mask=mask)

    def _sanitize_n_jobs(n_jobs: int | None) -> int:
        if n_jobs is None:
            return 1
        if n_jobs == 0:
            raise ValueError("n_jobs must be nonzero")
        if n_jobs < 0:
            cpus = (
                len(os.sched_getaffinity(0))
                if hasattr(os, "sched_getaffinity")
                else os.cpu_count() or 1
            )
            return max(cpus + 1 + n_jobs, 1)
        return n_jobs

    def _train(self, records: Iterable[gb_io.Record]) -> pyrodigal.TrainingInfo:
        sequences = []
        for i, record in enumerate(records):
            if i > 0:
                sequences.append("TTAATTAATTAA")
            sequences.append(str(record.sequence))
        if len(sequences) > 1:
            sequences.append("TTAATTAATTAA")
        return self.gene_finder.train("".join(sequences))

    def _process_record(
        self, record: gb_io.Record
    ) -> tuple[gb_io.Record, pyrodigal.Genes]:
        return record, self.gene_finder.find_genes(str(record.sequence))

    def find_genes(
        self,
        records: Iterable[gb_io.Record],
        min_gene_number: int = 0,
        progress: Callable[[gb_io.Record, int], None] | None = None,
        *,
        pool_factory: Union[Type[Pool], Callable[[int | None], Pool]] = ThreadPool,
    ) -> Iterator[Gene]:
        """Find all genes contained in a sequence of DNA records.

        Arguments:
            records (iterable of `~gb_io.Record`): An iterable of
                DNA records in which to find genes.
            min_gene_number int, only genes from genomes with more than
                this number of genes will be returned
            progress (callable, optional): A progress callback of signature
                ``progress(record, total)`` that will be called everytime a
                record has been processed successfully, with ``record`` being
                the `~gb_io.Record` instance, and ``total`` being
                the total number of records to process.

        Keyword Arguments:
            pool_factory (`type`): The callable for creating pools, defaults
                to the `multiprocessing.pool.ThreadPool` class, but
                `multiprocessing.pool.Pool` is also supported.

        Yields:
            `~Gene`: An iterator over all the genes found in
            the given records.

        """
        # detect the number of CPUs
        _cpus = self.n_jobs if self.n_jobs > 0 else max(1, len(os.sched_getaffinity(0)))
        _progress = (lambda x, y: None) if progress is None else progress

        # train first if needed
        if not self.metagenome:
            records = list(records)
            self._train(records)

        if _cpus == 1:

            def make_dummy(cpus):
                return DummyPool()

            pool_factory = make_dummy

        # run in parallel using a pool
        with pool_factory(_cpus) as pool:
            for record, genes in pool.imap(self._process_record, records):
                # only return genes from record with at least this number of genes
                if len(genes) >= min_gene_number:
                    _progress(record, len(genes))
                    for j, gene in enumerate(genes):
                        # wrap the protein into a Protein object
                        prot_seq = gene.translate()
                        protein = Protein(
                            id=f"{record.accession}_{j+1}",
                            seq=pyhmmer.easel.TextSequence(
                                sequence=prot_seq,
                            ).digitize(alphabet=pyhmmer.easel.Alphabet.amino()),
                            tt=gene.translation_table,
                        )
                        # wrap the gene into a Gene
                        yield Gene(
                            source=record,
                            start=min(gene.begin, gene.end)
                            - 1,  # to comply with 0-indexing
                            end=max(gene.begin, gene.end),
                            strand=Strand(gene.strand),
                            protein=protein,
                        )
