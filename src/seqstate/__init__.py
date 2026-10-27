from __future__ import annotations

import copy
import enum
from dataclasses import dataclass, field
from typing import Any

import gb_io
import pyhmmer
import pyrodigal

from ._version import __version__

__author__ = "Raymund Hackett <r.e.hackett@uni-konstanz.de>"
__license__ = "MIT"


class Strand(enum.IntEnum):
    """A flag to declare on which DNA strand a gene is located."""

    Coding = 1
    Reverse = -1

    @property
    def sign(self) -> str:
        """`str`: The strand as a single sign (``+`` or ``-``)."""
        return "+" if self is Strand.Coding else "-"

    @classmethod
    def from_sign(cls, s: str) -> Strand:
        """Create a Strand from '+' or '-'"""
        if s == "+":
            return cls.Coding
        elif s == "-":
            return cls.Reverse
        else:
            raise ValueError(f"Invalid strand sign: {s!r}")


@dataclass(frozen=True)
class Domain:
    """A conserved region within a protein.

    Attributes:
        name (`str`): The accession of the protein domain in the source HMM.
        start (`int`): The start coordinate of the domain within the protein
            sequence (first amino-acid at 0, inclusive).
        end (`int`): The end coordinate of the domain within the protein
            sequence (exclusive).
        hmm (`str`): The name of the HMM library this domain belongs to
            (e.g. ``Pfam``, ``Panther``).
        i_evalue (`float`): The independent e-value reported by ``hmmsearch``
            that measures how reliable the domain annotation is given
            the number of HMMs searched.
        pvalue (`float`): The p-value reported by ``hmmsearch`` that measure
            how likely the domain score is.
        score (`float` | `None`): The bitscore reported by ``hmmsearch``

    NOTE:
        start and end are zero-based indices with an inclusive start, exclusive end
        They describe a half-open interval.
        This is exactly like slices in python would behave but differs from common
        genomic sequence formats! gb-io also uses half-open intervals!
    """

    name: str
    start: int
    end: int
    hmm: str
    i_evalue: float
    pvalue: float
    score: float | None

    def __post_init__(self):
        assert isinstance(self.start, int)
        assert self.start >= 0
        assert isinstance(self.end, int)
        assert self.end >= 0
        assert self.start < self.end
        assert isinstance(self.i_evalue, float)
        assert self.i_evalue > 0
        assert isinstance(self.pvalue, float)
        assert self.pvalue > 0
        if self.score is not None:
            assert isinstance(self.score, float)
            self.score > 0


@dataclass(frozen=True)
class Protein:
    """A sequence of amino-acids translated from a gene.

    Attributes:
        id (`str`): The identifier of the protein.
        seq (`~pyhmmer.easel.DigitalSequence`): The sequence of amino-acids of this protein.
        tt (`int` | `None`): The translation table used. Pyrodigal uses either 4 or 11
        domains (`list` of `~gecco.model.Domain`): A list of domains found
            in the protein sequence.
    """

    id: str  # i could rely on id of digital sequence but its better to keep it seperate
    seq: pyhmmer.easel.DigitalSequence[pyhmmer.easel.AA]
    tt: int | None
    domains: tuple[Domain, ...] = field(default_factory=tuple)

    def __post_init__(self):
        assert isinstance(self.id, str)
        assert isinstance(self.seq, pyhmmer.easel.DigitalSequence), type(self.seq)
        assert self.seq.alphabet == pyhmmer.easel.Alphabet.amino()
        if self.tt is not None:
            assert isinstance(self.tt, int)
        assert isinstance(self.domains, tuple)
        assert all(isinstance(domain, Domain) for domain in self.domains)

    def __eq__(self, other):
        if not isinstance(other, Protein):
            return NotImplemented

        return (
            self.id == other.id
            and self.tt == other.tt
            and self.seq.sequence == other.seq.sequence
            and self.domains == other.domains
        )

    def __hash__(self):
        return hash(
            (
                self.id,
                self.tt,
                self.seq.sequence,
                self.domains,
            )
        )

    def with_domains(self, domains: tuple[Domain, ...]) -> Protein:
        """Return a copy of the protein with the domains replaced"""
        if self.domains:
            raise ValueError(
                "The protein already has domains. This would replace them!"
            )
        return Protein(
            id=self.id,
            seq=self.seq,
            tt=self.tt,
            domains=tuple(domains),
        )

    def disentangle(self) -> Protein:
        """For a protein with more than one domain, return a disentangled copy, else return self"""
        if len(self.domains) <= 1:
            return self

        # greedily pick the most significant domain, remove everything it overlaps,
        # then repeat on the remainder. It needs to be iterative
        # This version produces the correct result.
        # Note that the version in the featurizer is much faster (fully vectorized)
        # but fails when 3 domains overlap only by one domain in the middle
        # NOTE --> Favor this approach even though it is more costly!

        # sort the domains by pvalue to get a fixed order !
        domains = sorted(self.domains, key=lambda d: d.pvalue)
        kept = []
        while domains:
            # place the best p-value domain the the kept list
            best = domains[0]
            kept.append(best)
            # redefine domains list with the remaining domains
            domains = [
                d for d in domains[1:] if d.end < best.start or d.start > best.end
            ]
        kept.sort(key=lambda d: d.start)

        return Protein(
            id=self.id,
            seq=self.seq,
            tt=self.tt,
            domains=tuple(kept),
        )

    def to_faa(self) -> str:
        """
        Return the sequence in FASTA (.faa) format.

        Header format:
        >{self.id}
        Sequence is wrapped to 80 characters per line.
        """
        header = f">{self.id}"
        sequence = self.seq.textize().sequence
        wrapped = "\n".join(sequence[i : i + 80] for i in range(0, len(sequence), 80))

        return f"{header}\n{wrapped}"


@dataclass(frozen=True)
class Gene:
    """A nucleotide sequence coding a protein.

    Attributes:
        source (`~gb_io.Record`): The DNA sequence this gene was
            found in, as a gb_io record.
        start (`int`): The index of the leftmost nucleotide of the gene within
            the source sequence, independent of the strandedness. 0-based, inclusive
        end (`int`): The index of the nucleotide after the gene within
            the source sequence. 0-based, exclusive
        strand (`~Strand`): The strand where the gene is located.
        protein (`~Protein`): The protein translated from this gene.

    NOTE:
        start and end are zero-based indices with an inclusive start, exclusive end
        They describe a half-open interval.
        This is exactly like slices in python would behave but differs from common
        genomic sequence formats!
    """

    source: gb_io.Record
    start: int
    end: int
    strand: Strand
    protein: Protein

    def __post_init__(self):
        assert isinstance(self.source, gb_io.Record)
        assert isinstance(self.start, int)
        assert self.start >= 0
        assert isinstance(self.end, int)
        assert self.end >= 0
        assert self.start < self.end
        assert isinstance(self.strand, Strand)
        assert self.strand in {-1, 1}
        assert isinstance(self.protein, Protein)

    def __eq__(self, other):
        if not isinstance(other, Gene):
            return NotImplemented

        # Compare everything except 'source' because gb_io.Record has no equality
        return (
            self.start == other.start
            and self.end == other.end
            and self.strand == other.strand
            and self.protein == other.protein
        )

    def __hash__(self):
        # hash must match equality semantics
        return hash((self.start, self.end, self.strand, self.protein))

    def copy(self) -> Gene:
        return Gene(
            source=self.source,
            start=self.start,
            end=self.end,
            strand=self.strand,
            protein=self.protein,
        )

    def __copy__(self) -> Gene:
        # shallow copy is identical to copy()
        return self.copy()

    def __deepcopy__(self, memo: dict[int, Any]) -> Gene:
        if id(self) in memo:
            return memo[id(self)]

        new = Gene(
            source=copy.deepcopy(self.source, memo),
            start=self.start,
            end=self.end,
            strand=self.strand,
            protein=self.protein,  # immutable → reuse
        )

        memo[id(self)] = new
        return new

    def with_protein_domains(self, domains: tuple[Domain, ...]) -> Gene:
        """Return a copy of the gene with the protein domains replaced"""
        return Gene(
            start=self.start,
            end=self.end,
            strand=self.strand,
            source=self.source,
            protein=self.protein.with_domains(domains),
        )

    def disentangle(self) -> Gene:
        """For a gene with more than one domain, return a disentangled copy, else return self"""
        return Gene(
            start=self.start,
            end=self.end,
            strand=self.strand,
            source=self.source,
            protein=self.protein.disentangle(),
        )

    def to_feature(self) -> gb_io.Feature:
        location = gb_io.Range(start=self.start, end=self.end)
        if self.strand == Strand.Reverse:
            location = gb_io.Complement(location)

        # create qualifiers
        qualifiers: list[gb_io.Qualifier] = [
            # Add the translation table index
            # NOTE the translation tbl
            gb_io.Qualifier(
                key="transl_table",
                value=str(self.protein.tt),
            ),
            # Add the translated aa sequence
            gb_io.Qualifier(
                key="translation", value=self.protein.seq.textize().sequence
            ),
            # Add the id
            # should be id = f"{record.accession}_{i+1}"
            gb_io.Qualifier(
                key="locus_tag",
                value=self.protein.id,
            ),
            # Add pyrodigal + version as inference
            gb_io.Qualifier(
                key="inference",
                value=f"ab initio prediction: Pyrodigal:{pyrodigal.__version__}",
            ),
        ]

        feature = gb_io.Feature(
            kind="CDS",
            location=location,
            qualifiers=qualifiers,
        )

        return feature

    def to_faa_and_gff(self) -> tuple[str, str]:
        """
        Returns:
            faa_str: FASTA (.faa) formatted protein sequence
            gff_str: single GFF3 line for this CDS
        """
        faa_str = self.protein.to_faa()

        # GFF3
        gff_str = "\t".join(
            [
                str(self.source.accession),
                f"pyrodigal_v{pyrodigal.__version__}",
                "CDS",
                str(self.start + 1),  # reverse start again to compy with 1-based format
                str(self.end),
                ".",
                self.strand.sign,
                "0",
                f"ID={self.protein.id}",
            ]
        )

        return faa_str, gff_str
