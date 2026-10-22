from __future__ import annotations

import hashlib
import re
import urllib.request
from collections.abc import Callable
from contextlib import ExitStack
from dataclasses import dataclass
from itertools import chain
from pathlib import Path
from typing import Collection, Iterable, Iterator

try:
    from isal import igzip as gzip
except ImportError:
    import gzip

# compression lib was added to the standard library with python 3.14...
try:
    from compression import zstd  # ty:ignore[unresolved-import]
except ImportError as err:
    zstd = err  # ty:ignore[invalid-assignment]
# prior to python 3.14, using as optional dependency https://pypi.org/project/zstandard/
try:
    import zstandard  # ty:ignore[unresolved-import]
except ImportError as err:
    zstandard = err  # ty:ignore[invalid-assignment]

import pyhmmer
import rich.progress

from seqstate.utils import zopen


@dataclass(frozen=True)
class HMMSet:
    """A Hidden Markov Model definition to use with `pyHMMer`."""

    id: str
    version: str
    url: str | None
    path: str
    size: int | None = None  # number of hmms
    relabel_with: str | None = None
    md5: str | None = None

    def __post_init__(self) -> None:
        file_path = Path(self.path)

        if not file_path.exists():
            raise FileNotFoundError(f"HMM file not found: {file_path}")

        if self.md5 is not None:
            actual_md5 = self._compute_md5(file_path)
            if actual_md5 != self.md5:
                raise ValueError(
                    f"MD5 mismatch for {file_path}: "
                    f"expected {self.md5}, got {actual_md5}"
                )

    @staticmethod
    def _compute_md5(path: Path, chunk_size: int = 8192) -> str:
        hasher = hashlib.md5()
        with path.open("rb") as f:
            for chunk in iter(lambda: f.read(chunk_size), b""):
                hasher.update(chunk)
        return hasher.hexdigest()

    def relabel(self, domain: str) -> str:
        """Rename a domain obtained by this HMM to the right accession.

        This method can be used with HMM libraries that have separate HMMs
        for the same domain, such as Pfam.
        """
        if self.relabel_with is None:
            return domain
        before, after = re.match("^s/(.*)/(.*)/$", self.relabel_with).groups()  # type: ignore
        regex = re.compile(before)
        return regex.sub(after, domain)

    def to_HMMCollection(self) -> HMMCollection:
        return HMMCollection.from_path(Path(self.path))


class HMMCollection(Iterable[pyhmmer.plan7.HMM]):
    """
    A class for lazyly reading, downloading and writing HMM files.
    It lives here in case it is useful during tool setup
    """

    def __init__(self, hmms: Iterable[pyhmmer.plan7.HMM]):
        self.hmms = hmms

    def __iter__(self) -> Iterator[pyhmmer.plan7.HMM]:
        return iter(self.hmms)

    @classmethod
    def download_from_pfam(cls, version: str) -> HMMCollection:
        url = (
            "https://ftp.ebi.ac.uk/pub/databases/Pfam/releases/"
            f"Pfam{version}/Pfam-A.hmm.gz"
        )
        return cls.download_from_url(url)

    @classmethod
    def download_from_url(cls, url: str) -> HMMCollection:
        """Stream HMMs from a (gzipped) HMM file at `url`.

        Note:
            No data is requested until iterating over HMMCollection obj.
            network errors surface during iteration, not here.
        """
        return HMMCollection(HMMCollection._stream_hmms(url))

    @staticmethod
    def _stream_hmms(url: str) -> Iterator[pyhmmer.plan7.HMM]:
        # The ExitStack lives inside the generator so the connection and
        # files stay open for as long as the consumer is iterating.
        with ExitStack() as ctx:
            response = ctx.enter_context(urllib.request.urlopen(url))
            length = response.headers.get("Content-Length")
            dl = ctx.enter_context(
                rich.progress.wrap_file(
                    response,
                    total=int(length) if length else 0,
                    description="Downloading HMMs...",
                )
            )
            src = ctx.enter_context(zopen(dl))
            hmm_file = ctx.enter_context(pyhmmer.plan7.HMMFile(src))
            yield from hmm_file

    @staticmethod
    def _read_hmm_from_file(path: Path | str) -> Iterator[pyhmmer.plan7.HMM]:
        """
        Read HMM from file
        :param path: path to the HMM file
        :return: None
        """
        with ExitStack() as ctx:
            file = ctx.enter_context(zopen(path, mode="rb"))
            hmm_file = ctx.enter_context(pyhmmer.plan7.HMMFile(file))
            for hmm in hmm_file:
                yield hmm

    @classmethod
    def from_path(self, hmmdb_path: Path) -> HMMCollection:
        """
        Read HMMs from either a directory containing HMM files or from a single HMM file
        :param hmmdb_path: path to directory of HMMs or single HMM file
        :return: 'HMM_Loader'
        """
        if hmmdb_path.is_file():
            return HMMCollection(self._read_hmm_from_file(hmmdb_path))
        elif hmmdb_path.is_dir():
            return HMMCollection(
                chain.from_iterable(
                    self._read_hmm_from_file(f) for f in hmmdb_path.iterdir()
                )
            )
        else:
            raise ValueError(
                f"Encounterd path which isnt a file or a dir at {str(hmmdb_path)}"
            )

    @classmethod
    def from_file_with_paths(self, file_with_paths: Path) -> HMMCollection:
        """
        Read HMMs from a file with paths to HMM files
        :param file_with_paths: path to the file with paths
        :return: 'HMM_Loader'
        """
        if not file_with_paths.is_file():
            raise ValueError(f"No file {str(file_with_paths.resolve())} found!")

        with open(file_with_paths, "r") as _in:
            return HMMCollection(
                chain.from_iterable(self._read_hmm_from_file(f.strip()) for f in _in)
            )

    def _filter(self, predicate: Callable[[pyhmmer.plan7.HMM], bool]) -> HMMCollection:
        """Lazily keep only HMMs for which predicate(hmm) is true."""
        return HMMCollection(hmm for hmm in self.hmms if predicate(hmm))

    def filter_accessions(self, accessions: Collection[str]) -> HMMCollection:
        """
        Lazily keep only HMMs whose accession is in the supplied set

        Note:
            this checks `hmm.accession.split(".")[0]`
        """
        return self._filter(lambda hmm: hmm.accession.split(".")[0] in accessions)

    def write_to_h3m_file(self, output: Path):
        """Write the hmms to a single h3m binary file"""
        if output.suffix != ".h3m":
            output = output.with_name(output.name + ".h3m")
        output.parent.mkdir(parents=True, exist_ok=True)

        tmp = output.with_name(output.name + ".tmp")
        n = 0
        try:
            with open(tmp, "wb") as f:
                for hmm in self.hmms:
                    hmm.write(f, binary=True)
                    n += 1
            if n == 0:
                raise ValueError(
                    "No HMMs to write (collection empty or already consumed)"
                )
            tmp.replace(output)  # atomic on the same filesystem
        except BaseException:
            tmp.unlink(missing_ok=True)
            raise

        print(f"Wrote {n} HMMs to {str(output.resolve())}")

    def write_to_compressed_file(self, output: Path):
        """Write the hmms to a single z-standard compressed file"""
        if not output.name.endswith(".hmm.zst"):
            output = output.with_name(output.name + ".hmm.zst")
        output.parent.mkdir(parents=True, exist_ok=True)

        tmp = output.with_name(output.name + ".tmp")
        with ExitStack() as ctx:
            if not isinstance(zstd, ImportError):
                dst = ctx.enter_context(
                    zstd.open(
                        tmp,
                        mode="wb",
                        options={zstd.CompressionParameter.compression_level: 20},
                    )
                )
                print("jupiii")
            elif not isinstance(zstandard, ImportError):
                dst = ctx.enter_context(
                    zstandard.open(
                        tmp, mode="wb", cctx=zstandard.ZstdCompressor(level=20)
                    )
                )
            else:
                raise ImportError(
                    "Either 'zstandard' or 'compression' library required!"
                )

            n = 0
            try:
                for hmm in self.hmms:
                    hmm.write(dst, binary=False)
                    n += 1
                if n == 0:
                    raise ValueError(
                        "No HMMs to write (collection empty or already consumed)"
                    )
                tmp.replace(output)  # atomic on the same filesystem
            except BaseException:
                tmp.unlink(missing_ok=True)
                raise

        print(f"Wrote {n} HMMs to {str(output.resolve())}")
