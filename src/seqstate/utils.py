from __future__ import annotations

import contextlib
import io
import os
import sys
from contextlib import ExitStack
from pathlib import Path
from typing import BinaryIO, ContextManager, TextIO, overload

if sys.version_info >= (3, 14):
    from compression import zstd
else:
    # backports to previous python versions https://github.com/Rogdham/backports.zstd
    try:
        from backports import zstd  # ty:ignore[unresolved-import]
    except ImportError as err:
        zstd = err  # ty:ignore[invalid-assignment]
        try:
            import zstandard  # ty:ignore[unresolved-import]
        except ImportError as err:
            zstandard = err  # ty:ignore[invalid-assignment]


_BZ2_MAGIC = b"BZh"
_GZIP_MAGIC = b"\x1f\x8b"
_XZ_MAGIC = b"\xfd7zXZ"
_LZ4_MAGIC = b"\x04\x22\x4d\x18"
_ZSTD_MAGIC = b"\x28\xb5\x2f\xfd"


# for running either single or multithread
class DummyPool:

    def map(self, function, iterable):
        return list(map(function, iterable))

    def starmap(self, function, list_of_iterables):
        return list(function(*args) for args in list_of_iterables)

    def __enter__(self):
        return self

    def __exit__(self, exc_value, exc_type, traceback):
        return False


def is_parquet_file(path: Path) -> bool:
    # Check if the file exists
    if not path.is_file():
        raise FileNotFoundError(f"{path} does not exist.")

    with open(path, "rb") as file:
        peek = file.peek(4)

    return peek.startswith(b"PAR1")


def is_sqlite_file(path: Path) -> bool:
    if not path.is_file():
        raise FileNotFoundError(f"{path} does not exist.")

    with open(path, "rb") as f:
        header = f.read(16)

    return header == b"SQLite format 3\0"


@overload
def zopen(
    path,
    mode: str = "rb",
    encoding: None = ...,
    errors: None = ...,
    newline: None = ...,
) -> ContextManager[BinaryIO]: ...


@overload
def zopen(
    path,
    mode: str = "r",
    encoding: str | None = ...,
    errors: str | None = ...,
    newline: str | None = ...,
) -> ContextManager[TextIO]: ...


@contextlib.contextmanager  # ty:ignore[invalid-argument-type]
def zopen(path, mode="r", encoding=None, errors=None, newline=None):
    if mode not in {"r", "rb"}:
        raise ValueError("zopen supports only 'r' or 'rb'")

    with ExitStack() as ctx:
        # open or wrap existing file-like object
        if hasattr(path, "read"):
            raw = path
            if not isinstance(raw, io.BufferedReader):
                raw = io.BufferedReader(raw)
            raw = ctx.enter_context(raw)
        else:
            raw = ctx.enter_context(open(os.fspath(path), "rb"))
            raw = ctx.enter_context(io.BufferedReader(raw))

        peek = raw.peek(8)

        # compression detection
        if peek.startswith(_GZIP_MAGIC):
            try:
                from isal import igzip as gzip
            except ImportError:
                import gzip
            raw = ctx.enter_context(gzip.open(raw, "rb"))

        elif peek.startswith(_BZ2_MAGIC):
            import bz2

            raw = ctx.enter_context(bz2.open(raw, "rb"))

        elif peek.startswith(_XZ_MAGIC):
            import lzma

            raw = ctx.enter_context(lzma.open(raw, "rb"))

        elif peek.startswith(_LZ4_MAGIC):
            try:
                import lz4.frame  # type: ignore
            except ImportError as err:
                raise RuntimeError("LZ4 file but lz4 not installed") from err
            raw = ctx.enter_context(lz4.frame.open(raw))

        elif peek.startswith(_ZSTD_MAGIC):
            if sys.version_info >= (3, 14) or not isinstance(zstd, ImportError):
                raw = ctx.enter_context(zstd.open(raw, "rb"))
            elif not isinstance(zstandard, ImportError):
                dctx = zstandard.ZstdDecompressor()
                raw = ctx.enter_context(dctx.stream_reader(raw))
            else:
                raise ImportError(
                    "Either 'backports.zstd' or 'zstandard' or python 3.14 required!"
                )

        # text mode
        if mode == "r":
            raw = io.TextIOWrapper(
                raw,
                encoding=encoding or "utf-8",
                errors=errors,
                newline=newline,
            )

        yield raw
