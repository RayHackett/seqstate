"""Vocab and special symbol handling."""

from __future__ import annotations

import enum
import io
from abc import ABC, abstractmethod
from typing import Iterable, Iterator, TextIO, Type, TypeVar

import polars

from seqstate.featureizer import Featureizer


class Symbol(enum.IntEnum):
    """A special symbol to use for padding or masking."""

    PAD = 0
    UNKNOWN = 1
    MASK = 2
    CLS = 3


class SpecialLabel(enum.IntEnum):
    """A special unknown label to use for classification tasks"""

    UNKNOWN = 0


T = TypeVar("T", bound="_Index")


class _Index(ABC):
    def __init__(self, classes: Iterable[str] = ()) -> None:
        """Create a new vocabulary with the given classes."""
        self.classes_ = [str(i) for i in classes]
        if len(set(self.classes_)) != len(self.classes_):
            raise ValueError("Duplicate Tokens in Vocab")

        self.lookup_table = {x: i for i, x in enumerate(self.classes_)}

    @property
    @abstractmethod
    def default(self) -> int:
        """`int`: The default index for unknown symbols must be set in child class."""

    @classmethod
    def loads(cls: Type[T], data: str) -> T:
        """
        Load the vocab from a string.
        """
        lines = data.splitlines()
        return cls(lines)

    @classmethod
    def load(cls: Type[T], file: TextIO) -> T:
        """
        Load the vocab from a file.
        """
        lines = file.read().splitlines()
        return cls(lines)

    def __len__(self) -> int:
        """Return the number of classes in the vocabulary."""
        return len(self.classes_)

    def __reduce__(self):
        return type(self), (self.classes_,)

    def __hash__(self):
        return hash((tuple(self.classes_)))

    def dumps(self) -> str:
        """
        Dump the vocab to a string.
        """
        buffer = io.StringIO()
        self.dump(buffer)
        return buffer.getvalue()

    def dump(self, file: TextIO):
        """
        Write the vocab to a file object.
        """
        for token in self.classes_:
            file.write(token + "\n")

    def items(self) -> Iterator[tuple[str, int]]:
        """Iterate over the symbol/value pairs in the vocab."""
        for symbol in self.classes_:
            yield symbol, self.lookup_table[symbol]

    def lookup_tokens(self, indices: Iterable[int]) -> list[str]:
        """Lookup the class names corresponding to several indices."""
        return [self.classes_[i] for i in indices]

    def lookup_indices(self, tokens: Iterable[str]) -> list[int]:
        """Lookup the indices corresponding to several class names."""
        default = self.default
        return [self.lookup_table.get(str(x), default) for x in tokens]

    @abstractmethod
    def __getitem__(self, symbol: str) -> int:
        """Implement in child class"""

    @abstractmethod
    def from_iterator(
        cls,
        tokens: Iterable[str],
    ) -> _Index:
        """Implement in child class"""


class Labels(_Index):
    @classmethod
    def from_iterator(
        cls,
        tokens: Iterable[str],
    ) -> Labels:
        """Create a new vocabulary from an iterable of symbols.

        Feature names will be deduplicated and sorted, and special
        symbols will be added at the beginning of the class list.

        """
        symbols = sorted(SpecialLabel.__members__.values())
        symbol_names = {s.name for s in symbols}

        # Filter out anything that is a special symbol
        unique = sorted({t for t in tokens if t not in symbol_names})
        classes = [symbol.name for symbol in symbols] + unique
        return cls(classes)

    def __init__(self, classes: Iterable[str] = ()) -> None:
        """Create a new vocabulary with the given classes."""
        super().__init__(classes=classes)

        if [
            symbol.name for symbol in SpecialLabel.__members__.values()
        ] != self.classes_[: len(SpecialLabel.__members__)]:
            raise ValueError("First elements of Vocab do not match the special Tokens")

    def __eq__(self, other):
        if not isinstance(other, Labels):
            return NotImplemented

        return tuple(self.classes_) == tuple(other.classes_)

    def __getitem__(self, symbol: str) -> int:
        """Get the index corresponding to a symbol in the vocabulary."""
        if isinstance(symbol, SpecialLabel):
            symbol = symbol.name
        return self.lookup_table.get(symbol, self.default)

    @property
    def default(self) -> int:
        """`int`: The default index for unknown symbols."""
        return SpecialLabel.UNKNOWN

    @default.setter
    def default(self, index: int) -> None:
        if self.classes_[index] != SpecialLabel.UNKNOWN.name:
            raise ValueError("Invalid default index")

    @classmethod
    def build_labels(
        cls,
        features: polars.LazyFrame,
        genes: polars.LazyFrame,
        min_label_occurence: int,
        label_col: str = "domain",
    ) -> Labels:

        # I thought quite a bit about the order:
        # -> occurence filter comes last for vocab
        # Thus:

        # NOTE We annotate with genes
        # because currently the data split only splits the gene table
        # the feature table still contains ALL the genes.
        # We do feature selection only on the train set.
        # Otherwise we could skip the annotation step

        # TODO this assumes the labels have an i_evalue!!

        featureizer = (
            Featureizer(features=features)
            .apply_threshold_filtering(threshold_column="i_evalue")
            .select_n_features(n_features=1, threshold_column="i_evalue")
            .annotate(genes=genes)
        )

        assert isinstance(featureizer.features, polars.LazyFrame)

        # as a last step
        # filter out labels which occur less than min_label_occurence from Label vocab
        # In build_database these will be labeled as Symbol.UNKNOWN
        common_domains = (
            featureizer.features.group_by(label_col)
            .agg(polars.len())
            .filter(polars.col("len") >= min_label_occurence)
            .collect()[label_col]  # ty:ignore[not-subscriptable]
            .sort()
            .to_list()
        )

        return cls.from_iterator(common_domains)


class Vocab(_Index):
    """A simple vocabulary using a `dict` as a lookup table."""

    @classmethod
    def from_iterator(
        cls,
        tokens: Iterable[str],
    ) -> Vocab:
        """Create a new vocabulary from an iterable of symbols.

        Feature names will be deduplicated and sorted, and special
        symbols will be added at the beginning of the class list.

        """
        unique = sorted(set(tokens))
        symbols = sorted(Symbol.__members__.values())
        classes = [symbol.name for symbol in symbols] + unique
        return cls(classes)

    def __init__(self, classes: Iterable[str] = ()) -> None:
        """Create a new vocabulary with the given classes."""
        super().__init__(classes=classes)

        if [symbol.name for symbol in Symbol.__members__.values()] != self.classes_[
            : len(Symbol.__members__)
        ]:
            raise ValueError("First elements of Vocab do not match the special Tokens")

    def __eq__(self, other):
        if not isinstance(other, Vocab):
            return NotImplemented

        return tuple(self.classes_) == tuple(other.classes_)

    def __getitem__(self, symbol: str) -> int:
        """Get the index corresponding to a symbol in the vocabulary."""
        if isinstance(symbol, Symbol):
            symbol = symbol.name
        return self.lookup_table.get(symbol, self.default)

    @property
    def default(self) -> int:
        """`int`: The default index for unknown symbols."""
        return Symbol.UNKNOWN

    @default.setter
    def default(self, index: int) -> None:
        if self.classes_[index] != Symbol.UNKNOWN.name:
            raise ValueError("Invalid default index")
