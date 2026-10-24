from __future__ import annotations

import warnings
from pathlib import Path
from typing import IO, Literal

import polars
from rich.console import Console

from seqstate.featureizer import Featureizer
from seqstate.utils import is_parquet_file
from seqstate.vocab import Labels, Vocab


class GeneTable:
    """A Gene Table with all Gene information wrapping around a polars Frame"""

    schema = polars.Schema(
        {
            "sequence_id": polars.Categorical("sequence_id"),
            "protein_id": polars.Utf8,
            "start": polars.UInt32,
            "end": polars.UInt32,
            "strand": polars.Enum(["+", "-"]),
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

        if len(labels.unique(label_name).collect()) == 2:
            mapping = "1:1"
        else:
            mapping = "1:m"

        return LabeledGeneTable(
            table=(
                self.table.join(
                    other=labels,
                    on=["sequence_id", "protein_id"],
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
            "sequence_id": polars.Categorical("sequence_id"),
            "protein_id": polars.Utf8,
            "start": polars.UInt32,
            "end": polars.UInt32,
            "strand": polars.Enum(["+", "-"]),
            "label": polars.Categorical("label"),
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
            # Cast inner to a string!
            .with_columns(
                polars.col(label_column).cast(polars.Categorical("label"))
            ).rename({label_column: "label"})
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

        genes = genes.cast({label_column: polars.Categorical("label")})

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
            "sequence_id": polars.Categorical("sequence_id"),
            "protein_id": polars.Utf8,
            "start": polars.UInt32,
            "end": polars.UInt32,
            "strand": polars.Enum(["+", "-"]),
            "domain": polars.Categorical("domain"),
            "hmm": polars.Categorical("hmm"),
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
                .collect()["domain"]  # ty:ignore[not-subscriptable]
                .sort()
                .to_list()
            )
        else:
            domains = (
                featureizer.features.collect()["domain"]  # ty:ignore[not-subscriptable]
                .sort()
                .to_list()
            )

        return Vocab.from_iterator(domains), featureizer.features.collect()
