from __future__ import annotations

from typing import Generic, TypeVar, Union

import polars

T = TypeVar("T", polars.DataFrame, polars.LazyFrame)


class Featureizer(Generic[T]):
    def __init__(self, features: T):
        self.features: T = features

        if "protein_id" not in self.features.collect_schema().names():
            raise ValueError("Feature Table missing 'protein_id' column for Gene names")

    def apply_threshold_filtering(
        self,
        threshold_column: str,
        threshold: float = 1e-9,
    ) -> Featureizer[T]:

        if threshold_column not in self.features.collect_schema().names():
            raise ValueError(
                f"Feature Table missing {threshold_column} column for filtering features"
            )

        # NOTE you filter by pvalue for pfam domains
        # and you filter by "i_evalue" for kofam domains
        self.features = self.features.filter(polars.col(threshold_column) < threshold)

        return self

    def disentangle_domains(
        self,
        start_column: str = "domain_start",
        end_column: str = "domain_end",
        threshold_column: str = "pvalue",
    ) -> Featureizer[T]:

        schema_names = self.features.collect_schema().names()
        if start_column not in schema_names:
            raise ValueError(f"Feature Table missing {start_column} column")
        if end_column not in schema_names:
            raise ValueError(f"Feature Table missing {end_column} column")
        if threshold_column not in schema_names:
            raise ValueError(f"Feature Table missing {threshold_column} column")

        # NOTE
        # this fundamentally differs from the previous iterative model
        # in cases like 3 domains; A, B, C
        # where A and B overlap and B and C overlap but A and C do not.
        # if B is now less significant than either A or C this will remove the overlap
        # The previous pop implementation would keep A and C (if the order is right)

        self.features = (
            self.features.sort(["protein_id", start_column, threshold_column])
            # Track the furthest domain end seen so far among previous rows per protein.
            # If current domain start > cumulative max end of all prior rows per protein,
            # there's no overlap with anything before it.
            .with_columns(
                polars.col(end_column)
                .shift(1)
                .fill_null(polars.col(start_column) - 1)
                .cum_max()
                .over("protein_id")  # per protein
                .alias("_max_end_so_far")
            )
            # A new _group starts when current start > max end so far. These overlap
            .with_columns(
                (polars.col(start_column) > polars.col("_max_end_so_far"))
                .cast(polars.UInt32)
                .cum_sum()
                .over("protein_id")  # per protein
                .alias("_group")
            )
            # Within each overlap group, keep only the row with the min p-value
            .filter(
                polars.col(threshold_column)
                == polars.col(threshold_column).min().over(["protein_id", "_group"])
            )
            .drop(["_max_end_so_far", "_group"])
            .sort(["sequence_id", start_column])
        )

        return self

    def unique_features_per_gene(
        self,
        threshold_column="pvalue",
    ) -> Featureizer[T]:
        """
        Enforce that there is only one unique domain per protein.
        Pick the domain with the lowest threshold_column value.
        """

        if threshold_column not in self.features.collect_schema().names():
            raise ValueError(
                f"Feature Table missing {threshold_column} column for filtering features"
            )

        self.features = self.features.sort(
            by=["protein_id", "domain", "pvalue"], descending=False
        ).unique(subset=["protein_id", "domain"], keep="first")

        return self

    def select_n_features(
        self,
        n_features: int | None,
        threshold_column: str = "i_evalue",
    ) -> Featureizer[T]:

        if threshold_column not in self.features.collect_schema().names():
            raise ValueError(
                f"Feature Table missing {threshold_column} column for filtering features"
            )

        # Do not do feature selection. this selects all features for a gene.
        if n_features is None:
            pass

        elif n_features == 1:
            # select the feature with the lowest evalue per gene
            self.features = self.features.sort(["protein_id", threshold_column]).unique(
                subset="protein_id", keep="first"
            )

        elif n_features > 1:
            # select the n features with the lowest evalue per gene

            self.features = (
                self.features.sort(by=["protein_id", threshold_column])
                .group_by("protein_id", maintain_order=True)
                .head(n=n_features)
            )

            # The above cant use the streaming engine. maybe this is faster?
            # self.features = (
            #     self.features
            #     # sort features within genes oridinally with lowest domain first
            #     .with_columns(
            #         polars.col(threshold_column)
            #         .rank(method="ordinal")
            #         .over("protein_id")
            #         .alias("rank")
            #     )
            #     # only keep up to n features per gene
            #     .filter(polars.col("rank") <= domain_number).sort(
            #         ["protein_id", "rank"]
            #     )
            # )

        else:
            raise ValueError(
                f"domain_number should be a postive integer. Got {n_features}"
            )

        return self

    def annotate(
        self, genes: Union[polars.LazyFrame, polars.DataFrame]
    ) -> Featureizer[T]:

        if "protein_id" not in genes.collect_schema().names():
            raise ValueError("Gene Table missing 'protein_id' column for Gene names")

        if isinstance(genes, polars.DataFrame) and isinstance(
            self.features, polars.DataFrame
        ):
            self.features = genes.join(
                self.features,
                how="inner",
                # on=["sequence_id", "protein_id", "start", "end", "strand"],
                on="protein_id",  # this is a bit unsafer and allows changes in start/end
            )
        elif isinstance(genes, polars.LazyFrame) and isinstance(
            self.features, polars.LazyFrame
        ):
            self.features = genes.join(
                self.features,
                how="inner",
                # on=["sequence_id", "protein_id", "start", "end", "strand"],
                on="protein_id",  # this is a bit unsafer and allows changes in start/end
            )
        else:
            raise TypeError(
                f"Gene table is of type {type(genes)} while feature table is {type(self.features)}"
            )

        return self

    def rank_column(self, column_name: str) -> Featureizer[T]:
        if column_name not in self.features:
            raise ValueError(f"Column {column_name} not found in feature table!")
        self.features = (
            self.features
            # sort features within genes oridinally with lowest domain first
            .with_columns(
                polars.col(column_name)
                .rank(method="ordinal")
                .over("protein_id")
                .alias(f"{column_name}_rank")
            )
        )

        return self
