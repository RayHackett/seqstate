from __future__ import annotations

import importlib.resources
import os
import re
import sqlite3
import time
import warnings
from contextlib import ExitStack

try:
    from importlib.resources.abc import \
        TraversableResources  # ty:ignore[unresolved-import]
except ImportError:  # Python < 3.11
    from importlib.abc import TraversableResources

from pathlib import Path

import polars
import rich
import rich.progress
from rich.console import Console

from seqstate.featureizer import Featureizer
from seqstate.tables import FeatureTable, GeneTable, LabeledGeneTable
from seqstate.vocab import Labels, SpecialLabel, Symbol, Vocab

# strand conversion to integer
STRANDS = polars.DataFrame(
    [
        {"strand": "+", "strand_int": 1},
        {"strand": "-", "strand_int": -1},
    ],
    schema={
        "strand": polars.Enum(["+", "-"]),
        "strand_int": polars.Int8,
    },
)


def check_db_schema(
    db: TraversableResources | str | sqlite3.Connection,
    schema_file: TraversableResources,
) -> bool:
    """
    Check if the database has the correct schema.

    Notes
    - If `db` is a sqlite3.Connection, it is used directly and not closed.
    - If `db` is a filesystem path, a temporary read-only connection is opened.
    - If `db` is a Traversable, it is resolved via importlib.resources.as_file()
        (zero-copy when possible, temp extraction otherwise).
    """

    conn: sqlite3.Connection
    with ExitStack() as ctx:
        if isinstance(db, sqlite3.Connection):
            conn = db

        elif isinstance(db, (str, os.PathLike)):
            if str(db) == ":memory:":
                raise ValueError(
                    "In-memory SQLite databases require passing an active "
                    "sqlite3.Connection, not ':memory:' as a path."
                )

            path = Path(db)  # type: ignore
            if not path.is_file():
                raise FileNotFoundError(f"No database at {path} found!")

            conn = sqlite3.connect(f"file:{path.as_posix()}?mode=ro", uri=True)
            ctx.callback(conn.close)

        elif isinstance(db, TraversableResources):
            if not db.is_file():
                raise ValueError("Traversable database must be a file")

            path = ctx.enter_context(importlib.resources.as_file(db))
            conn = sqlite3.connect(
                f"file:{path.as_posix()}?mode=ro",
                uri=True,
            )
            ctx.callback(conn.close)

        else:
            raise TypeError(f"Unsupported db type: {type(db)!r}")

        with schema_file.open() as f:
            lines = f.readlines()

        cursor = conn.execute("SELECT sql FROM sqlite_schema WHERE type='table'")
        db_schema = [row[0] for row in cursor.fetchall() if row[0]]

    # read the schema file and split it on semicolons, removing comment lines
    clean_lines = [line.strip() for line in lines if not line.strip().startswith("--")]
    file_schema = " ".join(clean_lines).split(";")

    def normalize(block: str) -> str:
        block = re.sub(r"\s+", " ", block).replace("\n", " ").strip()

        # normalize case for keywords
        block = block.upper()

        # remove optional clauses that don't affect schema structure
        block = block.replace("CREATE TABLE IF NOT EXISTS", "CREATE TABLE")

        return block

    file_schema = [
        normalize(block)
        for block in file_schema
        if block.strip().startswith("CREATE TABLE")
    ]

    db_schema = [
        normalize(block)
        for block in db_schema
        if not block.startswith("CREATE TABLE sqlite")
    ]

    file_schema.sort()
    db_schema.sort()

    if db_schema != file_schema:
        print("Required Schema in Seqstate Module")
        print(file_schema)
        print("Schema in the Database")
        print(db_schema)

    return db_schema == file_schema


class GenomeDatabase:
    """
    An SQL GenomeDatabase. It uses the genome schema.

    NOTE
        This Database may have multiple occurences of the same domain for a given Gene.
        This is the desired behaviour for a MutliClass Model. For a Multilabel Model,
        the query planner in SQLMultilabelDataset enforces uniqueness.
        However this might slightly through of `num_annotated_genes` and
        `domain_occurences`

    """

    SCHEMA_FILE = "schema_genome.sql"

    @classmethod
    def schema(cls) -> TraversableResources:
        return importlib.resources.files(__package__).joinpath(cls.SCHEMA_FILE)

    def __init__(
        self,
        path: str | os.PathLike | Path,
        conn: sqlite3.Connection,
    ):
        self.path = path
        self._conn = conn

        if not check_db_schema(db=self._conn, schema_file=self.schema()):
            raise ValueError(
                f"The DB at {str(self.path)} is not in required genome schema!"
            )

    @property
    def conn(self) -> sqlite3.Connection:
        return self._conn

    @classmethod
    def open(cls, path: Path | str | os.PathLike) -> GenomeDatabase:
        if not Path(path).is_file():
            raise FileNotFoundError(f"File {str(path)} does not exist!")
        conn = sqlite3.connect(path)
        return cls(path, conn)

    def close(self) -> None:
        self._conn.close()

    def __enter__(self):
        return self

    def __exit__(self, *_):
        self.close()

    def __len__(self) -> int:
        return self._conn.execute("SELECT MAX(id) FROM Gene").fetchone()[0]

    @property
    def vocab_size(self) -> int:
        return self._conn.execute("SELECT COUNT(*) FROM Vocab").fetchone()[0]

    @property
    def num_annotated_genes(self) -> int:
        # total instances = number of unique annotated genes in domain table
        # for multilabel, each domain should only occur once per gene!
        return self._conn.execute(
            "SELECT COUNT(DISTINCT gene) FROM Domain;"
        ).fetchone()[0]

    @property
    def domain_occurences(self) -> dict[int, int]:
        # Count occurrences of each domain in domain table.
        # Technically this should be done only on the train set without the validation
        # The differences are likely very minor - order of magnitude is correct
        return dict(
            self._conn.execute(
                """
                SELECT domain, COUNT(*) 
                FROM Domain 
                GROUP BY domain
                """
            ).fetchall()
        )

    @property
    def padding_index(self) -> int:
        return self._conn.execute(
            "SELECT id FROM Vocab WHERE Vocab.name = 'PAD'"
        ).fetchone()[0]

    @property
    def max_special_index(self) -> int:
        specials = [i.__repr__() for i in Symbol.__members__.keys()]
        specials_str = ",".join(specials)
        query = f"""
            SELECT MAX(id)
            FROM Vocab
            WHERE Vocab.name IN ({specials_str});
        """
        return self._conn.execute(query).fetchone()[0]

    @property
    def vocab(self):
        query = "SELECT name FROM Vocab ORDER BY Vocab.id"
        return Vocab((name for (name,) in self._conn.execute(query)))

    @classmethod
    def _build_sql_db(
        cls,
        genes_df: polars.DataFrame,
        features_df: polars.DataFrame,
        vocab: Vocab,
        console: Console,
        output: str | os.PathLike,
        index_db: bool = True,
    ) -> GenomeDatabase:

        # sort genes and compute absolute index
        if "index" not in genes_df.columns:
            genes_df = genes_df.sort(by=["sequence_id", "start", "end"]).with_row_index(
                "index", offset=1
            )  # sql is offset by 1

        genes_df = (
            genes_df.join(
                STRANDS,
                left_on="strand",
                right_on="strand",
                how="inner",
            )
            .drop("strand")
            .rename({"strand_int": "strand"})
            .sort(by=["sequence_id", "start", "end"])
            .select(["index", "sequence_id", "protein_id", "start", "end", "strand"])
        )

        vocab_table = polars.DataFrame(
            {
                "domain": vocab.classes_,
                "vocab_id": vocab.lookup_indices(vocab.classes_),
            },
            schema={
                "domain": polars.Categorical("domain"),
                "vocab_id": polars.UInt16,
            },
        )
        features_df = (
            features_df.join(
                vocab_table, how="left", on="domain", maintain_order="left"
            )
            .fill_null(vocab.default)
            .select(["index", "vocab_id", "pvalue", "domain_start"])
        )

        console.print(f"Number of Genes: {genes_df.height}")
        console.print(f"Number of Contigs: {genes_df.n_unique('sequence_id')}")
        console.print(
            "Fraction of annotated Genes:",
            f"{features_df.n_unique('index')/genes_df.height}",
        )
        console.print(
            "Fraction of UNKNOWN Domains:",
            f"{features_df.filter(polars.col('vocab_id')==vocab.default).height/features_df.height}",
        )
        db = sqlite3.connect(output)

        # load database schema
        with cls.schema().open() as src:
            schema = src.read()

        with rich.progress.Progress(console=console) as progress:
            db.executescript(schema)
            db.executemany(
                "INSERT INTO Vocab VALUES(?, ?);",
                progress.track(
                    ((i, x) for (x, i) in vocab.items()),
                    total=len(vocab),
                    description="Vocab",
                ),
            )
            db.executemany(
                "INSERT INTO Gene VALUES(?, ?, ?, ?, ?, ?);",
                progress.track(
                    genes_df.iter_rows(), total=genes_df.height, description="Genes"
                ),
            )
            db.executemany(
                "INSERT INTO Domain VALUES(?, ?, ?, ?);",
                progress.track(
                    features_df.iter_rows(),
                    total=features_df.height,
                    description="Domains",
                ),
            )

        if index_db:
            start = time.time()
            console.print("Indexing Gene Table ...")
            db.execute("CREATE INDEX GeneContig ON Gene(contig);")
            console.print("Indexing Domain Table ...")
            db.execute("CREATE INDEX DomainGene ON Domain(gene);")
            db.execute("CREATE INDEX DomainPvalue ON Domain(pvalue);")
            db.execute("CREATE INDEX DomainStart ON Domain(start);")
            console.print(f"Indexing time: {time.time()-start}")

        db.commit()

        return cls(output, db)

    @classmethod
    def build_database(
        cls,
        genes: GeneTable,
        features: FeatureTable,
        vocab: Vocab,
        output: str | os.PathLike,
        index_db: bool = True,
        domain_number: int | None = None,
        only_unique_domains_per_gene: bool = True,
        console: Console | None = None,
    ) -> GenomeDatabase:
        """
        Create a genome sqlite3 database from a gene and domain table with a vocab.

        Arguments:
            genes: GeneTable
            features: FeatureTable
            vocab: Vocab instance
            output:
                - str | Path: path where the db will be created, or ":memory:" for in-memory
            index_db: `bool` Index the database
            domain_number: Optional, `int` how many domains per gene to save
            only_unique_domains_per_gene: `bool` save only unique domains per gene
            console: Optional rich.console.Console to print progress/debug messages

        Returns:
            GenomeDatabase

        """
        if console is None:
            console = Console()  # prints to stdout. to supress pass a dummy console

        if domain_number is not None and only_unique_domains_per_gene:
            warnings.warn(
                "Setting 'only_unique_domains_per_gene' implies multilabel while setting 'domain_number' implies multiclass!"
            )

        if isinstance(output, (str, os.PathLike)) and output != ":memory:":
            if Path(output).exists():
                raise ValueError(f"File already {str(output)!r} exists!")

            Path(output).parent.mkdir(parents=True, exist_ok=True)

        # sort genes and compute absolute index
        genes_df = genes.table.sort(by=["sequence_id", "start", "end"]).with_row_index(
            "index", offset=1
        )  # sql is offset by 1

        console.print("Loading Features")
        featureizer = Featureizer(features.table)
        if only_unique_domains_per_gene:
            console.print("Creating db with unique features per gene.")
            featureizer = Featureizer(features.table).unique_features_per_gene()
        else:
            console.print("Creating db with all features per gene.")
            featureizer = Featureizer(features.table)

        # join feature and gene table and filtering to below 1e-9
        features_df = (
            featureizer.apply_threshold_filtering(threshold_column="pvalue")
            .select_n_features(n_features=domain_number, threshold_column="pvalue")
            .annotate(genes_df)  # inner join maintaining gene order
            .features.collect()
        )

        console.print("Collecting Genes")
        genes_df = genes_df.collect()

        return cls._build_sql_db(
            genes_df=genes_df,
            features_df=features_df,
            vocab=vocab,
            console=console,
            output=output,
            index_db=index_db,
        )

    @classmethod
    def _merge_database(
        cls,
        db: sqlite3.Connection,
        src_path: Path,
        first: bool = False,
    ):
        """Merge the database located at src_path into the supplied db connection."""

        if not check_db_schema(db=src_path, schema_file=cls.schema()):
            raise ValueError(
                f"The DB at {str(src_path)} is not in required genome schema!"
            )

        cur = db.cursor()
        cur.execute("ATTACH DATABASE ? AS src", (str(src_path),))

        def vocab_fingerprint(
            cur: sqlite3.Cursor, schema: str
        ) -> list[tuple[str, str]]:
            return cur.execute(
                f"""
                SELECT id, name FROM {schema}.Vocab ORDER BY id
            """
            ).fetchall()

        # check if the voacbs really are identical
        if first:
            cur.execute(
                """
                INSERT INTO main.Vocab(id, name)
                SELECT id, name FROM src.Vocab ORDER BY id
            """
            )
        else:
            src_vocab = vocab_fingerprint(cur, "src")
            main_vocab = vocab_fingerprint(cur, "main")
            if src_vocab != main_vocab:
                raise RuntimeError(f"Vocab mismatch in {src_path}")

        cur.execute("PRAGMA foreign_keys = OFF")

        # get number of genes in table to calculate an offset
        # First one will be empty
        offset = cur.execute("SELECT COALESCE(MAX(id), 0) FROM main.Gene").fetchone()[0]

        # insert genes with shifted ids by the offset
        cur.execute(
            """
            INSERT INTO main.Gene(id, contig, name, start, end, strand)
            SELECT id + ?, contig, name, start, end, strand
            FROM src.Gene
        """,
            (offset,),
        )

        # insert domains with shifted gene ids
        cur.execute(
            """
            INSERT INTO main.Domain(gene, domain, pvalue, start)
            SELECT gene + ?, domain, pvalue, start
            FROM src.Domain
        """,
            (offset,),
        )

        cur.execute("COMMIT")
        cur.execute("DETACH DATABASE src")
        cur.execute("PRAGMA foreign_keys = ON")

    @classmethod
    def merge_databases(
        cls,
        output: Path,
        databases: list[Path],
        console: Console | None = None,
    ) -> GenomeDatabase:
        """Merge a list of databases in genome schema into a new larger one"""

        if console is None:
            console = Console()  # prints to stdout. to supress pass a dummy console

        for db_path in databases:
            if not Path(db_path).exists():
                raise FileNotFoundError(f"Database {str(db_path)} not found!")

        if isinstance(output, (str, os.PathLike)) and output != ":memory:":
            if Path(output).exists():
                raise ValueError(f"File already {str(output)!r} exists!")

        # load database schema
        with cls.schema().open() as src:
            schema = src.read()

        output.parent.mkdir(parents=True, exist_ok=True)

        with rich.progress.Progress(console=console) as progress:
            db = sqlite3.connect(output)
            db.executescript(schema)

            # For speed over safety:
            db.execute("PRAGMA journal_mode = WAL")
            db.execute("PRAGMA synchronous = OFF")

            for i, next_db_path in progress.track(
                enumerate(databases),
                total=len(databases),
                description="Merging databases ...",
            ):
                console.print(f"Merging {next_db_path}")
                cls._merge_database(db=db, src_path=Path(next_db_path), first=(i == 0))

        console.print("(Re)indexing Gene Table ...")
        start = time.time()
        db.execute("CREATE INDEX IF NOT EXISTS GeneContig ON Gene(contig);")
        console.print("(Re)indexing Domain Table ...")
        db.executescript(
            """
            CREATE INDEX IF NOT EXISTS DomainGene ON Domain(gene);
            CREATE INDEX IF NOT EXISTS DomainPvalue ON Domain(pvalue);
            CREATE INDEX IF NOT EXISTS DomainStart ON Domain(start);
        """
        )
        console.print(f"(Re)indexing time: {time.time()-start}")

        db.commit()

        return cls(output, db)


class LabeledGenomeDatabase(GenomeDatabase):
    """
    Inherits from Genome Database
    """

    SCHEMA_FILE = "schema_labeled_genome.sql"

    @classmethod
    def schema(cls) -> TraversableResources:
        return importlib.resources.files(__package__).joinpath(cls.SCHEMA_FILE)

    def __init__(self, path, conn):
        super().__init__(path, conn)

    @property
    def unknown_index(self) -> int:
        return self._conn.execute(
            "SELECT id FROM Label WHERE Label.name = 'UNKNOWN'"
        ).fetchone()[0]

    @property
    def max_special_label_index(self) -> int:
        specials = [i.__repr__() for i in SpecialLabel.__members__.keys()]
        specials_str = ",".join(specials)
        query = f"""
            SELECT MAX(id)
            FROM Label
            WHERE Label.name IN ({specials_str});
        """
        return self._conn.execute(query).fetchone()[0]

    @property
    def label_vocab(self):
        query = "SELECT id, name FROM Label ORDER BY id"
        return Vocab(self._conn.execute(query).fetchall())

    @classmethod
    def merge_databases(
        cls,
        output: Path,
        databases: list[Path],
        console: Console | None = None,
    ) -> GenomeDatabase:
        """Merging is not implemented for LabeledGenomeDatabases"""

        raise NotImplementedError(
            "Merging is not implemented for LabeledGenomeDatabases"
        )

    @classmethod
    def build_labeled_database(
        cls,
        labeled_genes: LabeledGeneTable,
        pfam_features: FeatureTable,
        pfam_vocab: Vocab,
        label_vocab: Labels,
        output: str | os.PathLike,
        domain_number: int | None = None,
        only_unique_domains_per_gene: bool = True,
        n_samples: int | None = None,
        seed: int = 42,
        console: Console | None = None,
    ) -> LabeledGenomeDatabase:
        """
        Create a labeled genome sqlite3 database from a gene and domain table with a vocab.

        Arguments:
            labeled_genes: LabeledGeneTable of a genetable with labels per gene
            pfam_features: FeatureTable of Pfam feature table
            pfam_vocab: Vocab instance of Pfams
            label_vocab: Labels instance of whatever the labels are
            output:
                - str | Path: path where the db will be created, or ":memory:" for in-memory
            domain_number: Optional, `int` how many domains per gene to save
            only_unique_domains_per_gene: `bool` save only unique domains per gene
            n_samples: Optional int for how many genes to select
            seed: Seed for how many genes to randomly select
            console: Optional rich.console.Console to print progress/debug messages

        Returns:
            LabeledGenomeDatabase with an open connection.
        """
        if console is None:
            console = Console()  # prints to stdout. to supress pass a dummy console

        if domain_number is not None and only_unique_domains_per_gene:
            warnings.warn(
                "Setting 'only_unique_domains_per_gene' implies multilabel while setting 'domain_number' implies multiclass!"
            )

        if isinstance(output, (str, os.PathLike)) and output != ":memory:":
            if Path(output).exists():
                raise ValueError(f"File already {str(output)!r} exists!")

            Path(output).parent.mkdir(parents=True, exist_ok=True)

        console.print("Loading Label Vocab")
        label_vocab_table = polars.LazyFrame(
            {
                "label": label_vocab.classes_,
                "label_vocab_id": label_vocab.lookup_indices(label_vocab.classes_),
            },
            schema={
                "label": polars.Categorical("label"),
                "label_vocab_id": polars.UInt16,
            },
        )

        # sort genes and compute absolute index
        console.print("Loading Labeled Genes")
        genes_df = (
            labeled_genes.table.join(
                label_vocab_table,
                how="left",
                on="label",
                maintain_order="left",
            )
            .fill_null(label_vocab.default)
            .sort(by=["sequence_id", "start", "end"])
            .with_row_index("index", offset=1)  # sql is offset by 1
        )

        if n_samples is not None:
            # sampling without replacement
            # excluding label_vocab.default
            genes_df = genes_df.filter(
                polars.col("label_vocab_id") != label_vocab.default
            ).sample(n=n_samples, seed=seed, with_replacement=False)

        console.print("Loading Features")
        if only_unique_domains_per_gene:
            featureizer = Featureizer(pfam_features.table).unique_features_per_gene()
        else:
            featureizer = Featureizer(pfam_features.table)

        # join feature and gene table and filtering to below 1e-9
        pfam_features_df = (
            featureizer.apply_threshold_filtering(threshold_column="pvalue")
            .select_n_features(n_features=domain_number, threshold_column="pvalue")
            .annotate(genes_df)  # inner join maintaining gene order
            .features.collect()
        )

        console.print("Collecting Genes")
        genes_df = genes_df.collect()

        genome_db = super()._build_sql_db(
            genes_df=genes_df,
            features_df=pfam_features_df,
            vocab=pfam_vocab,
            console=console,
            output=output,
            index_db=True,
        )

        # select the right columns
        label_df = genes_df.select(["index", "label_vocab_id"])

        console.print(
            "Fraction of labeled Genes:",
            f"{label_df.n_unique('index')/genes_df.height}",
        )
        console.print(
            "Fraction of UNKNOWN Labels:",
            f"{label_df.filter(polars.col('label_vocab_id')==label_vocab.default).height/label_df.height}",
        )

        db = genome_db.conn  # CREATE the extra tables!

        # build database
        with rich.progress.Progress(console=console) as progress:
            # db.executescript(schema)
            db.executemany(
                "INSERT INTO LabelVocab VALUES(?, ?);",
                progress.track(
                    ((i, x) for (x, i) in label_vocab.items()),
                    total=len(label_vocab),
                    description="Label Vocab",
                ),
            )
            db.executemany(
                "INSERT INTO Label Values(?, ?);",
                progress.track(
                    label_df.iter_rows(),
                    total=label_df.height,
                    description="Labels",
                ),
            )

        start = time.time()
        console.print("Indexing Label Table ...")
        db.execute("CREATE INDEX LabelGene ON Label(gene);")
        console.print(f"Indexing time: {time.time()-start}")

        db.commit()
        return cls(output, db)
