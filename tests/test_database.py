import io
import shutil
import tempfile
import unittest
from pathlib import Path

try:
    from importlib.resources import files as resource_files
except ImportError:
    from importlib_resources import files as resource_files  # type: ignore

from rich.console import Console

from seqstate.database import GenomeDatabase
from seqstate.tables import FeatureTable, GeneTable
from seqstate.vocab import Symbol
from tests import test_data


class test_GenomeDatabase(unittest.TestCase):
    @classmethod
    def setUpClass(cls):

        cls.genes = GeneTable.load(
            resource_files(test_data).joinpath(
                "preprocessing", "genetable.tsv"
            ),  # ty:ignore[invalid-argument-type]
        )
        cls.features = FeatureTable.load(
            resource_files(test_data).joinpath(
                "preprocessing", "featuretable.tsv"
            ),  # ty:ignore[invalid-argument-type]
        )

        cls.console = Console(file=io.StringIO())  # silent console
        # create a temp dir
        cls.tempdir = tempfile.TemporaryDirectory(delete=False)

        # # create a temp file inside that dir
        # cls.temp_db = tempfile.NamedTemporaryFile(
        #     dir=cls.tempdir.name,
        #     delete=True,
        #     mode="w",
        #     suffix=".sqlite3",
        # )

    @classmethod
    def tearDownClass(cls):
        # Cleanup the temporary directory and all contents
        shutil.rmtree(Path(tempfile.gettempdir(), cls.tempdir.name))

    def test_build_db(self):
        vocab, _ = self.features.build_vocab(
            console=self.console,
            genes=self.genes,
            only_unique_domains_per_gene=False,
            min_domain_occurence=11,  # has no effect here
        )
        db = GenomeDatabase.build_database(
            genes=self.genes,
            features=self.features,
            vocab=vocab,
            only_unique_domains_per_gene=False,
            output=Path(self.tempdir.name, "db.sqlite3"),
            console=self.console,
        )

        expected_domain_occurences = {
            4: 13,
            5: 13,
            6: 23,
        }

        self.assertEqual(len(db), self.genes.table.collect().height)
        self.assertEqual(db.vocab_size, len(vocab))
        self.assertEqual(db.vocab, vocab)
        self.assertEqual(db.num_annotated_genes, 25)
        self.assertEqual(db.padding_index, 0)
        self.assertEqual(db.max_special_index, 3)
        self.assertEqual(db.domain_occurences, expected_domain_occurences)
        self.assertEqual(len(db.domain_occurences) + len(Symbol), db.vocab_size)

        domain_len = db.conn.execute("SELECT COUNT(*) FROM Domain").fetchone()[0]
        self.assertEqual(domain_len, 49)

        db.close()

    def test_build_db_after_vocab(self):
        vocab, annotated_features = self.features.build_vocab(
            console=self.console,
            genes=self.genes,
            min_domain_occurence=5,  # below all of them
            only_unique_domains_per_gene=True,
        )

        db = GenomeDatabase._build_sql_db(
            genes_df=self.genes.table.collect(),
            features_df=annotated_features,
            vocab=vocab,
            console=self.console,
            output=Path(self.tempdir.name, "multilabel_db.sqlite3"),
            index_db=True,
        )
        db.close()

        self.assertEqual(len(vocab), 7)

        db = GenomeDatabase.open(Path(self.tempdir.name, "multilabel_db.sqlite3"))

        expected_domain_occurences = {
            4: 6,
            5: 13,
            6: 6,
        }
        self.assertEqual(len(db), self.genes.table.collect().height)
        self.assertEqual(db.vocab_size, len(vocab))
        self.assertEqual(db.vocab, vocab)
        self.assertEqual(db.num_annotated_genes, 25)
        self.assertEqual(
            db.num_annotated_genes,
            self.features.table.unique("protein_id").collect().height,
        )
        self.assertEqual(db.padding_index, 0)
        self.assertEqual(db.max_special_index, 3)
        self.assertEqual(db.domain_occurences, expected_domain_occurences)
        self.assertEqual(len(db.domain_occurences) + len(Symbol), db.vocab_size)

        # Test for the length of the domain table!
        domain_len = db.conn.execute("SELECT COUNT(*) FROM Domain").fetchone()[0]
        self.assertEqual(domain_len, db.num_annotated_genes)
        db.close()
