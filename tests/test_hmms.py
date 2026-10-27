import unittest
from importlib.resources import files as resource_files
from pathlib import Path

from seqstate.hmm_collection import HMMSet
from tests import test_data


class test_HMM(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        pass

    def test_h3m_input(self):
        hmm = HMMSet(
            id="test",
            version="0.1",
            url="Pfam37.1.first3",
            path=resource_files(test_data).joinpath(
                "preprocessing/selected_hmms.h3m",
            ),  # ty:ignore[invalid-argument-type]
            md5="ea1bf925ae8419d89b10644c8cb79395",
            size=3,
            relabel_with=r"s/(PF\d+)\.\d+/\1/",
        )
        self.assertEqual(hmm.id, "test")
        self.assertEqual(hmm.version, "0.1")
        self.assertEqual(hmm.url, "Pfam37.1.first3")
        self.assertEqual(hmm.md5, HMMSet._compute_md5(Path(hmm.path)))
        self.assertEqual(hmm.relabel("PF00001.1"), "PF00001")

        with self.assertRaises(ValueError):
            HMMSet(
                id="test",
                version="0.1",
                url="Pfam37.1.first3",
                path=resource_files(test_data).joinpath(
                    "preprocessing/selected_hmms.h3m",
                ),  # ty:ignore[invalid-argument-type]
                md5="fa1bf925ae8419d89b10644c8cb79395",
                relabel_with=r"s/(PF\d+)\.\d+/\1/",
            )
        with self.assertRaises(FileNotFoundError):
            HMMSet(
                id="test",
                version="0.1",
                url="Pfam37.1.first3",
                path=resource_files(test_data).joinpath(
                    "preprocessing/nofile.hmm.h3m",
                ),  # ty:ignore[invalid-argument-type]
                relabel_with=r"s/(PF\d+)\.\d+/\1/",
            )

    def test_zst_input(self):
        hmm = HMMSet(
            id="test",
            version="0.1",
            url=None,
            path=resource_files(test_data).joinpath(
                "preprocessing/6hmms.hmm.zst",
            ),  # ty:ignore[invalid-argument-type]
            md5="72ee2075eb3285f844ff90306f318905",
            size=6,
            relabel_with=r"s/(PF\d+)\.\d+/\1/",
        )
        self.assertEqual(
            hmm.md5, HMMSet._compute_md5(Path(hmm.path))
        ), HMMSet._compute_md5(Path(hmm.path))
