import unittest

import polars

from seqstate.featureizer import Featureizer


class test_Featureizer(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        pass

    def test_domain_disentangling(self):

        test_df = polars.DataFrame(
            [
                # no overlaps, all three should survive
                {
                    "sequence_id": "A",
                    "protein_id": "P1",
                    "domain_start": 10,
                    "domain_end": 50,
                    "pvalue": 0.01,
                },
                {
                    "sequence_id": "A",
                    "protein_id": "P1",
                    "domain_start": 100,
                    "domain_end": 150,
                    "pvalue": 0.05,
                },
                {
                    "sequence_id": "A",
                    "protein_id": "P1",
                    "domain_start": 200,
                    "domain_end": 250,
                    "pvalue": 0.001,
                },
                # three-way transitive overlap, only best pvalue (start=40) survives
                {
                    "sequence_id": "A",
                    "protein_id": "P2",
                    "domain_start": 10,
                    "domain_end": 50,
                    "pvalue": 0.05,
                },
                {
                    "sequence_id": "A",
                    "protein_id": "P2",
                    "domain_start": 40,
                    "domain_end": 90,
                    "pvalue": 0.001,
                },  # best
                {
                    "sequence_id": "A",
                    "protein_id": "P2",
                    "domain_start": 80,
                    "domain_end": 130,
                    "pvalue": 0.02,
                },
                # encompassing overlap, only best pvalue (start=50) survives
                {
                    "sequence_id": "A",
                    "protein_id": "P3",
                    "domain_start": 10,
                    "domain_end": 200,
                    "pvalue": 0.5,
                },
                {
                    "sequence_id": "A",
                    "protein_id": "P3",
                    "domain_start": 50,
                    "domain_end": 100,
                    "pvalue": 0.001,
                },  # best
                {
                    "sequence_id": "A",
                    "protein_id": "P3",
                    "domain_start": 80,
                    "domain_end": 150,
                    "pvalue": 0.01,
                },
                # encompassing overlap, only best pvalue (start=50) survives
                # NOTE this could be a bit controversial!
                {
                    "sequence_id": "A",
                    "protein_id": "P4",
                    "domain_start": 10,
                    "domain_end": 200,
                    "pvalue": 0.5,
                },
                {
                    "sequence_id": "A",
                    "protein_id": "P4",
                    "domain_start": 50,
                    "domain_end": 79,
                    "pvalue": 0.001,
                },  # best
                {
                    "sequence_id": "A",
                    "protein_id": "P4",
                    "domain_start": 80,
                    "domain_end": 150,
                    "pvalue": 0.01,
                },
                # encompassing overlap, only best pvalue (start=50) survives
                {
                    "sequence_id": "A",
                    "protein_id": "P5",
                    "domain_start": 10,
                    "domain_end": 200,
                    "pvalue": 0.001,
                },  # best
                {
                    "sequence_id": "A",
                    "protein_id": "P5",
                    "domain_start": 50,
                    "domain_end": 100,
                    "pvalue": 0.005,
                },
                {
                    "sequence_id": "A",
                    "protein_id": "P5",
                    "domain_start": 80,
                    "domain_end": 150,
                    "pvalue": 0.01,
                },
            ]
        )

        ft = Featureizer(features=test_df).disentangle_domains()

        expected_df = polars.DataFrame(
            [
                {
                    "sequence_id": "A",
                    "protein_id": "P1",
                    "domain_start": 10,
                    "domain_end": 50,
                    "pvalue": 0.01,
                },
                {
                    "sequence_id": "A",
                    "protein_id": "P1",
                    "domain_start": 100,
                    "domain_end": 150,
                    "pvalue": 0.05,
                },
                {
                    "sequence_id": "A",
                    "protein_id": "P1",
                    "domain_start": 200,
                    "domain_end": 250,
                    "pvalue": 0.001,
                },
                {
                    "sequence_id": "A",
                    "protein_id": "P2",
                    "domain_start": 40,
                    "domain_end": 90,
                    "pvalue": 0.001,
                },  # best
                {
                    "sequence_id": "A",
                    "protein_id": "P3",
                    "domain_start": 50,
                    "domain_end": 100,
                    "pvalue": 0.001,
                },  # best
                {
                    "sequence_id": "A",
                    "protein_id": "P4",
                    "domain_start": 50,
                    "domain_end": 79,
                    "pvalue": 0.001,
                },  # best
                {
                    "sequence_id": "A",
                    "protein_id": "P5",
                    "domain_start": 10,
                    "domain_end": 200,
                    "pvalue": 0.001,
                },  # best
            ]
        )

        self.assertTrue(expected_df.equals(ft.features.sort("protein_id")))
