import io
import shutil
import tempfile
import unittest
from contextlib import ExitStack
from importlib.resources import files as resource_files
from pathlib import Path

import gb_io
import polars
import pyhmmer
import pyrodigal
from pyhmmer.easel import TextSequence
from rich.console import Console

from seqstate import Domain, Gene, Protein, Strand
from seqstate.hmm_collection import HMMSet
from seqstate.states import (FeatureTable, GeneCollection, GeneTable,
                             SequenceCollection, SequencePathCollection)
from seqstate.vocab import Labels, Vocab
from tests import test_data


class test_Strand(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        pass

    def test_sign(self):
        self.assertEqual(Strand.Coding.sign, "+")
        self.assertEqual(Strand.Reverse.sign, "-")

    def test_from_sign(self):
        self.assertEqual(Strand.from_sign("+"), 1)
        self.assertEqual(Strand.from_sign("-"), -1)
        with self.assertRaises(ValueError):
            Strand.from_sign("x")


class test_Domain(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        pass

    def test_init(self):
        domain = Domain(
            name="PF00001",
            start=0,
            end=10,
            hmm="Pfam",
            i_evalue=1e-5,
            pvalue=1e-6,
            score=42.0,
        )
        self.assertEqual(domain.name, "PF00001")
        self.assertEqual(domain.start, 0)
        self.assertEqual(domain.end, 10)
        self.assertEqual(domain.hmm, "Pfam")
        self.assertAlmostEqual(domain.i_evalue, 1e-5)
        self.assertAlmostEqual(domain.pvalue, 1e-6)
        self.assertAlmostEqual(domain.score, 42)  # ty:ignore[no-matching-overload]


class test_Protein(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        alphabet = pyhmmer.easel.Alphabet.amino()
        seq = "A" * 100
        cls.seq = pyhmmer.easel.TextSequence(name="prot1", sequence=seq).digitize(
            alphabet
        )

    def make_domain(self, start, end, pvalue) -> Domain:
        return Domain(
            name="D",
            start=start,
            end=end,
            hmm="Pfam",
            i_evalue=1e-5,
            pvalue=pvalue,
            score=None,
        )

    def test_disentangle_no_domains(self):
        d1 = self.make_domain(0, 20, 0.05)
        p = Protein("p1", self.seq, 11, (d1,))
        self.assertIs(p.disentangle(), p)

    def test_disentangle_overlapping(self):
        d1 = self.make_domain(0, 20, 0.05)
        d2 = self.make_domain(10, 30, 0.01)

        p = Protein("p1", self.seq, 11, (d1, d2))
        disentangled = p.disentangle()

        self.assertEqual(len(disentangled.domains), 1)
        self.assertEqual(disentangled.domains[0], d2)

    def test_disentangle_non_overlapping(self):
        d1 = self.make_domain(0, 10, 0.05)
        d2 = self.make_domain(20, 30, 0.01)

        p = Protein("p1", self.seq, 11, (d1, d2))
        self.assertEqual(p.disentangle(), p)

    def test_to_faa_wrapping(self):
        p = Protein("p1", self.seq, 11)
        faa = p.to_faa()

        lines = faa.splitlines()
        self.assertEqual(lines[0], ">p1")
        self.assertTrue(all(len(line) <= 80 for line in lines[1:]))
        self.assertEqual("".join(lines[1:]), "A" * 100)


class test_Gene(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.seq_text = "A" * 100
        cls.seq = pyhmmer.easel.TextSequence(
            name="prot1", sequence=cls.seq_text
        ).digitize(pyhmmer.easel.Alphabet.amino())

        cls.record = gb_io.Record(accession="ACC123", sequence=b"ATG" * 10)

        cls.protein = Protein(id="p1", seq=cls.seq, tt=11)

    def test_init(self):
        gene = Gene(
            source=self.record,
            start=0,
            end=100,
            strand=Strand.Coding,
            protein=self.protein,
        )
        self.assertEqual(gene.source, self.record)
        self.assertEqual(gene.start, 0)
        self.assertEqual(gene.end, 100)
        self.assertEqual(gene.strand, Strand.Coding)
        self.assertEqual(gene.protein, self.protein)

    def test_disentangle(self):
        d1 = Domain("D", 0, 20, "Pfam", 1e-5, 0.05, None)
        d2 = Domain("D", 10, 30, "Pfam", 1e-5, 0.01, None)

        protein = Protein("p1", self.seq, 11, (d1, d2))
        gene = Gene(self.record, 0, 100, Strand.Coding, protein)

        new_gene = gene.disentangle()
        self.assertEqual(len(new_gene.protein.domains), 1)

    def test_to_feature(self):
        # create qualifiers
        qualifiers = [
            gb_io.Qualifier(
                key="transl_table",
                value=str(11),
            ),
            gb_io.Qualifier(key="translation", value=str("A" * 100)),
            gb_io.Qualifier(
                key="locus_tag",
                value="p1",
            ),
            gb_io.Qualifier(
                key="inference",
                value=f"ab initio prediction: Pyrodigal:{pyrodigal.__version__}",
            ),
        ]
        made_feature = gb_io.Feature(
            kind="CDS",
            location=gb_io.Between(1, 100),
            qualifiers=qualifiers,
        )
        gene = Gene(
            source=self.record,
            start=0,
            end=100,
            strand=Strand.Coding,
            protein=self.protein,
        )
        # TODO qualifers need an equality comparison
        # self.assertEqual(gene.to_feature(), made_feature)

    def test_to_faa_and_gff(self):
        gene = Gene(
            source=self.record,
            start=0,
            end=100,
            strand=Strand.Coding,
            protein=self.protein,
        )
        faa, gff = gene.to_faa_and_gff()

        lines = faa.split()
        self.assertEqual(lines[0], ">p1")
        self.assertEqual(len(lines), 3)
        self.assertEqual(lines[1], "A" * 80)

        fields = gff.split("\t")
        self.assertEqual(fields[0], "ACC123")
        self.assertEqual(fields[1], f"pyrodigal_v{pyrodigal.__version__}")
        self.assertEqual(fields[2], "CDS")
        self.assertEqual(fields[3], "1")
        self.assertEqual(fields[4], "100")
        self.assertEqual(fields[5], ".")
        self.assertEqual(fields[6], "+")
        self.assertEqual(fields[7], "0")
        self.assertEqual(fields[8], "ID=p1")


class test_SequencePathCollection(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.paths1 = [
            resource_files(test_data).joinpath(
                "preprocessing", "GCF_036527415.1.contigs1+2.fna"
            ),
        ]
        cls.paths2 = [
            resource_files(test_data).joinpath(
                "preprocessing", "GCF_036527415.1.contigs1+2.fna"
            ),
            resource_files(test_data).joinpath(
                "preprocessing", "GCF_036527415.1.contigs3_with_features.gbff.gz"
            ),
        ]
        cls.paths3 = [
            resource_files(test_data).joinpath(
                "preprocessing", "GCF_036527415.1.contigs1+2.fna"
            ),
            resource_files(test_data).joinpath(
                "preprocessing", "GCF_036527415.1.contigs1+2.fna"
            ),
        ]
        cls.paths4 = [
            resource_files(test_data).joinpath(
                "preprocessing", "GCF_000013265.1_part.faa"
            ),
            resource_files(test_data).joinpath(
                "preprocessing", "GCF_000013265.1_part.gff"
            ),
        ]
        cls.console = Console(file=io.StringIO())  # silent console

    def test_format_guessing(self):
        self.assertEqual(
            SequencePathCollection.guess_sequences_format(
                self.paths2[0]  # ty:ignore[invalid-argument-type]
            ),
            "fna",
        )
        self.assertEqual(
            SequencePathCollection.guess_sequences_format(
                self.paths2[1]  # ty:ignore[invalid-argument-type]
            ),
            "genbank",
        )
        self.assertEqual(
            SequencePathCollection.guess_sequences_format(
                self.paths4[0]  # ty:ignore[invalid-argument-type]
            ),
            "faa",
        )
        self.assertEqual(
            SequencePathCollection.guess_sequences_format(
                self.paths4[1]  # ty:ignore[invalid-argument-type]
            ),
            "gff",
        )

    def test_init(self):
        paths = self.paths3
        col1 = SequencePathCollection(
            paths=paths,  # ty:ignore[invalid-argument-type]
            console=self.console,
        )
        self.assertEqual(col1[0], paths[0])
        self.assertEqual(len(col1), 2)

    def test_fna(self):
        col = SequencePathCollection(
            paths=self.paths1, console=self.console  # ty:ignore[invalid-argument-type]
        ).load_files()
        self.assertEqual(len(list(col)), 2)

    def test_fna_with_genbank(self):
        col = SequencePathCollection(
            paths=self.paths2, console=self.console  # ty:ignore[invalid-argument-type]
        ).load_files()
        self.assertEqual(len(list(col)), 3)

    def test_gff(self):
        col = SequencePathCollection(
            format="faa+gff",
            paths=self.paths4,  # ty:ignore[invalid-argument-type]
            console=self.console,
        ).load_files()
        self.assertEqual(len(list(col)), 1)

        col2 = SequencePathCollection(
            paths=self.paths4,  # ty:ignore[invalid-argument-type]
            console=self.console,
        ).load_files()
        self.assertEqual(len(list(col2)), 1)

    # def test_only_faa(self):
    #     col = SequencePathCollection(
    #         paths=[self.paths4[0]],  # ty:ignore[invalid-argument-type]
    #         console=self.console,
    #     ).load_files()

    #     col.split_and_dump(Path("/home/ray/Documents/seqstate/file.faa"), compressed=False)

    #     x = col.extract_genes()
    #     print(x)

    def test_errors(self):
        with self.assertRaises(ValueError):
            SequencePathCollection(
                paths=self.paths3,  # ty:ignore[invalid-argument-type]
                console=self.console,
            ).load_files()

        with self.assertRaises(FileNotFoundError):
            paths = [Path("path1"), Path("path2")]
            SequencePathCollection(paths=paths, console=self.console).load_files()

    def test_min_record_length(self):
        pass


class test_SequenceCollection(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.record = gb_io.Record(
            sequence=b"ATCGTGGT",
            name="seq1",
            length=8,
            molecule_type="DNA",
            accession="acc1",
        )
        cls.record2 = gb_io.Record(
            sequence=b"GGGATTGtttaa",
            name="seq2",
            length=12,
            molecule_type="DNA",
            accession="acc2",
        )
        cls.console = Console(file=io.StringIO())  # silent console

    def test_init(self):
        col1 = SequenceCollection(
            sequences=[self.record, self.record2],
            console=self.console,
        )
        self.assertEqual(col1.sequences, [self.record, self.record2])
        self.assertEqual(col1.console, self.console)

    def test_iter(self):
        col1 = SequenceCollection(sequences=[self.record, self.record2])
        self.assertEqual(list(col1), [self.record, self.record2])

    # def test_split_and_dump_roundtrip(self):
    #     # TODO requires equality comparison for gb_io.Record objects
    #     col1 = SequenceCollection(sequences=[self.record, self.record2])
    #     with tempfile.TemporaryDirectory() as tmpdir:
    #         out_path = Path(tmpdir, "output.gbk")
    #         col1.split_and_dump(output=out_path, compressed=False, batch_size=None)
    #         new = SequencePathCollection(
    #             paths=[out_path],
    #             console=self.console,
    #         ).load_files()
    #         self.assertEqual(list(col1.sequences), list(new.sequences))


class test_GeneCollection(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.paths1 = [
            resource_files(test_data).joinpath(
                "preprocessing", "GCF_036527415.1.contigs1+2.fna"
            ),
        ]
        cls.paths2 = [
            resource_files(test_data).joinpath(
                "preprocessing", "GCF_036527415.1.contigs3_with_features.gbff.gz"
            ),
        ]

        cls.console = Console(file=io.StringIO())  # silent console

        # create a temp dir
        cls.tempdir = tempfile.TemporaryDirectory(delete=False)

        # create a temp file inside that dir
        cls.tempfile = tempfile.NamedTemporaryFile(
            dir=cls.tempdir.name,
            delete=False,
            mode="w",
            suffix=".gbk",
        )

    @classmethod
    def tearDownClass(cls):
        # Cleanup the temporary directory and all contents
        cls.tempfile.close()
        shutil.rmtree(Path(tempfile.gettempdir(), cls.tempdir.name))

    def test_rounttrip_via_featured_seq_collection(self):
        # TODO!
        # warnings because of invalid CDS sequences which are not multiples of 3
        # with self.assertWarns(UserWarning):
        gene_collection = (
            SequencePathCollection(
                self.paths2,  # ty:ignore[invalid-argument-type]
                console=self.console,
            )
            .load_files()
            .extract_genes(
                call_genes=False,
                cds_feature="CDS",
                streaming=False,
            )
        )
        featured_seq_col = gene_collection.to_featured_sequence_collection(
            is_ordered=True
        )
        # roundtrip via featured sequence collection
        new = featured_seq_col.extract_genes(
            call_genes=False,
            cds_feature="CDS",
            streaming=False,
        )
        self.assertEqual(list(gene_collection.genes), list(new.genes))

    def test_roundtrip_via_io_faa(self):
        gene_collection = (
            SequencePathCollection(
                self.paths1,  # ty:ignore[invalid-argument-type]
                console=self.console,
            )
            .load_files()
            .extract_genes(
                call_genes=True,
                streaming=False,
            )
        )
        featured_seq_col = gene_collection.to_featured_sequence_collection(
            is_ordered=True
        )

        # roundtrip via split_and_dump
        featured_seq_col.split_and_dump(
            output=Path(self.tempfile.name), compressed=False, batch_size=None
        )
        new = (
            SequencePathCollection(
                paths=[Path(self.tempfile.name)],
                console=self.console,
            )
            .load_files(streaming=False)
            .extract_genes(call_genes=False, streaming=False)
        )
        self.assertEqual(list(gene_collection.genes), list(new.genes))

        # roundtrip via dump
        outdir = output = Path(self.tempdir.name)
        output = Path(outdir, "file")
        gene_collection.dump(
            output=output,
        )
        new = (
            SequencePathCollection(
                paths=[Path(outdir, "file.gff"), Path(outdir, "file.faa")],
                console=self.console,
                format="faa+gff",
            )
            .load_files(streaming=False)
            .extract_genes(
                call_genes=False,
                cds_feature="CDS",
                streaming=False,
            )
        )
        self.assertEqual(list(gene_collection.genes), list(new.genes))

    def test_roundtrip_via_io_gbk(self):
        gene_collection = (
            SequencePathCollection(
                self.paths2,  # ty:ignore[invalid-argument-type]
                console=self.console,
            )
            .load_files()
            .extract_genes(
                call_genes=True,
                streaming=False,
            )
        )
        # roundtrip via dump
        outdir = output = Path(self.tempdir.name)
        output = Path(outdir, "file")
        gene_collection.dump(
            output=output,
        )
        new = (
            SequencePathCollection(
                paths=[Path(outdir, "file.gff"), Path(outdir, "file.faa")],
                console=self.console,
                format="faa+gff",
            )
            .load_files(streaming=False)
            .extract_genes(
                call_genes=False,
                cds_feature="CDS",
                streaming=False,
            )
        )
        self.assertEqual(list(gene_collection.genes), list(new.genes))


class test_AnnotatedGeneCollection(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.paths1 = [
            resource_files(test_data).joinpath(
                "preprocessing", "GCF_036527415.1.contigs1+2.fna"
            ),
        ]
        cls.console = Console(file=io.StringIO())  # silent console
        cls.gene_col = (
            SequencePathCollection(
                cls.paths1,  # ty:ignore[invalid-argument-type]
                console=cls.console,
            )
            .load_files()
            .extract_genes(
                call_genes=True,
                streaming=False,
            )
        )
        cls.hmm = HMMSet(
            id="test",
            version="0.1",
            url="Pfam37.1.first3",
            path=resource_files(test_data).joinpath(
                "preprocessing", "selected_hmms.h3m"
            ),  # ty:ignore[invalid-argument-type]
            md5="ea1bf925ae8419d89b10644c8cb79395",
            relabel_with=r"s/(PF\d+)\.\d+/\1/",
            size=3,
        )
        # create a temp dir
        cls.tempdir = tempfile.TemporaryDirectory(delete=False)

        # create a temp file inside that dir
        cls.temp_tsv_file = tempfile.NamedTemporaryFile(
            dir=cls.tempdir.name,
            delete=False,
            mode="w",
            suffix=".tsv",
        )
        cls.temp_tsv_file2 = tempfile.NamedTemporaryFile(
            dir=cls.tempdir.name,
            delete=False,
            mode="w",
            suffix=".tsv",
        )

    @classmethod
    def tearDownClass(cls):
        # Cleanup the temporary directory and all contents
        cls.temp_tsv_file.close()
        cls.temp_tsv_file2.close()
        shutil.rmtree(Path(tempfile.gettempdir(), cls.tempdir.name))

    def test_genetable(self):
        self.gene_col.create_gene_table().dump(self.temp_tsv_file.name, format="tsv")

        with ExitStack() as ctx:
            f = ctx.enter_context(
                open(
                    resource_files(test_data).joinpath(
                        "preprocessing", "genetable.tsv"
                    ),
                    "r",
                )  # ty:ignore[no-matching-overload]
            )
            f2 = ctx.enter_context(open(self.temp_tsv_file.name, "r"))

            content, content2 = f.read(), f2.read()
            self.assertMultiLineEqual(content, content2)

    def test_genetable_load(self):
        tbl = GeneTable.load(self.temp_tsv_file.name)
        expected = GeneTable.load(
            resource_files(test_data).joinpath(
                "preprocessing", "genetable.tsv"
            ),  # ty:ignore[invalid-argument-type]
        ).table.collect()

        self.assertTrue(tbl.table.collect().equals(expected))

    def test_annotation(self):
        annot = self.gene_col.annotate_domains(
            hmms=[self.hmm],
            pvalue_filter=1e-9,
        )

        self.assertEqual(len(annot), 414)

        annot.create_feature_table().dump(self.temp_tsv_file2.name, format="tsv")
        with ExitStack() as ctx:
            f = ctx.enter_context(
                open(
                    resource_files(test_data).joinpath(
                        "preprocessing", "featuretable.tsv"
                    ),
                    "r",
                )  # ty:ignore[no-matching-overload]
            )
            f2 = ctx.enter_context(open(self.temp_tsv_file2.name, "r"))

            content, content2 = f.read(), f2.read()
            self.assertMultiLineEqual(content, content2)

    def test_annotation_with_multidomains(self):
        seq = """MRVLKFGGTSVANAERFLRVADILESNARQGQVATVLSAPAKITNHLVAMIEKTISGQDALPNISDAERIFAELLTGLAA
            AQPGFPLAQLKTFVDQEFAQIKHVLHGISLLGQCPDSINAALICRGEKMSIAIMAGVLEARGHNVTVIDPVEKLLAVGHY
            LESTVDIAESTRRIAASRIPADHMVLMAGFTAGNEKGELVVLGRNGSDYSAAVLAACLRADCCEIWTDVDGVYTCDPRQV
            PDARLLKSMSYQEAMELSYFGAKVLHPRTITPIAQFQIPCLIKNTGNPQAPGTLIGASRDEDELPVKGISNLNNMAMFSV
            SGPGMKGMVGMAARVFAAMSRARISVVLITQSSSEYSISFCVPQSDCVRAERAMQEEFYLELKEGLLEPLAVTERLAIIS
            VVGDGMRTLRGISAKFFAALARANINIVAIAQGSSERSISVVVNNDDATTGVRVTHQMLFNTDQVIEVFVIGVGGVGGAL
            LEQLKRQQSWLKNKHIDLRVCGVANSKALLTNVHGLNLENWQEELAQAKEPFNLGRLIRLVKEYHLLNPVIVDCTSSQAV
            ADQYADFLREGFHVVTPNKKANTSSMDYYHQLRYAAEKSRRKFLYDTNVGAGLPVIENLQNLLNAGDELMKFSGILSGSL
            SYIFGKLDEGMSFSEATTLAREMGYTEPDPRDDLSGMDVARKLLILARETGRELELADIEIEPVLPAEFNAEGDVAAFMA
            NLSQLDDLFAARVAKARDEGKVLRYVGNIDEDGVCRVKIAEVDGNDPLFKVKNGENALAFYSHYYQPLPLVLRGYGAGND
            VTAAGVFADLLRTLSWKLGV
        """.replace(
            "\n", ""
        ).replace(
            " ", ""
        )
        protein = Protein(
            id="ecoli_gene_1",
            seq=TextSequence(sequence=seq).digitize(pyhmmer.easel.Alphabet.amino()),
            tt=11,
        )
        gene_col = GeneCollection(
            genes=[
                Gene(
                    source=gb_io.Record(sequence=b"ATCG", accession="ecoli_gene_1"),
                    start=1,
                    end=10,
                    strand=Strand(1),
                    protein=protein,
                ),
            ],
            console=self.console,
        )
        annot = gene_col.annotate_domains(
            # TODO replace with proper path downloaded during setup
            hmms=[
                HMMSet(
                    "Pfam",
                    version="37.1",
                    size=23794,
                    url=None,
                    path=resource_files(test_data).joinpath(
                        "preprocessing",
                        "6hmms.hmm.zst",  # ty:ignore[too-many-positional-arguments]
                    ),  # ty:ignore[invalid-argument-type]
                    relabel_with=r"s/(PF\d+)\.\d+/\1/",
                )
            ],
            pvalue_filter=1e-9,
        )
        all = annot.create_feature_table().table.collect()
        self.assertEqual(all.height, 8)

        tangle = annot.disentangle_domains().create_feature_table().table.collect()
        self.assertEqual(tangle.height, 5)

        self.assertEqual(
            tangle["domain"].to_list(),
            ["PF00696", "PF22468", "PF22468", "PF03447", "PF00742"],
        )

    def test_whitelist(self):
        whitelist = sorted(["PF00535", "PF16640"])  # "PF00037"
        annot = self.gene_col.annotate_domains(
            hmms=[self.hmm],
            pvalue_filter=1e-11,
            whitelist=whitelist,
        )
        table = annot.create_feature_table().table.collect()
        self.assertEqual(
            sorted(table["domain"].unique().to_list()),
            whitelist,
        )
        self.assertEqual(len(table.filter(polars.col("pvalue") > 1e-11)), 0)

    def test_featuretable_load(self):
        tbl = FeatureTable.load(self.temp_tsv_file2.name)
        expected = FeatureTable.load(
            resource_files(test_data).joinpath(
                "preprocessing", "featuretable.tsv"
            ),  # ty:ignore[invalid-argument-type]
        ).table.collect()

        self.assertTrue(tbl.table.collect().equals(expected))

    def test_build_vocab_from_featuretable(self):
        tbl = FeatureTable.load(self.temp_tsv_file2.name)
        vocab, _ = tbl.build_vocab(console=self.console)
        self.assertEqual(vocab, Vocab.from_iterator(["PF00535", "PF16640", "PF00037"]))

        vocab, _ = tbl.build_vocab(console=self.console, min_domain_occurence=7)
        self.assertEqual(vocab, Vocab.from_iterator(["PF00535"]))

        genes = GeneTable.load(
            resource_files(test_data).joinpath(
                "preprocessing", "genetable.tsv"
            ),  # ty:ignore[invalid-argument-type]
        )
        vocab, _ = tbl.build_vocab(console=self.console, genes=genes)
        self.assertEqual(vocab, Vocab.from_iterator(["PF00535", "PF16640", "PF00037"]))

    def test_build_labels(self):
        genes = GeneTable.load(
            resource_files(test_data).joinpath(
                "preprocessing", "genetable.tsv"
            ),  # ty:ignore[invalid-argument-type]
        )
        # boolean but not 1:1 mapping
        label_df = polars.LazyFrame(
            {
                "protein_id": [
                    "NZ_JAQOLM010000012.1_1",
                    "NZ_JAQOLM010000012.1_1",
                    "NZ_JAQOLM010000012.1_2",
                    "NZ_JAQOLM010000012.1_3",
                    "NZ_JAQOLM010000012.1_3",
                ],
                "test_label": [str(i) for i in [True, False, True, True, False]],
            },
            schema={
                "protein_id": polars.Utf8,
                "test_label": polars.Categorical("label"),
            },
        ).with_columns(
            sequence_id=polars.lit("NZ_JAQOLM010000012.1").cast(
                polars.Categorical("sequence_id")
            )
        )

        with self.assertRaises(polars.exceptions.ComputeError):
            lbl_tbl = genes.label(label_df, label_name="test_label")
            labels = lbl_tbl.build_label_vocab()

        # classic string multilabel case
        label_df = polars.LazyFrame(
            {
                "protein_id": [
                    "NZ_JAQOLM010000012.1_1",
                    "NZ_JAQOLM010000012.1_1",
                    "NZ_JAQOLM010000012.1_2",
                    "NZ_JAQOLM010000012.1_3",
                    "NZ_JAQOLM010000012.1_3",
                ],
                "test_label": ["cls1", "cls2", "cls1", "cls3", "cls4"],
            },
            schema={
                "protein_id": polars.Utf8,
                "test_label": polars.Categorical("label"),
            },
        ).with_columns(
            sequence_id=polars.lit("NZ_JAQOLM010000012.1").cast(
                polars.Categorical("sequence_id")
            )
        )
        lbl_tbl = genes.label(label_df, label_name="test_label")
        labels = lbl_tbl.build_label_vocab()
        self.assertEqual(labels, Labels.from_iterator(["cls1", "cls2", "cls3", "cls4"]))


if __name__ == "__main__":
    unittest.main()
