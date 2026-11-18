import gzip
import io
from pathlib import Path
import random
import shutil
import subprocess
import sys
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "workflow/scripts"))
from species_composition import AMBIG, NONE, RAW, TOTAL, align_hits, classify, queries, render


class SpeciesCompositionTests(unittest.TestCase):
    def test_reconstructed_units_and_raw_reads_have_separate_denominators(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            raw = root / "raw.fq"
            raw.write_text("".join(f"@r{i}\nACGTACGT\n+\nIIIIIIII\n" for i in range(3)))
            table = root / "table.gz"
            with gzip.open(table, "wt") as out:
                out.write("read_id\traw_length\tinsert_start\tinsert_end\tinsert_length\tis_adapter_dimer\tassignment_status\tscope_id\n"
                          "r0\t8\t0\t4\t4\t0\tassigned\tlib01\n"
                          "r0\t8\t4\t8\t4\t0\tunassigned\t\n"
                          "r1\t8\t4\t4\t0\t1\tunassigned\t\n")
            self.assertEqual(queries([raw], table, io.StringIO(), 3), ["lib01", "Unassigned index", RAW])
            with self.assertRaises(ValueError):
                queries([raw], table, io.StringIO(), 4)

    def test_score_margin_and_low_quality_hits(self):
        panel = ["hg38", "mm39"]
        strong = ("hg38", 60, True, ("chr1", 0))
        mixed = [strong, ("mm39", 55, True, None)]
        separated = [strong, ("mm39", 54, True, None)]
        weak = [("hg38", 90, False, None)]
        data = classify(["lib01"] * 3, [mixed, separated, weak], panel, "hg38")
        self.assertEqual(data["rows"][0]["counts"], {AMBIG: 1, "hg38": 1, NONE: 1})
        self.assertEqual(data["rows"][0]["windows"], {"hg38": 1})
        self.assertIn("Human (expected)", render(data))
        self.assertEqual(classify([], [], panel, "hg38")["rows"][0]["n"], 0)

    @unittest.skipUnless(shutil.which("bowtie2-build") and shutil.which("bowtie2"), "Bowtie 2 required")
    def test_known_mixture_including_shared_sequence(self):
        rng = random.Random(27)
        seq = lambda: "".join(rng.choices("ACGT", k=160))
        host, other, shared, novel = [seq() for _ in range(4)]
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            indexes = {}
            for species, unique in [("hg38", host), ("mm39", other)]:
                ref = root / (species + ".fa")
                ref.write_text(f">one\n{unique}\n>shared\n{shared}\n")
                indexes[species] = root / species
                subprocess.run(["bowtie2-build", "--quiet", ref, indexes[species]], check=True)
            query = root / "mix.fa"
            query.write_text("".join(f">{i}\n{s}\n" for i,s in enumerate([host, other, shared, novel, host])))
            hits = align_hits(query, 5, indexes, "bowtie2", 1)
            data = classify(["lib01", "lib02", "lib01", "lib02", RAW], hits, indexes, "hg38")
            rows = {r["scope"]: r for r in data["rows"]}
            self.assertEqual(rows[TOTAL]["counts"], {"hg38": 1, "mm39": 1, AMBIG: 1, NONE: 1})
            self.assertEqual(rows[RAW]["counts"], {"hg38": 1})
            self.assertEqual(rows["lib02"]["n"], 2)
            self.assertIn("Human (expected)", render(data))
