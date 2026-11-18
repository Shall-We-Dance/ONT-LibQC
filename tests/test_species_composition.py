import gzip
import json
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

    def test_visible_bars_use_each_assigned_library_as_denominator(self):
        hit = ("hg38", 60, True, ("chr1", 0))
        scopes = ["lib01", "lib01", "lib02", "Unassigned index", "Ambiguous index", RAW]
        data = classify(scopes, [[hit], [], [hit], [], [], []], ["hg38", "mm39"], "hg38")
        chart, details = render(data).split("<details>")
        self.assertEqual(chart.count("class='species-bar'"), 2)
        self.assertIn("lib01, Human (expected) 50.0%", chart)
        self.assertIn("lib02, Human (expected) 100.0%", chart)
        for scope in (TOTAL, RAW, "Unassigned index", "Ambiguous index"):
            self.assertNotIn(scope, chart)
            self.assertIn(scope, details)
        empty = classify([], [], ["hg38", "mm39"], "hg38")
        self.assertIn("No inserts assigned to a sublibrary", render(empty))

    def test_refresh_moves_existing_section_and_navigation_once(self):
        from report_refresh import refresh_species_layout
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            directory = root / "data/species_composition"
            directory.mkdir(parents=True)
            data = classify(["lib01"], [[]], ["hg38", "mm39"], "hg38")
            (directory / "sample.json").write_text(json.dumps(data))
            sequencing = "<section id='sequencing'></section>"
            mixture = "<section id='sublibrary-composition'><h2>Sublibrary composition</h2></section>"
            nav = "<nav><a href='#sequencing'>Sequencing</a><a href='#species-composition'>Species</a></nav>"
            for previous in (mixture, ""):
                page = nav + previous + sequencing + render(data)
                new = refresh_species_layout(root, "sample", page)
                self.assertIn(previous + render(data) + sequencing, new)
                self.assertIn("<nav><a href='#species-composition'>", new)
                self.assertEqual(refresh_species_layout(root, "sample", new), new)

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
