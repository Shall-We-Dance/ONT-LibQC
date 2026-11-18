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

    def test_microc_end_queries_respect_the_contact_window(self):
        seq = "A" * 30 + "G" * 40 + "T" * 30
        with tempfile.TemporaryDirectory() as tmp:
            raw = Path(tmp) / "raw.fq"
            raw.write_text(f"@r0\n{seq}\n+\n{'I'*100}\n@raw\nTTGCAA\n+\nIIIIII\n")
            table = Path(tmp) / "table.gz"
            with gzip.open(table, "wt") as out:
                out.write("read_id\traw_length\tinsert_start\tinsert_end\tinsert_length\tis_adapter_dimer\tassignment_status\tscope_id\n"
                          "r0\t100\t0\t100\t100\t0\tassigned\tlib01\n")
            fasta = io.StringIO()
            scopes = queries([raw], table, fasta, 2, leg_length=30)
            self.assertEqual(scopes, ["lib01", RAW])
            self.assertEqual(fasta.getvalue().splitlines(), [">0", "A"*30, ">1", "A"*30, ">2", "TTGCAA"])

    def test_microc_counts_inserts_once_and_keeps_conflicting_ends_ambiguous(self):
        left = [("dm6", 80, True, ("chr2", 4))]
        right = [("dm6", 70, True, ("chr3", 1))]
        mouse = [("mm39", 68, True, ("chr1", 2))]
        tied = left + [("mm39", 77, True, None)]
        hits = [left, right, [], left, left, mouse, tied, left, [], [], mouse, []]
        result = classify(["lib01"] * 5 + [RAW], hits, ["dm6", "mm39"], "dm6", leg_length=100)
        self.assertEqual(result["rows"][0]["n"], 5)
        self.assertEqual(result["rows"][0]["counts"], {"dm6": 2, AMBIG: 2, NONE: 1})
        self.assertEqual(result["rows"][0]["windows"], {"dm6": 2})
        self.assertEqual(result["rows"][-1]["counts"], {"mm39": 1})
        self.assertEqual(result["end_evidence"], dict(concordant=1, one_end=1, discordant=1, ambiguous=1, neither=1))
        chart, details = render(result).split("<details>")
        self.assertIn("<b>40.0%</b>", chart)
        self.assertIn("up to 100 bp per end", details)

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
        self.assertEqual(details.count("class='species-counts'"), 2)
        self.assertIn("<th scope='row'>Screened</th><td>2</td><td>1</td>", details)
        self.assertIn("<td><span>1</span> <small>50.0%</small></td>", details)
        empty = classify([], [], ["hg38", "mm39"], "hg38")
        self.assertIn("<th scope='row'>Human (expected)</th><td>—</td><td>—</td>", render(empty))
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
