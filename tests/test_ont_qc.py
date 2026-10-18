import gzip
import importlib.util
import math
import tempfile
import unittest
from unittest import mock
from pathlib import Path


SCRIPT = Path(__file__).parents[1] / "workflow" / "scripts" / "ont_qc.py"
SPEC = importlib.util.spec_from_file_location("ont_qc", SCRIPT)
ont_qc = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(ont_qc)


class ScopeProfileTests(unittest.TestCase):
    def test_raw_q_are_per_molecule_and_insert_lengths_are_per_unit(self):
        header = ["read_id", "raw_length", "mean_q", "unit_id", "scope_id",
                  "assignment_status", "insert_length"]
        rows = [
            ["r1", 500, 30, "r1|u1", "lib01", "assigned", 100],
            ["r1", 500, 30, "r1|u2", "lib01", "assigned", 120],
            ["r2", 700, 34, "r2|u1", "lib01", "assigned", 200],
            ["r3", 900, 20, "r3|u1", "", "unassigned", 300],
        ]
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "reads.tsv.gz"
            with gzip.open(path, "wt") as handle:
                handle.write("\t".join(header) + "\n")
                for row in rows:
                    handle.write("\t".join(map(str, row)) + "\n")
            profile = ont_qc.read_scope_profiles(path)["lib01"]
        self.assertEqual(sorted(profile["raw_lengths"]), [500, 700])
        self.assertEqual(profile["insert_lengths"], [100, 120, 200])
        self.assertEqual(profile["mean_read_q"], 32)


class CollisionOETests(unittest.TestCase):
    def test_clustered_reads_have_larger_collision_oe(self):
        sizes = {"chr1": 100000}
        spread = [{"chrom": "chr1", "start": i * 10000, "end": i * 10000 + 100} for i in range(10)]
        clustered = [{"chrom": "chr1", "start": 100, "end": 200} for _ in range(10)]
        self.assertGreater(ont_qc.bin_stats(clustered, sizes, 10000)["pair_oe"],
                           ont_qc.bin_stats(spread, sizes, 10000)["pair_oe"])


class IndexPairClusteringTests(unittest.TestCase):

    def units(self, pair_counts):
        units = []
        number = 0
        for (i5, i7), count in pair_counts:
            for _ in range(count):
                number += 1
                units.append({"unit_id": f"r{number}|u1", "read_id": f"r{number}",
                              "i5": i5, "i7": i7, "p5_detected": True,
                              "p7_detected": True, "is_adapter_dimer": False})
        return units

    def test_one_edit_ont_variants_merge_into_dominant_pair(self):
        units = self.units([(("CGAGAGAA", "GGAAGAGA"), 100),
                            (("CAGAGAA", "GGAAGAGA"), 12),
                            (("CAGAGAGAA", "GGAAGAGA"), 4)])
        result, _ = ont_qc.cluster_index_pairs(
            units, minimum_cluster_units=2, minimum_cluster_fraction=.01)
        self.assertEqual(result["mode"], "single_library")
        self.assertEqual(len(result["pairs"]), 1)
        self.assertEqual(result["pairs"][0]["i5"], "CGAGAGAA")
        self.assertEqual(result["assigned_informative_units"], 116)
        self.assertEqual(result["top_pair_fraction"], 1.0)

    def test_two_bp_i5_indels_merge_only_into_dominant_pair(self):
        # HL238-like ONT errors: CGTATC is a two-base deletion and
        # CGTATCTCAC is a two-base insertion relative to CGTATCTC.
        units = self.units([(("CGTATCTC", "CAAGGTAC"), 1803),
                            (("CGTATC", "CAAGGTAC"), 100),
                            (("CGTATCTCAC", "CAAGGTAC"), 88),
                            (("CGTATCT", "CAAGGTAC"), 35)])
        result, assignments = ont_qc.cluster_index_pairs(units)
        self.assertEqual(len(result["pairs"]), 1)
        self.assertEqual(result["pairs"][0]["i5"], "CGTATCTC")
        self.assertEqual(result["assigned_informative_units"], len(units))
        self.assertTrue(all(row["scope_id"] == "lib01" for row in assignments))

    def test_similarly_supported_two_edit_pairs_remain_multiplexed(self):
        units = self.units([(("CGAGAGAA", "GGAAGAGA"), 60),
                            (("CGAGAA", "GGAAGAGA"), 50)])
        result, _ = ont_qc.cluster_index_pairs(
            units, minimum_cluster_units=2, minimum_cluster_fraction=.01)
        self.assertEqual(result["mode"], "multiplexed_library")
        self.assertEqual(len(result["pairs"]), 2)

    def test_distant_pairs_remain_multiplexed(self):
        units = self.units([(("CGAGAGAA", "GGAAGAGA"), 60),
                            (("TCTCTCTC", "ACACACAC"), 40)])
        result, _ = ont_qc.cluster_index_pairs(
            units, minimum_cluster_units=2, minimum_cluster_fraction=.01)
        self.assertEqual(result["mode"], "multiplexed_library")
        self.assertEqual(len(result["pairs"]), 2)


class KaryotypeTests(unittest.TestCase):
    def test_karyoploter_receives_only_mapq_pass_primary_chromosome_midpoints(self):
        sizes = {"chr1": 1000, "chr2": 500, "chrM": 100}
        records = {
            "r1": {"chrom": "chr1", "start": 0, "end": 10, "mapq": 20},
            "r2": {"chrom": "chr1", "start": 190, "end": 210, "mapq": 35},
            "r3": {"chrom": "chr2", "start": 490, "end": 500, "mapq": 30},
            "lowq": {"chrom": "chr1", "start": 300, "end": 320, "mapq": 19},
            "mito": {"chrom": "chrM", "start": 20, "end": 30, "mapq": 40},
        }
        captured = {}

        def fake_run(command, **kwargs):
            captured["command"] = command
            captured["reads"] = Path(command[2]).read_text()
            captured["chromosomes"] = Path(command[3]).read_text()
            captured["genes"] = Path(command[4]).read_text()
            Path(command[5]).write_text("<svg/>")
            return ont_qc.subprocess.CompletedProcess(command, 0, "", "")

        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory) / "karyotype.svg"
            genes = Path(directory) / "genes.bed"
            genes.write_text("chr1\t0\t100\tgene1\n")
            with mock.patch.object(ont_qc.subprocess, "run", side_effect=fake_run):
                result = ont_qc.write_karyoploter(
                    records, sizes, ["chr1", "chr2"], 20, 10, "Rscript",
                    "karyotype_plot.R", genes, output, "Total library", "test")
            self.assertTrue(output.exists())
        self.assertEqual(result, {"mapped_reads": 3, "window_size_bp": 100})
        self.assertEqual(captured["reads"].splitlines(), [
            "chromosome\tposition", "chr1\t6", "chr1\t201", "chr2\t496"])
        self.assertNotIn("chrM", captured["chromosomes"])
        self.assertEqual(captured["genes"], "chr1\t0\t100\tgene1\n")
        self.assertEqual(captured["command"][7], "100")


class ReportLayoutTests(unittest.TestCase):
    def test_traffic_overview_is_technical_only(self):
        rendered = ont_qc.traffic_html({
            "scope_name": "Total library",
            "sequencing_yield_call": "pass",
            "adapter_dimer_call": "pass",
            "biological_qc_call": "low_quality",
        })
        self.assertIn("Sequencing", rendered)
        self.assertIn("Library quality", rendered)
        self.assertNotIn("Biological", rendered)

    def test_mycoplasma_unavailable_result_is_red_and_labeled_inconclusive(self):
        value = "INCONCLUSIVE_LOW_SCREEN_DEPTH"
        self.assertEqual(ont_qc.metric_light_class("mycoplasma_status", value, {}), "red")
        rendered = ont_qc.scoped_value_html("mycoplasma_status", value, {})
        self.assertIn("light red", rendered)
        self.assertIn("low screen depth", rendered)
        self.assertNotIn("not configured", rendered)

    def test_all_sample_overview_excludes_assay_specific_metrics(self):
        rendered = ont_qc.summary_overview_svg([{
            "sample": "S1", "assay": "microc", "raw_reads": 1000,
            "mean_read_q": 31.0, "complete_read_rate": 0.8,
            "raw_n50_bp": 1200, "insert_median_bp": 180,
            "primary_mapping_rate": 0.9, "adapter_dimer_rate": 0.01,
            "collision_pair_z": 99, "valid_contact_rate": 0.8,
            "biological_qc_call": "high_quality",
        }], ont_qc.DEFAULT_REFERENCE_RANGES)
        self.assertIn("ONT sequencing and NGS library QC", rendered)
        self.assertNotIn("Collision Z", rendered)
        self.assertNotIn("Biological QC", rendered)


if __name__ == "__main__":
    unittest.main()
