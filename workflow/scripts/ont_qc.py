#!/usr/bin/env python3
"""ONT-LibQC CLI and backwards-compatible public helper imports."""
from qc_definitions import profile as metric_profile
import json
import argparse
import subprocess
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from library_sequences import (
    FAMILIES,
    NEX_R1,
    NEX_R2,
    P5,
    P7,
    TRU_R1,
    TRU_R2,
    edit_distance,
    fastq_records,
    open_text,
    rc
)
from qc_analysis import (
    add_indexes,
    approximate_matches,
    auto_mitochondrial_chromosomes,
    base_unit_id,
    bin_stats,
    build_sublibrary_gene_profile,
    cigar_ref_span,
    cluster_index_pairs,
    cmd_metrics,
    cmd_mycoplasma,
    cmd_prepare,
    cmd_spikein,
    cmd_track,
    contact_call,
    contact_metrics,
    contributing_molecules,
    external_frip,
    find_units,
    five_prime,
    full_length_template_identity,
    full_library_profile,
    index_architecture,
    index_close,
    index_consensus,
    index_pair_distance,
    interval_overlaps,
    json_dump,
    length_summary,
    load_gene_annotation,
    matches_in_window,
    mean_error_q,
    merge_intervals,
    microc_contacts,
    mycoplasma_call,
    n50,
    nearest_after,
    nearest_before,
    parse_expected_index_pairs,
    quantile,
    read_chrom_sizes,
    read_index_assignments,
    read_numbers,
    read_scope_profiles,
    resolve_chromosome_sets,
    rna_mapping_metrics,
    sam_primary_records,
    scope_qc,
    scope_structure_metrics,
    scoped_mycoplasma,
    single_index_p5_hit,
    svg_bars,
    svg_hist,
    svg_profile_placeholder,
    wilson,
    write_karyoploter,
    write_library_size_assets,
    write_sparse_cpm_bedgraph,
    zero_insert_template
)
from report_metrics import (
    DEFAULT_REFERENCE_RANGES,
    METRIC_META,
    PERCENT_METRICS,
    REPORT_CSS,
    REPORT_TOOLTIP_SCRIPT,
    assay_reference_ranges,
    call_class,
    call_color,
    display_value,
    fmt,
    is_percent_metric,
    merged_reference_ranges,
    metric_light_class,
    metric_missing_reason,
    metric_reading_html,
    range_status,
    reference_interval_text,
    scope_metric_cell,
    scoped_table_html,
    scoped_value_html,
    summary_overview_svg,
    svg_contact_evidence,
    svg_contact_funnel,
    traffic_html
)
from report_render import (
    cmd_report,
    ecoli_calibration_note_html,
    extension_report_html,
    identity_panel_html,
    library_structure_html,
    pair_label,
    profile_scope_html,
    rna_screening_html
)

def parser():
    from scrna_qc import add_commands as add_scrna_commands
    p = argparse.ArgumentParser()
    sub = p.add_subparsers(dest="command", required=True)
    add_scrna_commands(sub)
    x = sub.add_parser("prepare")
    x.add_argument("--sample", required=True); x.add_argument("--assay", choices=["chromatin", "microc", "rnaseq", "clip"], required=True)
    x.add_argument("--fastq", nargs="+", required=True); x.add_argument("--inserts", required=True); x.add_argument("--legs", required=True)
    x.add_argument("--leg1", required=True); x.add_argument("--leg2", required=True)
    x.add_argument("--read-table", required=True); x.add_argument("--index-assignments", required=True); x.add_argument("--summary", required=True)
    x.add_argument("--raw-lengths", required=True); x.add_argument("--insert-lengths", required=True)
    x.add_argument("--leg-length", type=int, default=100); x.add_argument("--max-unit", type=int, default=5000)
    x.add_argument("--adapter-dimer-max-bp", type=int, default=10)
    x.add_argument("--ultrashort-insert-max-bp", type=int, default=50)
    x.add_argument("--adapter-dimer-raw-min-bp", type=int, default=120); x.add_argument("--adapter-dimer-raw-max-bp", type=int, default=170)
    x.add_argument("--adapter-dimer-template-identity", type=float, default=.80)
    x.add_argument("--single-index-pair-fraction", type=float, default=.90)
    x.add_argument("--minimum-index-informative-units", type=int, default=20)
    x.add_argument("--minimum-index-cluster-units", type=int, default=20)
    x.add_argument("--minimum-index-cluster-fraction", type=float, default=.01)
    x.add_argument("--index-variant-max-edit", type=int, default=2,
                   help="Maximum total i5+i7 edit distance for ONT index-variant rescue")
    x.add_argument("--index-variant-dominance-ratio", type=float, default=5.0,
                   help="Minimum support ratio of a dominant pair over a rescued variant")
    x.add_argument("--expected-index-pairs", default="")
    x.add_argument("--library-structure", default="{}", help="JSON library layout with fixed-length UMI positions")
    x.set_defaults(func=cmd_prepare)
    x = sub.add_parser("metrics")
    x.add_argument("--sample", required=True); x.add_argument("--assay", choices=["chromatin", "microc", "rnaseq", "clip"], required=True)
    x.add_argument("--bam", required=True); x.add_argument("--samtools", default="samtools"); x.add_argument("--adapter-summary", required=True)
    x.add_argument("--read-table", required=True); x.add_argument("--index-assignments", required=True)
    x.add_argument("--assembly", default="unspecified"); x.add_argument("--mycoplasma", required=True); x.add_argument("--spike-in", required=True)
    x.add_argument("--chrom-sizes", required=True); x.add_argument("--primary", default="auto"); x.add_argument("--mito", default="auto")
    x.add_argument("--aligner", default="unknown")
    x.add_argument("--mapq", type=int, default=metric_profile()["mapq"])
    x.add_argument("--collision-bin-size", type=int, default=metric_profile()["collision_qc_bin_size"])
    x.add_argument("--microc-cis-pass", type=float, default=None, help=argparse.SUPPRESS); x.add_argument("--microc-cis-fail", type=float, default=None, help=argparse.SUPPRESS)
    x.add_argument("--microc-long-range-cis-min-bp", type=int, default=metric_profile()["microc_long_range_cis_min_bp"])
    x.add_argument("--microc-minimum-contacts", type=int, default=metric_profile()["microc_minimum_contacts"])
    x.add_argument("--microc-short-inward-bp", type=int, default=metric_profile()["microc_short_inward_bp"])
    x.add_argument("--microc-minimum-query-coverage", type=float, default=metric_profile()["microc_minimum_query_coverage"])
    x.add_argument("--dimer-pass-max", type=float, default=.05)
    x.add_argument("--dimer-warning-max", type=float, default=.10)
    x.add_argument("--external-peaks", default="")
    x.add_argument("--minimum-raw", type=int, default=metric_profile()["minimum_raw_reads"]); x.add_argument("--minimum-units", type=int, default=metric_profile()["minimum_complete_units"])
    x.add_argument("--ligation", default="", help="Full-insert split-alignment ligation JSON for Micro-C")
    x.add_argument("--metrics", required=True)
    x.add_argument("--fragments", required=True); x.add_argument("--pairs", required=True)
    x.set_defaults(func=cmd_metrics)
    x = sub.add_parser("spikein")
    x.add_argument("--fastq", required=True); x.add_argument("--output", required=True)
    x.add_argument("--enabled", action="store_true"); x.add_argument("--index", default="")
    x.add_argument("--target-regex", default=""); x.add_argument("--label", default="spike-in")
    x.add_argument("--bowtie2", default="bowtie2"); x.add_argument("--mapq", type=int, default=20)
    x.add_argument("--minimum-aligned-bp", type=int, default=80)
    x.add_argument("--minimum-query-coverage", type=float, default=.80)
    x.add_argument("--minimum-identity", type=float, default=.85)
    x.add_argument("--minimum-screened-reads", type=int, default=500)
    x.set_defaults(func=cmd_spikein)
    x = sub.add_parser("mycoplasma")
    x.add_argument("--fastq", required=True); x.add_argument("--output", required=True)
    x.add_argument("--enabled", action="store_true"); x.add_argument("--index", default=""); x.add_argument("--target-regex", default="mycoplasma_target|mycoplasmopsis|mycoplasmoides|metamycoplasma|mesomycoplasma|acholeplasma|spiroplasma|ureaplasma|mollicutes")
    x.add_argument("--bowtie2", default="bowtie2"); x.add_argument("--mapq", type=int, default=20)
    x.add_argument("--minimum-aligned-bp", type=int, default=80); x.add_argument("--minimum-query-coverage", type=float, default=.80); x.add_argument("--minimum-identity", type=float, default=.85)
    x.add_argument("--minimum-screened-reads", type=int, default=500)
    x.set_defaults(func=cmd_mycoplasma)
    x = sub.add_parser("track")
    x.add_argument("--bam", required=True); x.add_argument("--output", required=True)
    x.add_argument("--samtools", default="samtools"); x.add_argument("--converter", default="bedGraphToBigWig")
    x.add_argument("--chrom-sizes", required=True); x.add_argument("--primary", default="auto"); x.add_argument("--mapq", type=int, default=20)
    x.set_defaults(func=cmd_track)
    x = sub.add_parser("report")
    x.add_argument("--reuse-plot-assets", action="store_true", help="Reuse existing karyotype and sublibrary metaplot assets. Remove assets affected by input/config changes first; other plots and tables are refreshed.")
    x.add_argument("--config", required=True); x.add_argument("--out", required=True); x.set_defaults(func=cmd_report)
    return p


if __name__ == "__main__":
    arguments = parser().parse_args()
    arguments.func(arguments)
