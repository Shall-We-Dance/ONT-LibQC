"""Pure metric formatting and report presentation resources."""
import html
from pathlib import Path
from qc_definitions import profile as metric_profile
from qc_reference_guide import context as reference_context
from report_library_structure import family_label

def call_color(call):
    call = call or ""
    if "pass" in call or "high_quality" in call:
        return "#16835D"
    if "low_quality" in call or "borderline" in call or "inconclusive" in call or "warning" in call:
        return "#C77A18"
    return "#C23B3B"


METRIC_META = {
    "collision_contributing_molecules": ("Contributing ONT molecules", "Distinct original ONT molecules contributing eligible chromatin fragments; units from one molecule are dependent."),
    "biological_qc_call": ("Cis-ratio evidence", "Confidence-aware endpoint cis screen only; cannot establish complete Micro-C library success."),
    "outer_complete_unit_rate": ("Both outer handles detected", "Reconstructed units with both outer P5 and P7 handles / all reconstructed units."),
    "outer_complete_units": ("Units with both outer handles", "Both outer P5/P7 handles detected, beyond the inner anchors required for reconstruction."),
    "index_assignment_rate": ("All-unit index assignment", "Assigned reconstructed units / all reconstructed units, including units missing index information."),
    "index_ambiguous_fraction": ("Ambiguous index fraction", "Units tied between index pairs / all reconstructed units."),
    "index_unassigned_fraction": ("Unassigned index fraction", "Units without an assignment / all reconstructed units."),
    "unreconstructed_probable_dimer_reads": ("Unresolved dimer candidates", "Unreconstructed raw reads matching a zero-insert template; excluded from confirmed burden because true short inserts can match."),
    "short_unmapped_units": ("Too short for two legs", "Non-dimer Micro-C inserts below 50 bp, excluded from mapping denominators."),
    "raw_reads": ("Raw reads", "ONT molecules supplied to the pipeline before adapter interpretation."),
    "sequencing_yield_call": ("ONT sequencing QC", "Technical sequencing-yield call based on the configured minimum raw-read and reconstructed-unit requirements."),
    "mean_read_q": ("Mean read Q", "Mean across raw molecules of -10 log10(mean base error probability), using all bases; differs from arithmetic mean base Q and basecaller-specific trimming."),
    "raw_median_bp": ("Raw median length", "Median full ONT molecule length before splitting."),
    "raw_n50_bp": ("Raw read N50", "Length above which half of all sequenced bases occur."),
    "insert_median_bp": ("Insert median length", "Median length of adapter-bounded NGS library inserts reconstructed from ONT molecules."),
    "complete_read_rate": ("Inner-anchor reconstruction rate", "Total: raw molecules containing both inner adapters; outer P5/P7 handles are not required. Sublibrary: contribution to all sample raw molecules, not a reconstruction success rate."),
    "library_units": ("Library units", "Number of adapter-bounded library units recovered from ONT molecules."),
    "units_per_complete_read": ("Units per complete read", "Number of inner-adapter-bounded units per reconstructed ONT molecule."),
    "mappable_library_units": ("Mappable units", "Non-dimer query units actually emitted for mapping; Micro-C also requires an insert of at least 50 bp (two nonoverlapping legs of at least 25 bp)."),
    "adapter_dimer_max_bp": ("Adapter-dimer cutoff", "Maximum adapter-bounded insert length counted as adapter dimer."),
    "adapter_dimer_units": ("Adapter-dimer units", "Units called by a canonical central gap or pair-specific zero-insert full-molecule template."),
    "canonical_adapter_dimer_units": ("Canonical adapter dimers", "Units with a central insert gap at or below the configured canonical cutoff."),
    "probable_adapter_dimer_units": ("Structure-matched adapter dimers", "Additional units matching their i5–i7-specific zero-insert full-molecule template."),
    "raw_structure_dimer_rate": ("Raw structure-dimer rate", "Fraction of raw ONT molecules supported as zero-insert dimers by reconstructed library structure."),
    "adapter_dimer_rate": ("Adapter-dimer burden", "Total: raw molecules with a reconstructed dimer-like unit / all raw molecules. Sublibrary: assigned raw molecules with a dimer-like unit / assigned raw molecules. Unreconstructed template candidates are excluded."),
    "adapter_dimer_unit_rate": ("Reconstructed-unit dimer rate", "Fraction of reconstructed library units called as zero-insert dimers; raw-molecule burden is used for the traffic-light call."),
    "adapter_dimer_call": ("Adapter-dimer QC", "Screening call based on adapter-dimer burden and available unit count."),
    "concatemer_rate": ("Concatemer rate", "Total: fraction of all raw molecules with multiple reconstructed units. Sublibrary: fraction of its assigned raw molecules with multiple units, including units from other scopes; shared molecules may contribute to multiple scopes."),
    "scope_unit_fraction": ("Sublibrary share of all units", "Assigned units in this sublibrary / all reconstructed sample units, including unassigned and ambiguous units."),
    "adapter_family": ("Adapter family", "Best-supported NGS-compatible P5/P7 sequence family."),
    "index_architecture": ("Index architecture", "Detected single-index, dual-index, mixed or unresolved architecture."),
    "i5": ("i5 consensus", "Consensus i5 index sequence normalized to P5-to-P7 orientation."),
    "i5_purity": ("i5 purity", "Fraction supporting the dominant i5 consensus among informative units."),
    "i7": ("i7 consensus", "Consensus i7 index sequence normalized to P5-to-P7 orientation."),
    "i7_purity": ("i7 purity", "Fraction supporting the dominant i7 consensus among informative units."),
    "aligner": ("Aligner", "Reference aligner selected for this assay."),
    "mycoplasma_status": ("Mycoplasma signal screen", "Opportunistic competitive-alignment signal screen; NOT_DETECTED_AT_THIS_DEPTH does not exclude contamination."),
    "mycoplasma_screened_reads": ("Mycoplasma-screened inserts", "Adapter-trimmed full inserts competitively aligned against host, eight common cell-culture Mollicutes and bacterial decoys."),
    "mycoplasma_qualifying_reads": ("Mycoplasma insert count", "Unique inserts passing target-contig, MAPQ, aligned-length, coverage and identity filters."),
    "mycoplasma_qualifying_loci": ("Mycoplasma locus count", "Distinct target genome 100-bp windows supported by qualifying inserts."),
    "mycoplasma_qualifying_fraction": ("Mycoplasma signal fraction", "Qualifying target reads divided by all screened adapter-trimmed inserts; a zero shallow result is not proof of absence."),
    "spike_in_status": ("Spike-in candidate screen", "Automatic competitive reference screen. Foreign-sequence evidence does not prove intentional spike-in; experimental metadata are required."),
    "spike_in_label": ("Candidate species", "Species with qualifying foreign-sequence assignments; same-species spike-in cannot be separated from host without specific barcodes."),
    "ligation_screened_units": ("Ligation-screened units", "All adapter-bounded non-dimer units, including unresolved cases. This is not the number of independent raw ONT molecules."),
    "ligation_resolved_units": ("Units with resolved block structure", "Query coverage >=80%, minimum aligned block 40 bp, MAPQ >=20 and <255, identity >=85%, and no unexplained query gap >30 bp."),
    "ligation_resolved_fraction": ("Resolved structure fraction", "Resolved units / all ligation-screened units."),
    "ligation_single_block_units": ("One continuous genomic block", "No resolved ligation junction. Adjacent or collinear ligated fragments can be indistinguishable from one block; this is not proof of unligated DNA."),
    "ligation_two_units": ("2-fragment units", "Units with two resolved genomic blocks; evidence-supported minimum fragment count."),
    "ligation_three_units": ("3-fragment units", "Units with three resolved genomic blocks; evidence-supported minimum fragment count."),
    "ligation_four_units": ("4-fragment units", "Units with four resolved genomic blocks; evidence-supported minimum fragment count."),
    "ligation_five_plus_units": ("≥5-fragment units", "Units with at least five resolved genomic blocks; short or repetitive fragments can remain undetected."),
    "ligation_unresolved_units": ("Unresolved ligation structure", "Insufficient unique mapping or query coverage, unexplained gaps, or large indels; kept in the denominator."),
    "ligation_multifragment_fraction": ("Multi-fragment evidence / all units", "Units with >=2 resolved genomic blocks / all screened library units, including unresolved cases in the denominator."),
    "ligation_multifragment_fraction_resolved": ("Multi-fragment evidence / resolved units", "Units with >=2 resolved blocks / units with resolved structure. Conditional on successful structure resolution."),
    "spike_in_screened_reads": ("Spike-in-screened inserts", "Adapter-trimmed inserts included in the host plus spike-in competitive alignment."),
    "spike_in_reads": ("Spike-in-assigned reads", "Unique inserts assigned to target spike-in contigs after MAPQ, alignment-length, coverage and identity filters."),
    "spike_in_fraction": ("Spike-in fraction", "Spike-in-assigned inserts divided by all inserts screened in the competitive reference."),
    "primary_mapping_rate": ("Primary-alignment rate", "Any mapped primary alignment / emitted query legs; primary alignment is distinct from configured primary chromosomes. MAPQ filtering is applied separately."),
    "usable_fragments": ("Usable fragments", "MAPQ-passing chromatin fragments on primary chromosomes."),
    "usable_fragment_rate": ("Usable fragment rate", "Usable chromatin fragments divided by mappable library units."),
    "mitochondrial_fraction": ("Mitochondrial fraction", "Mitochondrial / (configured primary nuclear + mitochondrial) MAPQ-pass alignments; excludes other contigs."),
    "collision_bin_size": ("Collision bin size", "Genomic bin size used for the primary collision statistic."),
    "collision_oe": ("Collision O/E (depth normalized)", "Observed same-bin pairs / expected pairs under uniform sampling of configured nuclear bases. Expected value is depth independent; 1 is random background. Low-depth estimates remain noisy; reference, bin size and mappability still matter. Not target-specific enrichment."),
    "collision_fragments": ("Collision observations", "MAPQ-pass nuclear chromatin fragments used for concentration."),
    "collision_observed_pairs": ("Observed same-bin pairs", "Sum over bins of n_j choose 2. Rare/zero collision counts at shallow depth may be inconclusive."),
    "collision_null_expected_pairs": ("Expected same-bin pairs", "N choose 2 times the uniform-base pair probability, corrected for chromosome terminal-bin widths."),
    "both_legs_mapped": ("Both legs mapped", "Micro-C units with a primary alignment for both inferred ligation legs."),
    "both_legs_mapping_rate": ("Both-leg mapping rate", "Fraction of mappable Micro-C units with both legs mapped."),
    "mapq_pass_pairs": ("MAPQ-pass pairs", "Both-leg Micro-C pairs passing MAPQ (excluding 255), query-coverage and primary-chromosome filters."),
    "mapq_pass_pair_rate": ("MAPQ-pass pair rate", "MAPQ-pass pairs divided by pairs with both legs mapped."),
    "same_fragment_pairs": ("Short inward pairs excluded", "Heuristic based on strand-aware aligned 5-prime separation; can exclude genuine local contacts and does not prove a shared fragment."),
    "same_fragment_filter_pass_pairs": ("Post-artifact-filter contacts", "MAPQ-pass pairs retained after the short inward-facing same-fragment filter."),
    "same_fragment_filter_pass_rate": ("Artifact-filter retention", "Post-artifact-filter contacts divided by MAPQ-pass pairs."),
    "unique_valid_contacts": ("Unique endpoint contacts", "Screened endpoint pairs after exact 5-prime coordinate deduplication within physical index scope; not junction-validated ligation events."),
    "valid_contact_rate": ("Valid-contact rate", "Unique valid contacts divided by mappable reconstructed library units."),
    "contact_deduplication_retention": ("Contact deduplication retention", "Unique valid contacts divided by post-artifact-filter contacts."),
    "duplicate_contact_fraction": ("Duplicate contact fraction", "Fraction of filtered Micro-C contacts removed as coordinate duplicates."),
    "cis_contacts": ("Cis contacts", "Unique valid Micro-C contacts whose legs map to the same chromosome."),
    "trans_contacts": ("Trans contacts", "Unique valid Micro-C contacts whose legs map to different chromosomes."),
    "cis_ratio": ("Cis ratio", "Fraction of unique valid Micro-C contacts mapping within the same chromosome."),
    "cis_ratio_ci95_low": ("Cis ratio CI low", "Lower bound of the Wilson 95% confidence interval for cis ratio."),
    "cis_ratio_ci95_high": ("Cis ratio CI high", "Upper bound of the Wilson 95% confidence interval for cis ratio."),
    "long_range_cis_ci95_low": ("Long-range cis CI low", "Lower Wilson 95% bound for long-range cis contacts divided by all unique contacts."),
    "long_range_cis_ci95_high": ("Long-range cis CI high", "Upper Wilson 95% bound for long-range cis contacts divided by all unique contacts."),
    "cis_lt1kb": ("Cis contacts <1 kb", "Unique cis endpoint contacts separated by less than 1 kb."),
    "cis_1_10kb": ("Cis contacts 1–10 kb", "Unique cis endpoint contacts separated by at least 1 kb and less than 10 kb."),
    "cis_10_100kb": ("Cis contacts 10–100 kb", "Unique cis endpoint contacts separated by at least 10 kb and less than 100 kb."),
    "cis_ge100kb": ("Cis contacts ≥100 kb", "Unique cis endpoint contacts separated by at least 100 kb."),
    "long_range_cis_min_bp": ("Long-range cis threshold", "Minimum genomic separation used to define a long-range cis contact."),
    "long_range_cis_contacts": ("Long-range cis contacts", "Unique valid cis contacts at or above the configured genomic-separation threshold."),
    "long_range_cis_fraction": ("Long-range cis fraction", "Long-range cis contacts divided by all unique valid contacts."),
    "long_range_within_cis_fraction": ("Long-range fraction within cis", "Long-range cis contacts divided only by unique valid cis contacts."),
}

PERCENT_METRICS = {
    "complete_read_rate", "adapter_dimer_rate", "raw_structure_dimer_rate", "concatemer_rate", "i5_purity", "i7_purity",
    "primary_mapping_rate",
    "usable_fragment_rate", "mitochondrial_fraction",
"both_legs_mapping_rate", "duplicate_contact_fraction",
    "cis_ratio", "cis_ratio_ci95_low", "cis_ratio_ci95_high", "long_range_cis_ci95_low", "long_range_cis_ci95_high", "mycoplasma_qualifying_fraction",
    "spike_in_fraction", "ligation_resolved_fraction", "ligation_multifragment_fraction", "ligation_multifragment_fraction_resolved", "mapq_pass_pair_rate", "same_fragment_filter_pass_rate",
    "valid_contact_rate", "contact_deduplication_retention",
    "long_range_cis_fraction", "long_range_within_cis_fraction",
    "adapter_dimer_unit_rate", "outer_complete_unit_rate", "index_assignment_rate", "index_ambiguous_fraction", "index_unassigned_fraction", "scope_unit_fraction",
}

DEFAULT_REFERENCE_RANGES = {
    "raw_reads": {"min": 500, "warning_min": 250},
    "mean_read_q": None,
    "complete_read_rate": {"min": .50, "max": 1.0, "warning_min": .30, "type": "percent"},
    "adapter_dimer_rate": {"min": 0.0, "max": .05, "warning_max": .10, "type": "percent"},
    "primary_mapping_rate": {"min": .50, "max": 1.0, "warning_min": .30, "type": "percent"},
    "usable_fragment_rate": {"min": .30, "max": 1.0, "warning_min": .15, "type": "percent"},
    "mitochondrial_fraction": {"min": 0.0, "max": .20, "warning_max": .40, "type": "percent"},
    "both_legs_mapping_rate": {"min": .50, "max": 1.0, "warning_min": .30, "type": "percent"},
    "cis_ratio": None,
}


def fmt(value):
    if value is None:
        return "NA"
    if isinstance(value, float):
        return f"{value:.4g}"
    return str(value)


def is_percent_metric(key, ranges):
    return key in PERCENT_METRICS or (ranges.get(key) or {}).get("type") == "percent"


def display_value(key, value, ranges):
    if value == "inconclusive_low_yield":
        return "low"
    if key == "adapter_family" and value is not None:
        return family_label(value)
    if key == "mycoplasma_status" and value == "NOT_DETECTED_AT_THIS_DEPTH":
        return "NOT_DETECTED"
    if value is None:
        return "NA"
    if is_percent_metric(key, ranges) and isinstance(value, (int, float)):
        if key in {"mycoplasma_qualifying_fraction", "spike_in_fraction"} and 0 < value < .001:
            return f"{value * 100:.3g}%"
        return f"{value:.1%}"
    if key.endswith("_bp") or key == "collision_bin_size":
        return f"{value:,.0f} bp" if isinstance(value, (int, float)) else str(value)
    if isinstance(value, int):
        return f"{value:,}"
    if isinstance(value, float):
        return f"{value:.4g}"
    return str(value)


def reference_interval_text(key, ranges, metric=None):
    spec = ranges.get(key) or {}
    minimum, maximum = spec.get("min"), spec.get("max")
    if minimum is not None and maximum is not None:
        return f"{display_value(key, minimum, ranges)}–{display_value(key, maximum, ranges)}"
    if minimum is not None:
        return f"≥{display_value(key, minimum, ranges)}"
    if maximum is not None:
        return f"≤{display_value(key, maximum, ranges)}"
    return reference_context(key, metric)[0]


def merged_reference_ranges(cfg):
    ranges = {key: dict(value) for key, value in DEFAULT_REFERENCE_RANGES.items() if value is not None}
    for key, value in (cfg.get("qc_reference_ranges") or {}).items():
        if value is None:
            ranges.pop(key, None)
        else:
            ranges[key] = {**ranges.get(key, {}), **value}
    ranges.pop("cis_ratio", None)
    return ranges


def range_status(key, value, ranges):
    spec = ranges.get(key)
    if not spec or not isinstance(value, (int, float)):
        return None
    if spec.get("min") is not None and value < spec["min"]:
        return "low"
    if spec.get("max") is not None and value > spec["max"]:
        return "high"
    return "normal"


def metric_light_class(key, value, ranges):
    if key in {"sequencing_yield_call", "adapter_dimer_call", "biological_qc_call"}:
        status = str(value or "").lower()
        if not status:
            return "gray"
        if "pass" in status or "high_quality" in status:
            return "green"
        if any(word in status for word in ("warning", "borderline", "inconclusive")):
            return "amber"
        return "red"
    if key == "spike_in_status":
        status = str(value or "").upper()
        if status == "NOT_CONFIGURED" or not status:
            return "gray"
        if status.startswith("INCONCLUSIVE"):
            return "amber"
        if status in {"CANDIDATE_SIGNAL", "TRACE"}:
            return "amber"
        if status in {"DETECTED", "NOT_DETECTED"}:
            return "gray"
        return "red"
    if key == "mycoplasma_status":
        status = str(value or "").upper()
        return "green" if status.startswith("NOT_DETECTED") else "red"
    if key in {"mycoplasma_qualifying_reads", "mycoplasma_qualifying_loci", "mycoplasma_qualifying_fraction"}:
        return "green" if isinstance(value, (int, float)) and value == 0 else "red"
    status = range_status(key, value, ranges)
    if status is None:
        return "gray"
    if status == "normal":
        return "green"
    spec = ranges[key]
    if status == "low" and spec.get("warning_min") is not None and value >= spec["warning_min"]:
        return "amber"
    if status == "high" and spec.get("warning_max") is not None and value <= spec["warning_max"]:
        return "amber"
    return "red"


def assay_reference_ranges(cfg, assay):
    ranges = merged_reference_ranges(cfg)
    if assay in {"rnaseq", "clip"}:
        for key in ("primary_mapping_rate", "usable_fragment_rate", "mitochondrial_fraction"):
            if key not in cfg.get("qc_reference_ranges", {}):
                ranges[key] = None
    return ranges


def summary_overview_svg(metrics, ranges, assay_ranges=None):
    """Technical quick overview: ONT sequencing and NGS library shape only."""
    columns = [
        ("Sample", 30, None),
        ("Assay", 130, None),
        ("Raw reads", 220, "raw_reads"),
        ("Mean Q", 365, "mean_read_q"),
        ("Reconstructed", 495, "complete_read_rate"),
        ("Raw N50", 650, "raw_n50_bp"),
        ("Insert median", 820, "insert_median_bp"),
        ("Mapping", 1000, "primary_mapping_rate"),
        ("Adapter dimer", 1165, "adapter_dimer_rate"),
    ]
    sample_extra = max(0, max((len(m["sample"]) for m in metrics), default=0) * 8 - 85)
    columns = [(label, x + (sample_extra if i else 0), key) for i, (label, x, key) in enumerate(columns)]
    width, row_h = 1370 + sample_extra, 52
    height = 135 + row_h * len(metrics)
    light_colors = {"green": "#16835D", "amber": "#E2A12A", "red": "#C23B3B", "gray": "#9BA8B2"}
    svg = [
        f"<svg xmlns='http://www.w3.org/2000/svg' viewBox='0 0 {width} {height}' role='img' aria-label='All samples technical QC summary' style='font-family:Arial,Helvetica,sans-serif'>",
        "<rect width='100%' height='100%' fill='#F5F8FA'/>",
        "<text x='30' y='40' font-size='28' font-weight='700' fill='#132A3E'>ONT sequencing and NGS library QC — all samples</text>",
        "<text x='30' y='65' font-size='13' fill='#607080'>Quick overview contains technical sequencing and library type/shape metrics only.</text>",
    ]
    for label, x, _ in columns:
        svg.append(f"<text x='{x}' y='102' font-size='12' font-weight='700' fill='#536475'>{html.escape(label)}</text>")
    for index, metric in enumerate(metrics):
        current_ranges = (assay_ranges or {}).get(metric["assay"], ranges)
        y = 132 + index * row_h
        if index % 2 == 0:
            svg.append(f"<rect x='18' y='{y-23}' width='{width-36}' height='{row_h}' rx='8' fill='white'/>")
        svg.extend([
            f"<text x='30' y='{y+5}' font-size='14' font-weight='700' fill='#18324A'>{html.escape(metric['sample'])}</text>",
            f"<text x='{130+sample_extra}' y='{y+5}' font-size='12' fill='#607080'>{html.escape({'chromatin':'Chromatin','microc':'Micro-C','rnaseq':'RNA-seq · beta','clip':'CLIP-seq · beta'}.get(metric['assay'],metric['assay']))}</text>",
        ])
        for _, x, key in columns[2:]:
            value = metric.get(key)
            color = light_colors["amber" if key == "adapter_dimer_rate" and str(metric.get("adapter_dimer_call", "")).startswith("inconclusive") else metric_light_class(key, value, current_ranges)]
            svg.append(f"<circle cx='{x+7}' cy='{y}' r='7' fill='{color}'/>")
            svg.append(f"<text x='{x+20}' y='{y+5}' font-size='11.5' fill='#34495E'>{html.escape(display_value(key, value, ranges))}</text>")
    svg.append(f"<text x='30' y='{height-9}' font-size='10' fill='#71808D'>Assay-specific evidence is available in each sample report.</text></svg>")
    return "".join(svg)


REPORT_CSS = Path(__file__).with_name("report_base.css").read_text()

REPORT_CSS += Path(__file__).with_name("report_modern.css").read_text()


REPORT_TOOLTIP_SCRIPT = "<script>" + Path(__file__).with_name("report_tooltips.js").read_text() + "</script>"



def svg_contact_funnel(metric, title):
    stages = [
        ("Mappable units", metric.get("mappable_library_units")),
        ("Both legs mapped", metric.get("both_legs_mapped")),
        ("MAPQ pass", metric.get("mapq_pass_pairs")),
        ("Artifact filter pass", metric.get("same_fragment_filter_pass_pairs")),
        ("Unique endpoint contacts", metric.get("unique_valid_contacts")),
    ]
    width, height, left, top, bar_width, row_height = 900, 400, 235, 62, 600, 58
    values = [int(value or 0) for _, value in stages]
    peak = max(values[:1] + [1])
    elements = [
        f"<svg xmlns='http://www.w3.org/2000/svg' viewBox='0 0 {width} {height}' style='font-family:Arial,Helvetica,sans-serif'>",
        "<rect width='100%' height='100%' fill='white'/>",
        f"<text x='{width/2}' y='27' text-anchor='middle' font-size='17'>{html.escape(title)}</text>",
    ]
    previous = None
    colors = ["#356A91", "#3E7C8A", "#4A8D79", "#5A9C69", "#16835D"]
    for index, ((label, _), value, color) in enumerate(zip(stages, values, colors)):
        y = top + index * row_height
        length = bar_width * value / peak
        retention = None if previous in (None, 0) else value / previous
        suffix = "" if retention is None else f" · {retention:.1%} of previous"
        elements.extend([
            f"<text x='{left-14}' y='{y+18}' text-anchor='end' font-size='13' fill='#34495E'>{html.escape(label)}</text>",
            f"<rect x='{left}' y='{y}' width='{bar_width}' height='27' rx='6' fill='#E7ECEF'/>",
            f"<rect x='{left}' y='{y}' width='{length:.1f}' height='27' rx='6' fill='{color}'/>",
            f"<text x='{left+8}' y='{y+43}' font-size='12' fill='#34495E'>{value:,}{html.escape(suffix)}</text>",
        ])
        previous = value
    threshold = metric.get("long_range_cis_min_bp")
    long_fraction = metric.get("long_range_cis_fraction")
    note = "Long-range cis: NA"
    if threshold is not None and long_fraction is not None:
        note = f"Long-range cis (≥{threshold:,} bp / all unique valid contacts): {long_fraction:.1%}"
    elements.append(f"<text x='{left}' y='{height-18}' font-size='12' fill='#536475'>{html.escape(note)}</text></svg>")
    return "".join(elements)


def svg_contact_evidence(metric, title):
    width, height, left, span = 900, 340, 205, 630
    parts = [f"<svg xmlns='http://www.w3.org/2000/svg' viewBox='0 0 {width} {height}' style='font-family:Arial,Helvetica,sans-serif'>",
             "<rect width='100%' height='100%' fill='white'/>",
             f"<text x='450' y='26' text-anchor='middle' font-size='17'>{html.escape(title)}</text>",
             f"<text x='450' y='47' text-anchor='middle' font-size='11'>n={metric.get('unique_valid_contacts',0):,} unique endpoint contacts; Wilson 95% intervals</text>"]
    for i, (label, key, low, high) in enumerate([
        ("Cis / all contacts", "cis_ratio", "cis_ratio_ci95_low", "cis_ratio_ci95_high"),
        ("Long cis / all contacts", "long_range_cis_fraction", "long_range_cis_ci95_low", "long_range_cis_ci95_high")]):
        y = 85 + i*75
        parts.append(f"<text x='{left-12}' y='{y+4}' text-anchor='end' font-size='12'>{label}</text>")
        parts.append(f"<line x1='{left}' x2='{left+span}' y1='{y}' y2='{y}' stroke='#dce3e8' stroke-width='9'/>")
        value = metric.get(key)
        if value is not None and metric.get(low) is not None:
            x, xlo, xhi = [left+span*v for v in [value, metric[low], metric[high]]]
            parts.append(f"<line x1='{xlo:.2f}' x2='{xhi:.2f}' y1='{y}' y2='{y}' stroke='#245d7e' stroke-width='3'/><circle cx='{x:.2f}' cy='{y}' r='5' fill='#245d7e'/>")
            parts.append(f"<text x='{left}' y='{y+24}' font-size='12'>{value:.1%} [{metric[low]:.1%}, {metric[high]:.1%}]</text>")
        else:
            parts.append(f"<text x='{left}' y='{y+24}' font-size='12'>NA: no eligible contacts</text>")
    for i in range(6):
        parts.append(f"<text x='{left+i*span/5}' y='215' text-anchor='middle' font-size='11'>{i*20}%</text>")
    parts.append("<text x='450' y='289' text-anchor='middle' font-size='11'>Cis ratio alone does not validate ligation; short inward exclusion is heuristic.</text>")
    parts.append("<text x='450' y='311' text-anchor='middle' font-size='11'>Intervals assume independent contacts; ONT recovery biases are not represented.</text></svg>")
    return "".join(parts)


def call_class(call):
    call=call or ""
    if not call or "not_classified" in call or "NOT_CONFIGURED" in call: return "gray"
    if "pass" in call or "high_quality" in call: return "green"
    if any(x in call for x in ("warning","low_quality","inconclusive","insufficient","borderline")): return "amber"
    return "red"


def traffic_html(scope, ranges=None, assay=None):
    ranges=DEFAULT_REFERENCE_RANGES if ranges is None else ranges
    assay=assay or scope.get("assay")
    label=scope.get("scope_name",scope.get("scope_id","Total"))
    seq_call=scope.get("sequencing_yield_call")
    shape_keys=["adapter_dimer_rate","primary_mapping_rate"]
    if assay == "chromatin":
        shape_keys.extend(["usable_fragment_rate", "mitochondrial_fraction"])
    if scope.get("scope_id", "total") == "total": shape_keys.append("complete_read_rate")
    evaluated=[metric_light_class(k,scope.get(k),ranges) for k in shape_keys if scope.get(k) is not None and ranges.get(k)]
    if str(seq_call or "").startswith(("inconclusive", "insufficient")) or scope.get("raw_reads",0) < scope.get("metric_settings",metric_profile())["minimum_raw_reads"] or scope.get("library_units",0) < scope.get("metric_settings",metric_profile())["minimum_complete_units"]:
        shape_call="insufficient_evidence"
    elif "red" in evaluated:shape_call="fail"
    elif "amber" in evaluated:shape_call="warning"
    elif evaluated:shape_call="pass"
    else:shape_call=scope.get("adapter_dimer_call") or "insufficient_evidence"
    quality_tip="Combined adapter-dimer burden, primary mapping and total-library adapter recovery checks."
    if assay == "chromatin":
        quality_tip+=" Includes usable fragment rate and mitochondrial fraction."
    flagged=[f"{METRIC_META[k][0]}: {display_value(k, scope[k], ranges)} ({'outside range' if metric_light_class(k,scope[k],ranges)=='red' else 'borderline'})."
             for k in shape_keys if scope.get(k) is not None and ranges.get(k) and metric_light_class(k,scope[k],ranges) in {"red", "amber"}]
    quality_tip+=" Not a biological enrichment or ligation success certificate."
    if flagged: quality_tip+=" "+" ".join(flagged)
    items=[("Sequencing",seq_call,"Yield screen; inspect read quality and molecule recovery below."),("Library quality",shape_call,quality_tip)]
    return f"<section class='scope-card'><h3>{html.escape(label)}</h3>"+''.join(f"<div class='traffic' title='{html.escape(tip,quote=True)}'><span class='light {call_class(v)}'></span><b>{html.escape(k)}</b><em>{html.escape(display_value(k, v, ranges).replace('_',' ').capitalize())}</em></div>" for k,v,tip in items)+"</section>"


def metric_reading_html(rendered, light, state):
    unavailable = rendered in {"Not attributable", "Sample-level"}
    classes = " metric-unavailable" if unavailable else ""
    marker = "" if light == "gray" else f"<span class='light {light}' aria-hidden='true'></span>"
    value = "—" if unavailable else rendered.replace("NOT_CONFIGURED", "Not configured")
    if "_" in value and value.isupper():
        value = value.replace("_", " ").capitalize()
    return (f"<div class='metric-cell'><div class='metric-reading metric-tip metric-status{classes}' tabindex='0' data-tip='{html.escape(state, quote=True)}'"
            f" aria-label='{html.escape(rendered + ': ' + state, quote=True)}'>"
            f"{marker}<span>{html.escape(value)}</span></div></div>")


def scoped_value_html(key, value, ranges, descriptive=False):
    rendered = display_value(key, value, ranges)
    light = "gray" if descriptive else metric_light_class(key, value, ranges)
    if value == "NOT_CONFIGURED":
        state = "Not configured"
    elif descriptive:
        state = "Descriptive only"
    elif key == "mycoplasma_status":
        status = str(value or "").upper()
        if status.startswith("INCONCLUSIVE"):
            state = "Inconclusive: low screen depth"
        elif status.startswith("NOT_DETECTED"):
            state = "No signal detected at this depth"
        else:
            state = "Mycoplasma signal detected" if light == "red" else "Not configured"
    elif key in {"mycoplasma_qualifying_reads", "mycoplasma_qualifying_loci", "mycoplasma_qualifying_fraction"}:
        state = "Mycoplasma signal detected" if light == "red" else "No signal detected at this depth"
    elif key == "spike_in_status":
        state = "Candidate foreign-sequence evidence; intentional spike-in requires experimental metadata" if str(value) in {"CANDIDATE_SIGNAL", "TRACE"} else "No candidate signal at this depth" if value == "NOT_DETECTED" else "Reference screening result; experimental spike-in intent cannot be inferred"
    else:
        state = {"green":"Within range", "amber":"Borderline", "red":"Outside range", "gray":"Descriptive only"}[light]
    return metric_reading_html(rendered, light, state)



def metric_missing_reason(scope, key):
    reason = scope.get("metric_unavailable_reasons", {}).get(key)
    if reason:
        return reason
    if key.startswith("external_frip") and scope.get("external_frip_status") == "NOT_CONFIGURED":
        return "Not configured"
    if key.startswith("spike_in_") and scope.get("spike_in_status") == "NOT_CONFIGURED":
        return "Not configured"
    if key.startswith("mycoplasma_") and scope.get("mycoplasma_status") == "NOT_CONFIGURED":
        return "Not configured: Mycoplasma competitive screening is disabled."
    if key.startswith("mycoplasma_") and scope.get("mycoplasma_status") == "REQUIRES_SCREEN_RERUN":
        return "Rerun screen: per-insert supporting evidence is missing from the legacy result."
    return "Insufficient observations to estimate this metric."


def scope_metric_cell(scope, key, ranges):
    value = scope.get(key)
    if key in {"cis_ratio", "long_range_cis_fraction", "valid_contact_rate"} and value is not None:
        denominator=scope.get("unique_valid_contacts",0) if key!='valid_contact_rate' else scope.get('mappable_library_units',0)
        numerator=scope.get("cis_contacts",0) if key=='cis_ratio' else scope.get('long_range_cis_contacts',0) if key=='long_range_cis_fraction' else scope.get('unique_valid_contacts',0)
        text=f"{value:.1%}"
        if key!='valid_contact_rate':
            prefix='cis_ratio' if key=='cis_ratio' else 'long_range_cis'
            lo,hi=scope.get(prefix+'_ci95_low'),scope.get(prefix+'_ci95_high')
            if lo is not None and hi is not None:text+=f" [{lo:.1%}, {hi:.1%}]"
        text+=f" · {numerator:,}/{denominator:,}"
        evidence=scope.get('contact_evidence_status','')
        tip='Wilson 95% interval assumes independent contacts; count sufficiency does not establish biological success.'
        if evidence.startswith('insufficient'):tip='Insufficient evidence: too few valid contacts. '+tip
        return f"<td>{metric_reading_html(text,'amber' if evidence.startswith('insufficient') else 'gray',tip)}</td>"
    if key == "collision_oe" and value is not None:
        tip=f"{scope.get('collision_bin_size',10000)/1000:g} kb bins; {scope.get('collision_observed_pairs',0):,} observed pairs; {scope.get('collision_null_expected_pairs',0):.4g} expected under uniform-base sampling. Concentration is not target-specific enrichment."
        return f"<td>{metric_reading_html(display_value(key,value,ranges),'gray',tip)}</td>"
    if key.startswith("mycoplasma_") and value is not None and scope.get("mycoplasma_screened_reads") == 0:
        return f"<td>{metric_reading_html(display_value(key,value,ranges), 'red', 'No inserts screened; mycoplasma cannot be assessed.')}</td>"
    if key == "mycoplasma_status" and scope.get("mycoplasma_screened_reads", 0) and scope.get("mycoplasma_qualifying_fraction") == 0:
        value = "NOT_DETECTED_AT_THIS_DEPTH"
    if value is None:
        reason = metric_missing_reason(scope, key)
        label = "NOT_CONFIGURED" if reason.startswith("Not configured") else "Not attributable" if reason.startswith("Not attributable") else "Sample-level" if reason.startswith("Sample-level") else "Not estimable"
        missing_light='red' if key.startswith('mycoplasma_') else 'gray'
        return f"<td>{metric_reading_html(label, missing_light, reason)}</td>"
    descriptive = ((key == 'complete_read_rate' and scope.get('scope_id') != 'total')
                   or (key == 'cis_ratio' and str(scope.get('biological_qc_call', '')).startswith('inconclusive'))
                   or (key == 'adapter_dimer_rate' and str(scope.get('adapter_dimer_call', '')).startswith('inconclusive')))
    return f"<td>{scoped_value_html(key, value, ranges, descriptive)}</td>"

def scoped_table_html(scopes, keys, ranges):
    heads=''.join(f"<th>{html.escape(x.get('scope_name',x.get('scope_id','Total')))}</th>" for x in scopes)
    rows=[]
    for key in keys:
        label,desc=METRIC_META.get(key,(key,f"Pipeline QC metric: {key}."))
        interval=html.escape(reference_interval_text(key,ranges,scopes[0] if scopes else None))
        _,basis,detail=reference_context(key,scopes[0] if scopes else None)
        if ranges.get(key):
            basis="Local screening setting"
            detail="Numeric range from the effective pipeline configuration, not a universal assay standard."
            for boundary,boundary_label in [("warning_min","Lower warning boundary"),("warning_max","Upper warning boundary")]:
                if ranges[key].get(boundary) is not None:detail+=f" {boundary_label}: {display_value(key,ranges[key][boundary],ranges)}."
            if key=="complete_read_rate":detail+=" The threshold applies to total-library reconstruction. A sublibrary value is its contribution to all sample raw molecules, not a separate reconstruction success rate."
            if key in {"adapter_dimer_rate","primary_mapping_rate","usable_fragment_rate","both_legs_mapping_rate","mitochondrial_fraction","cis_ratio"}:
                detail+=" These ONT measurements and denominators require calibration against matched libraries and controls."
        interval=f"<span class='metric-tip reference-context' tabindex='0' data-tip='{html.escape(basis+': '+detail,quote=True)}'>{interval}</span>"
        cells=''.join(scope_metric_cell(x, key, ranges) for x in scopes)
        rows.append(f"<tr><td><span class='metric-tip' tabindex='0' data-tip='{html.escape(desc,quote=True)}'>{html.escape(label)} <b>ⓘ</b></span></td>{cells}<td class='reference-interval'>{interval}</td></tr>")
    return f"<div class='table-scroll' tabindex='0' aria-label='QC metrics; scroll horizontally for additional columns'><table class='qc-metrics'><tr><th>Metric</th>{heads}<th>Reference</th></tr>{''.join(rows)}</table></div>"


