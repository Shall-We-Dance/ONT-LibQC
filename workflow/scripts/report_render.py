"""Report assembly using retained analysis results and presentation helpers."""
from profile_data import save_profile, write_data
import collections
import html
import json
import datetime
import os
import sys
from pathlib import Path
from qc_definitions import profile as metric_profile, configured_params
from qc_plot_metrics import svg_size_trace
from qc_reference_guide import reference_notes_html
from report_palette import scope_color
from report_brand import brand_bar, copy_mark, favicon
from report_export import add_pdf_export
from report_provenance import RELEASE_VERSION, analysis_record, analysis_label
from report_regions import region_counts as annotation_region_counts, svg_region_bars, load_te_bed, prepare_background
from species_composition import render as render_species
from library_sequences import FAMILIES
from qc_analysis import build_sublibrary_gene_profile, read_numbers, read_scope_profiles, svg_bars, sam_primary_records, base_unit_id, json_dump, read_index_assignments, resolve_chromosome_sets, read_chrom_sizes, write_library_size_assets, write_karyoploter, load_gene_annotation
from report_metrics import COLLISION_ROWS, display_value, fmt, traffic_html, assay_reference_ranges, svg_contact_funnel, METRIC_META, scoped_table_html, REPORT_CSS, svg_contact_evidence, merged_reference_ranges, summary_overview_svg, REPORT_TOOLTIP_SCRIPT, metric_reading_html

def library_structure_html(metric, adapter_summary=None, unit_rows=None, rna_summary=None):
    from report_library_structure import library_structure_html as render_structure
    return render_structure(metric, adapter_summary, unit_rows, adapter_templates=FAMILIES, rna_summary=rna_summary)


def pair_label(pair):
    return f"{pair.get('name',pair.get('scope_id','library'))}: i5 {pair.get('i5') or '∅'} / i7 {pair.get('i7') or '∅'} ({(pair.get('support_fraction') or 0):.1%})"


def identity_panel_html(metric, pairing, pairs, unit_rows=None):
    from report_library_structure import family_label
    families = collections.defaultdict(collections.Counter)
    for row in unit_rows or []:
        if row.get('assignment_status') == 'assigned' and row.get('scope_id') and row.get('family'):
            families[row['scope_id']][row['family']] += 1
    mode=str(pairing.get("mode","unresolved"))
    mode_label=mode.replace("_"," ").title()
    top=pairing.get("top_pair_fraction")
    top_text="NA" if top is None else f"{top:.1%}"
    pair_rows=[]
    sub_ids = {scope['scope_id'] for scope in metric.get('sub_libraries', [])}
    for pair in pairs:
        matched = families.get(pair.get('scope_id'), {})
        adapter = ' / '.join(family_label(name) for name in sorted(matched, key=lambda name: (-matched[name], name))) or 'Unresolved'
        sid = pair.get('scope_id')
        row = f'<tr data-sublibrary="{html.escape(sid, quote=True)}">' if sid in sub_ids else '<tr>'
        pair_rows.append(row + "".join([
            f"<td><strong>{html.escape(str(pair.get('name',pair.get('scope_id','library'))))}</strong></td>",
            f"<td>{html.escape(adapter)}</td>",
            f"<td><code>{html.escape('Not detected' if pair.get('i5') is None else pair['i5'] or 'No index')}</code></td>",
            f"<td><code>{html.escape('Not detected' if pair.get('i7') is None else pair['i7'] or 'No index')}</code></td>",
            f"<td>{(pair.get('support_fraction') or 0):.1%}</td>"]) + "</tr>")
    if not pair_rows:
        pair_rows=["<tr><td colspan='5'>No confident index combination resolved</td></tr>"]
    return ("<div class='identity-panel'><section class='identity-card'><h3>Library identity</h3><div class='identity-grid'>"
            f"<div class='identity-item'><small>Reference assembly</small><strong>{html.escape(str(metric.get('reference_assembly','unspecified')))}</strong></div>"
            f"<div class='identity-item'><small>Assay</small><strong>{html.escape({'chromatin':'Chromatin','microc':'Micro-C'}.get(metric.get('assay'),str(metric.get('assay','unknown'))))}</strong></div>"
            f"<div class='identity-item'><small>Library mode</small><span class='mode-badge'>{html.escape(mode_label)}</span></div>"
            f"<div class='identity-item'><small>Top-pair percentage</small><strong>{html.escape(top_text)}</strong></div></div>"
            "<div class='cutoff-note'>Percentages are among index-informative units."
            f"<details><summary>Library classification</summary><p>A dominant i5–i7 pair is single_library at ≥{pairing.get('single_library_fraction', .9):.0%}, or when it is the only resolved pair. See overall assignment and ambiguous/unassigned fractions below.</p></details></div></section>"
            "<section class='identity-card'><h3>Physical indexes</h3><div class='pair-table-scroll' tabindex='0' aria-label='Physical indexes'><table class='pair-table'><tr><th>Library</th><th>Adapter</th><th>i5</th><th>i7</th><th title='Among index-informative units'>Percentage</th></tr>"
            + "".join(pair_rows) + "</table></div></section></div>")



def profile_scope_html(label, raw_name, insert_name, region_name, karyotype_name, gene_name=None,
                       *, full_library_assets, mitochondrial_html=""):
    def card(name, alt, extra_class=''):
        name = html.escape(name, quote=True)
        alt = html.escape(alt, quote=True)
        return f"<figure class='plot-card {extra_class}'><a href='assets/{name}' target='_blank' rel='noopener'><img src='assets/{name}' alt='{alt}'></a><figcaption class='figure-link'><a href='assets/{name}' target='_blank' rel='noopener'>Open figure ↗</a></figcaption></figure>"
    assets = full_library_assets
    info = assets['model_summary']
    offsets = info['adapter_offsets_bp']
    if len(offsets) == 1:
        definition = f"Full library = insert + {offsets[0]} bp adapter/index sequence."
    elif offsets:
        definition = "Full library = insert + adapter/index sequence; added lengths: " + ', '.join(map(str, offsets)) + " bp."
    else:
        definition = "Full-library length could not be modeled from the available adapter/index information."
    note = definition + f" {info['modeled_units']:,} modeled units."
    if info['excluded_units']:
        note += f" {info['excluded_units']:,} units lack the adapter/index information needed for this model."
    primary = card(assets['weighted'], 'ONT full-library profile (length-weighted)', 'length-weighted-profile')
    counts = card(assets['count'], 'ONT full-library profile (count-based)', 'full-library-count-card')
    counts += card(raw_name, 'Raw molecule size (count-based)') + card(insert_name, 'Insert size (count-based)')
    analysis = card(region_name, 'Genomic regions', 'region-composition') if region_name else ''
    if gene_name:
        analysis += card(gene_name, 'CPM-normalized TSS to TES profile')
    downloads = ' · '.join(f"<a href='assets/{html.escape(assets[key], quote=True)}'>{title}</a>" for key, title in
                           [('profiles', 'Download both profiles'), ('units', 'Modeled unit lengths'), ('model', 'Model details')])
    return (f"<section class='profile-scope'><h3>{html.escape(label)}</h3>"
            f"<p class='muted size-model-note'>{html.escape(note)} The main curve weights each unit by its full-library length.</p>"
            f"<div class='plots size-plots primary-size-plots'>{primary}</div>"
            f"<p class='muted size-profile-downloads'>{downloads}</p>"
            "<details class='count-based-details'><summary>Count-based profiles · full library, raw reads and inserts</summary>"
            "<div class='report-explanation'><ul><li><strong>Full library:</strong> one count per reconstructed unit. Length includes the insert and expected complete adapters/indexes, including undetected outer adapters.</li><li><strong>Raw ONT molecule:</strong> one count per FASTQ read at its observed length.</li><li><strong>Insert:</strong> one count per recovered insert.</li></ul><p>A read with several library units contributes several full-library counts but one raw-read count. Each curve is scaled to its own peak of 100.</p></div>"
            f"<div class='plots size-plots count-size-plots'>{counts}</div></details>" +
            (f"<div class='plots analysis-plots'>{analysis}</div>" if analysis else '') +
            "<details class='chromosome-detail'><summary>Chromosome coverage · auxiliary diagnostic</summary>"
            "<p>Inspect broad coverage imbalance and isolated concentration. Gene density and shallow coverage alone cannot establish target-specific enrichment or copy-number changes.</p>"
            f"<img src='assets/{html.escape(karyotype_name, quote=True)}' alt='Mapped reads by chromosome karyotype'>{mitochondrial_html}</details></section>")


def rna_screening_html(scopes):
    """Keep screening results visible when RNA replaces the generic assay panel."""
    keys = ["mycoplasma_status", "mycoplasma_screened_reads", "mycoplasma_qualifying_reads",
            "mycoplasma_qualifying_fraction", "spike_in_status", "spike_in_screened_reads", "spike_in_reads", "spike_in_fraction"]
    return ("<section class='report-section' id='contamination'><h2>Contamination screening</h2>"
            + scoped_table_html(scopes, keys, {})
            + "<p class='muted'>Not configured = not run. No signal detected applies only to the screened inserts. Foreign-sequence hits may include spike-in.</p></section>")


def ecoli_calibration_note_html(screen):
    """Explain reagent carry-over only for a detected E. coli candidate signal."""
    candidates=screen.get("candidates") or {}
    if not any((key in {"e_coli", "ecoli", "escherichia_coli"}
                or str(candidate.get("label", "")).lower().startswith(("escherichia coli", "e. coli")))
               and candidate.get("status") in {"CANDIDATE_SIGNAL", "TRACE"}
               and (candidate.get("reads") or 0) > 0
               for key,candidate in candidates.items()):
        return ""
    return ("<aside class='note ecoli-calibration-note'><p><strong>E. coli and enzyme carry-over.</strong> "
            "E. coli DNA may carry over from pA/G-MNase (CUT&amp;RUN) or pA/G-Tn5 (CUT&amp;Tag) preparations "
            "and can serve as a spike-in-like calibration reference. Use comparable cell inputs and the same enzyme batch and amount; "
            "confirm sufficient E. coli counts before normalization. Detection alone does not establish the source.</p>"
            "<p>References: <a href='https://doi.org/10.7554/eLife.46314' target='_blank' rel='noopener noreferrer'>"
            "Meers et al., eLife (2019)</a> · "
            "<a href='https://doi.org/10.1038/s41467-019-09982-5' target='_blank' rel='noopener noreferrer'>"
            "Kaya-Okur et al., Nature Communications (2019)</a>.</p></aside>")


def extension_report_html(out, sample, scopes, metric):
    sections=[]
    if metric["assay"] in {"rnaseq", "clip"}:
        align_path = out / "qc/rna_alignment" / f"{sample}.json"
        if align_path.exists():
            evidence = json.loads(align_path.read_text())
            peak = evidence.get("peak_child_rss_kib", 0) / 1024**2
            sections.append(f"<details class='metric-group' id='rna-alignment'><summary>Alignment settings and provenance</summary><p>STARlong {html.escape(str(evidence.get('version', '')))} · maximum insert {evidence.get('maximum_insert_bp', 0):,} bp · peak mapper/sorter RSS {peak:.2f} GiB. BAM sorting uses a separate 256 MiB buffer.</p><p><a href='../qc/rna_alignment/{sample}.json'>Settings and alignment counts</a> · <a href='../qc/rna_preparation/{sample}.json'>Technical preparation evidence</a></p></details>")
    screen_path=out/"qc"/"spike_in"/f"{sample}.json"
    if screen_path.exists():
        screen=json.loads(screen_path.read_text())
        if screen.get("mode") == "auto":
            sublibraries = metric.get("sub_libraries", [])
            if sublibraries:
                from report_spike_matrix import spike_matrix_html
                sections.append(spike_matrix_html(screen, sublibraries, sample, ecoli_calibration_note_html(screen)))
            else:
                rows=[]
                for key,candidate in screen.get("candidates",{}).items():
                    label=html.escape(candidate["label"]);status=candidate["status"]
                    light="amber" if status in {"CANDIDATE_SIGNAL","TRACE"} else "gray"
                    fraction=display_value("spike_in_fraction",candidate.get("fraction"),{})
                    tooltip="Same species as the host; a spike-in cannot be distinguished" if status=="HOST_NOT_DISTINGUISHABLE" else "No spike-in was declared. This is sequence evidence only; a positive signal may have another source."
                    rows.append(f"<tr><td>{label}</td><td>{html.escape(candidate['role'].replace('_',' '))}</td><td>{fmt(candidate.get('reads'))}</td><td>{fraction}</td><td>{fmt(candidate.get('loci'))}</td><td>{metric_reading_html(status,light,tooltip)}</td></tr>")
                declaration="Spike-in was declared in the sample configuration." if screen.get("declared_spike_in") else "No spike-in was declared in the sample configuration."
                tip=declaration+" Candidate signals are exploratory foreign-sequence evidence; they do not establish intentional spike-in. Host, candidate species and microbial decoys compete in one reference. Same-species spike-in and unlisted synthetic standards are not identifiable by this panel."
                scoped_screen=scoped_table_html(scopes,["spike_in_status","spike_in_screened_reads","spike_in_reads","spike_in_fraction"],{}) if len(scopes)>1 else ""
                sections.append("<section class='report-section' id='spike-screen'><h2>Spike-in candidate screen <span class='metric-tip info-icon' tabindex='0' data-tip='"+html.escape(tip,quote=True)+"'>?</span></h2><p class='muted'>"+fmt(screen.get("screened_reads"))+" screened inserts · "+("Spike-in declared" if screen.get("declared_spike_in") else "No spike-in declared")+"</p><div class='table-scroll'><table><tr><th>Reference candidate</th><th>Role</th><th>Insert count</th><th>Fraction</th><th title='Distinct 100-bp reference windows with qualifying alignments'>Locus count</th><th>Result</th></tr>"+''.join(rows)+"</table></div>"+ecoli_calibration_note_html(screen)+scoped_screen+"</section>")

    return ''.join(sections)


def render_sample_report(sample, cfg, out, report_params, prior_provenance, annotation_cache, region_background_cache, args):
    """Render one sample and return its report metrics; caches belong to this run."""
    report_tools=cfg.get("tools", {})
    report_samtools=report_tools.get("samtools", "samtools")
    report_rscript=report_tools.get("rscript", "Rscript")
    karyotype_script=Path(__file__).with_name("karyotype_plot.R")
    promoter_window=int(report_params.get("annotation_promoter_window_bp", 2000))
    report_mapq=int(report_params.get("mapq", 20))
    karyotype_bins=int(report_params.get("karyotype_bins_per_longest_chromosome", 250))
    with open(out/"qc"/"metrics"/f"{sample}.json") as handle: metric=json.load(handle)
    ranges = assay_reference_ranges(cfg, metric["assay"])
    if metric.get("assay") == "chromatin":
        metric["biological_qc_call"] = "not_classified"
        metric["classification_note"] = "Chromatin enrichment is not classified. Collision O/E is depth normalized (random baseline 1). Sampling uncertainty remains at low depth. TSS–TES profiles are CPM normalized and visual only."
    raw=read_numbers(out/"qc"/"lengths"/f"{sample}.raw.tsv"); inserts=read_numbers(out/"qc"/"lengths"/f"{sample}.insert.tsv")
    scope_profiles=read_scope_profiles(out/"qc"/"adapters"/f"{sample}.reads.tsv.gz", include_total=True)
    with open(out/"qc"/"adapters"/f"{sample}.summary.json") as handle:
        adapter_summary=json.load(handle)
    metric["scope_id"]="total"; metric["scope_name"]="Total library"
    metric["scope_unit_fraction"] = 1.0 if metric["library_units"] else None
    scopes=[metric]+metric.get("sub_libraries",[])
    if metric.get("assay") == "chromatin":
        for scope in scopes:
            scope["biological_qc_call"] = "not_classified"
    reference_key=cfg["samples"][sample].get("reference", cfg.get("default_reference", "dm6"))
    report_reference=cfg["references"][reference_key]
    report_assembly=report_reference.get("assembly", reference_key)
    report_gene_bed=out/"reference"/f"{reference_key}.genes.bed"
    report_chrom_sizes=read_chrom_sizes(report_reference["chrom_sizes"])
    primary_arg=report_reference.get("primary_chromosomes", "auto")
    primary_arg=primary_arg if isinstance(primary_arg, str) else ",".join(primary_arg)
    mito_arg=report_reference.get("mitochondrial_chromosomes", "auto")
    mito_arg=mito_arg if isinstance(mito_arg, str) else ",".join(mito_arg)
    report_primary, report_mito=resolve_chromosome_sets(report_chrom_sizes, primary_arg, mito_arg)
    report_annotation=report_reference.get("annotation", "")
    annotation_key=(report_annotation, tuple(report_chrom_sizes.items()), tuple(report_primary), promoter_window)
    report_promoters, report_genes = {}, {}
    if metric["assay"] == "chromatin":
        if annotation_key not in annotation_cache:
            annotation_cache[annotation_key]=load_gene_annotation(report_annotation, report_chrom_sizes, promoter_window)
        report_promoters, report_genes=annotation_cache[annotation_key]
    region_extra={}
    if metric["assay"] == "chromatin":
        te_path=report_reference.get('te_bed')
        background_key=(annotation_key,te_path)
        if background_key not in region_background_cache:
            te=load_te_bed(te_path,report_chrom_sizes)
            region_background_cache[background_key]=(te,prepare_background(report_promoters,report_genes,report_primary,report_chrom_sizes,te))
        te,background=region_background_cache[background_key]
        region_extra=dict(sizes=report_chrom_sizes,te=te,background=background,te_source=report_reference.get('te_source'))
    report_assignments={}
    scope_allowed={}
    if metric.get("sub_libraries"):
        report_assignments=read_index_assignments(out/"qc"/"adapters"/f"{sample}.index_assignments.tsv.gz")
        for sub_scope in metric["sub_libraries"]:
            scope_allowed[sub_scope["scope_id"]]={unit_id for unit_id, row in report_assignments.items() if row.get("scope_id") == sub_scope["scope_id"] and row.get("assignment_status") == "assigned" and row.get("is_adapter_dimer") != "1"}
    report_records=sam_primary_records(str(out/"bam"/f"{sample}.bam"), report_samtools)
    scope_records={"total": report_records}
    for sub_scope in metric.get("sub_libraries", []):
        allowed=scope_allowed.get(sub_scope["scope_id"], set())
        scope_records[sub_scope["scope_id"]]={qname: rec for qname, rec in report_records.items() if base_unit_id(qname) in allowed}
    region_metrics=annotation_region_counts(report_records, report_promoters, report_genes, report_primary, report_mapq, promoter_window, **region_extra) if metric["assay"] == "chromatin" else {}
    sub_region_metrics={}
    for sub_scope in metric.get("sub_libraries", []):
        scoped_records=scope_records[sub_scope["scope_id"]]
        if metric["assay"] == "chromatin":
            sub_region_metrics[sub_scope["scope_id"]]=annotation_region_counts(scoped_records, report_promoters, report_genes, report_primary, report_mapq, promoter_window, **region_extra)
    gene_profile_assets={}
    if metric["assay"] == "chromatin":
        gene_profile_source=out/"plots"/"gene_profile"/f"{sample}.TSS_TES.svg"
        gene_profile_asset=out/"reports"/"assets"/f"{sample}.TSS_TES.svg"
        gene_profile_asset.write_bytes(gene_profile_source.read_bytes())
        gene_profile_assets["total"]=gene_profile_asset
        matrix = out/'qc/gene_profile'/f'{sample}.matrix.gz'
        if matrix.is_file():
            save_profile(matrix, matrix.with_name(f'{sample}.profile.tsv'), matrix.with_name(f'{sample}.profile.json'))
        if metric.get("sub_libraries"):
            genes_path=report_gene_bed
            profile_threads=int(report_params.get("gene_profile_sublibrary_threads", min(16, os.cpu_count() or 2)))
            profile_upstream=int(report_params.get("gene_profile_upstream_bp", 5000))
            profile_downstream=int(report_params.get("gene_profile_downstream_bp", 5000))
            profile_body=int(report_params.get("gene_profile_body_bp", 5000))
            profile_bin=int(report_params.get("gene_profile_bin_size", 100))
            for sub_scope in metric["sub_libraries"]:
                scope_id=sub_scope["scope_id"]
                cached_asset=out/"reports"/"assets"/f"{sample}.{scope_id}.TSS_TES.svg"
                if getattr(args, "reuse_plot_assets", False) and cached_asset.exists():
                    gene_profile_assets[scope_id]=cached_asset
                    continue
                gene_profile_assets[scope_id]=build_sublibrary_gene_profile(
                    out, sample, scope_id, sub_scope.get("scope_name", scope_id),
                    out/"bam"/f"{sample}.bam", genes_path, scope_allowed.get(scope_id, set()),
                    report_samtools, report_tools.get("bedGraphToBigWig", "bedGraphToBigWig"),
                    report_tools.get("computeMatrix", "computeMatrix"), report_tools.get("plotProfile", "plotProfile"),
                    report_mapq, report_reference["chrom_sizes"], report_chrom_sizes, report_primary,
                    profile_upstream, profile_downstream, profile_body, profile_bin, profile_threads)
    profile_sections=[]
    for scope in scopes:
        scope_id=scope["scope_id"]; label=scope.get("scope_name",scope_id)
        profile={"raw_lengths":raw,"insert_lengths":inserts} if scope_id=="total" else scope_profiles.get(scope_id, {"raw_lengths":[],"insert_lengths":[]})
        raw_asset=out/"reports"/"assets"/f"{sample}.{scope_id}.raw_lengths.svg"
        insert_asset=out/"reports"/"assets"/f"{sample}.{scope_id}.insert_lengths.svg"
        region_asset=out/"reports"/"assets"/f"{sample}.{scope_id}.regions.svg"
        karyotype_asset=out/"reports"/"assets"/f"{sample}.{scope_id}.karyotype.svg"
        trace_max=float(report_params.get("size_trace_max_bp", 2000))
        trace_bandwidth=float(report_params.get("size_trace_bandwidth", .045))
        trace_color = scope_color(scope_id) if scope_id != "total" else None
        raw_asset.write_text(svg_size_trace(profile["raw_lengths"],f"{label}: raw ONT molecule size",trace_max,trace_bandwidth, color=trace_color))
        insert_asset.write_text(svg_size_trace(profile["insert_lengths"],f"{label}: adapter-bounded insert size",trace_max,trace_bandwidth, color=trace_color))
        full_library_assets=write_library_size_assets(out/"reports"/"assets", f"{sample}.{scope_id}", label,
            scope_profiles.get(scope_id, {}).get("unit_rows", []), adapter_summary, trace_max, trace_bandwidth, color=trace_color)
        if metric["assay"] == "chromatin":
            region_metric=region_metrics if scope_id == "total" else sub_region_metrics[scope_id]
            region_asset.write_text(svg_region_bars(region_metric, f"{label}: genomic regions", scope_color(scope_id)))
            write_data(region_asset.with_suffix(".json"), "genomic_regions", region_metric)
        else:
            region_asset.unlink(missing_ok=True)
        if not getattr(args, "reuse_plot_assets", False) or not karyotype_asset.exists():
            write_karyoploter(
                scope_records.get(scope_id, {}), report_chrom_sizes, report_primary, report_mapq,
                karyotype_bins, report_rscript, karyotype_script, report_gene_bed,
                karyotype_asset, f"{label}: gene and mapped-read density",
                metric.get("reference_assembly", "custom genome"),
                gene_window_bp=int(report_params.get("karyotype_gene_window_bp", 250000)),
                mitochondrial=report_mito)
        profile_sections.append(profile_scope_html(label, raw_asset.name, insert_asset.name,
                                                   region_asset.name if metric["assay"] == "chromatin" else None, karyotype_asset.name,
                                                   gene_profile_assets[scope_id].name if scope_id in gene_profile_assets else None,
                                                   full_library_assets=full_library_assets,
                                                   mitochondrial_html=""))
    rna_summary = None
    if metric['assay'] in {'rnaseq', 'clip'}:
        rna_summary_path = out / "qc/rna" / sample / "summary.json"
        if rna_summary_path.exists():
            rna_summary = json.loads(rna_summary_path.read_text())
    from report_profiles import render_profiles_card
    microc = None
    if metric["assay"] == "microc":
        from report_microc_layout import collect_ligation_assets, prepare_microc_profiles
        collect_ligation_assets(out, sample, scopes)
        for scope in scopes:
            asset=out/"reports"/"assets"/f"{sample}.{scope['scope_id']}.contact_funnel.svg"
            asset.write_text(svg_contact_funnel(scope,f"{scope.get('scope_name',scope['scope_id'])}: valid-contact funnel"))
            evidence=out/"reports"/"assets"/f"{sample}.{scope['scope_id']}.contact_evidence.svg"
            evidence.write_text(svg_contact_evidence(scope, f"{scope.get('scope_name',scope['scope_id'])}: contact evidence"))
            
            distance_asset=out/"reports"/"assets"/f"{sample}.{scope['scope_id']}.contact_distances.svg"
            distance_asset.write_text(svg_bars(["<1 kb", "1-10 kb", "10-100 kb", ">=100 kb", "trans"],
                [scope.get(k,0) for k in ["cis_lt1kb","cis_1_10kb","cis_10_100kb","cis_ge100kb","trans_contacts"]],
                f"{scope.get('scope_name',scope['scope_id'])}: contact counts (not P(s))"))
        microc = prepare_microc_profiles(sample, scopes, out/'qc', out/'reports/assets', 'assets/', '../qc')
    profiles_html = render_profiles_card(sample, scopes, profile_sections, out/'reports/assets', out/'qc/gene_profile',
                                         rna_summary=rna_summary, microc=microc)
    microc_nav = ''.join(f"<a href='#{key}'>{label}</a>" for key, label in
                         [('pairtools', 'Pairtools'), ('ligation', 'Ligation'), ('mnase', 'MNase')]
                         if microc and key in microc.section_ids)
    pairing=metric.get("index_pairing",{}); pairs=pairing.get("pairs",[])
    if pairing.get("mode")=="single_library" and pairs: pairs=[max(pairs,key=lambda x:x.get("informative_support",0))]
    top=''.join(traffic_html(x, ranges, metric["assay"]) for x in scopes)
    from report_sublibrary import render_sublibrary_composition
    mixture_html = render_sublibrary_composition(metric)
    ont_keys=["raw_reads","sequencing_yield_call","mean_read_q","raw_median_bp","raw_n50_bp","complete_read_rate"]
    library_keys=(["scope_unit_fraction"] if metric.get("sub_libraries") else [])+["outer_complete_units","outer_complete_unit_rate","index_assignment_rate","index_ambiguous_fraction","index_unassigned_fraction","unreconstructed_probable_dimer_reads","adapter_dimer_unit_rate","concatemer_rate","insert_median_bp","library_units","units_per_complete_read","mappable_library_units","canonical_adapter_dimer_units","probable_adapter_dimer_units","adapter_dimer_rate","raw_structure_dimer_rate","adapter_dimer_call","primary_mapping_rate"]
    screen_keys=["mycoplasma_status","mycoplasma_screened_reads","mycoplasma_qualifying_reads","mycoplasma_qualifying_loci","mycoplasma_qualifying_fraction"]
    if metric["assay"]=="chromatin":
        type_specific_keys=[*COLLISION_ROWS,"collision_contributing_molecules","usable_fragments","usable_fragment_rate","mitochondrial_fraction"]+screen_keys
        type_specific_title="Chromatin QC"
    elif metric["assay"] in {"rnaseq", "clip"}:
        type_specific_keys=["usable_fragments","usable_fragment_rate","mitochondrial_fraction"]+screen_keys
        type_specific_title="RNA alignment QC · beta"
    else:
        type_specific_keys=["unique_valid_contacts","cis_ratio","long_range_cis_fraction","valid_contact_rate"]+screen_keys
        type_specific_title="Micro-C QC"
    ont_table=scoped_table_html(scopes,ont_keys,ranges)
    library_key_metrics=(["scope_unit_fraction"] if metric.get("sub_libraries") else [])+["adapter_dimer_call","adapter_dimer_rate","insert_median_bp","outer_complete_unit_rate","index_assignment_rate","primary_mapping_rate","mappable_library_units"]
    library_primary_table=scoped_table_html(scopes,library_key_metrics,ranges)
    library_table=scoped_table_html(scopes,[key for key in library_keys if key not in library_key_metrics],ranges)
    type_specific_table=scoped_table_html(scopes,type_specific_keys,ranges)
    funnel_details = ''
    if metric['assay']=='microc':
        detail_keys=["short_unmapped_units","both_legs_mapped","both_legs_mapping_rate","mapq_pass_pairs","mapq_pass_pair_rate","same_fragment_pairs","same_fragment_filter_pass_pairs","same_fragment_filter_pass_rate","contact_deduplication_retention","cis_contacts","trans_contacts","long_range_cis_contacts","long_range_within_cis_fraction","cis_lt1kb","cis_1_10kb","cis_10_100kb","cis_ge100kb"]
        funnel_details = "<details class='metric-group'><summary>Contact counts and filtering details</summary>"+scoped_table_html(scopes,detail_keys,{})+"</details>"
    metric['analysis_provenance']=analysis_record(metric, prior_provenance)
    definition_label=analysis_label(metric['analysis_provenance'], report_version=RELEASE_VERSION)
    from custom_regions import report_html as custom_region_report
    custom_html=custom_region_report(out,sample,report_assembly,report_chrom_sizes,report_primary,scope_records,scopes,cfg['samples'][sample],Path(__file__).resolve().parents[2],report_mapq) if metric['assay']=='chromatin' else ""
    species_file = out / "qc/species_composition" / f"{sample}.json"
    species_html = render_species(json.loads(species_file.read_text())) if cfg.get("species_composition", {}).get("indexes") else ""
    rna_html = ""
    if metric['assay'] in {'rnaseq', 'clip'}:
        from rna_report import render_rna_section
        rna_html = render_rna_section(out / "qc/rna" / sample / "summary.json", f"../qc/rna/{sample}",
                                      scope_labels={s['scope_id']: s.get('scope_name', s['scope_id']) for s in scopes})
    screen_html = rna_screening_html(scopes) if rna_html else ""
    assay_html = "" if rna_html else (f"<section class='report-section' id='assay-qc'><h2>{type_specific_title}</h2>{type_specific_table}<details><summary>Assay interpretation</summary><p>{html.escape(metric.get('classification_note',''))}</p></details>" +
                funnel_details + "</section>")
    report=f"""<!doctype html><html lang='en'><head><meta charset='utf-8'><meta name='viewport' content='width=device-width,initial-scale=1'><title>{sample} · ONT-LibQC</title>{favicon("assets/ont-libqc-mark.svg")}<style id="report-base-style">{REPORT_CSS}</style></head><body>
{brand_bar("assets/ont-libqc-mark.svg")}
<header class='report-hero'><h1 class='sample-name'>{html.escape(sample)}</h1><div class='hero-meta'><span>{html.escape({'chromatin':'Chromatin','microc':'Micro-C','rnaseq':'RNA-seq · beta','clip':'CLIP-seq · beta'}.get(metric['assay'],metric['assay']))}</span><span>{html.escape(report_assembly)}</span><span>{metric['raw_reads']:,} raw molecules</span><span class="report-version" title="Report release and original analysis date">{html.escape(definition_label)}</span></div></header>
<nav class='report-nav'>{"<a href='#species-composition'>Species</a>" if species_html else ""}<a href='#sequencing'>Sequencing</a><a href='#library-shape'>Library quality</a><a href='#profiles'>{"Size &amp; Micro-C QC" if microc else "Size &amp; profiles"}</a>{"<a href='#rna-qc'>RNA QC · beta</a><a href='#contamination'>Contamination</a>" if rna_html else ""}{"" if rna_html else "<a href='#assay-qc'>Assay metrics</a>"}{"<a href='#custom-regions'>Custom regions</a>" if custom_html else ""}{microc_nav}{"<a href='#spike-screen'>Spike-in screen</a>" if metric.get('spike_in_mode')=='auto' or json.loads((out/'qc/spike_in'/f'{sample}.json').read_text()).get('mode')=='auto' else ''}<a href='../summary/index.html'>All samples</a></nav>
<div class='traffic-grid'>{top}</div>
{mixture_html}
{species_html}
<section class='report-section' id='sequencing'><h2>ONT sequencing</h2>{ont_table}</section>
<section class='report-section' id='library-shape'><h2>Library quality</h2>{identity_panel_html(metric,pairing,pairs,scope_profiles.get("total", {}).get("unit_rows", []))}
{library_structure_html(metric, adapter_summary, scope_profiles.get("total", {}).get("unit_rows", []), rna_summary)}
<div id='library-key-metrics'>{library_primary_table}</div>
<details class='metric-group' id='library-extra-metrics'><summary>All structure and index metrics</summary>{library_table}</details>
{reference_notes_html()}
<details><summary>Size and adapter interpretation</summary><p>The 120–170 bp raw-size window is supportive only. Dimer-like units have a ≤{fmt(metric.get('adapter_dimer_max_bp'))}-bp central gap. Unreconstructed template matches remain candidates; resolved 11–50 bp inserts are retained. Full-library sizes include expected outer adapters and inferred physical indexes, including units with undetected outer handles. The main profile weights each recovered unit by its modeled full-library length; count-based profiles are expandable. Curves are conditional on sequencing and adapter recovery and do not measure TapeStation fluorescence or NGS molarity.</p></details>
</section>
{profiles_html}
{rna_html}
{screen_html}
{assay_html}
{custom_html}{extension_report_html(out, sample, scopes, metric)}
<footer class='report-footer'><details><summary>Methods &amp; interpretation</summary><div class='report-explanation'><ul>{"<li><strong>Chromatin alignment:</strong> Local alignment; ≥30 aligned bp and ≥80% insert coverage. Sequence-supported terminal foldbacks are trimmed only after end-to-end remapping confirms the same locus. MAPQ threshold: ≥"+str(report_mapq)+".</li>" if metric.get('alignment_method') else ''}<li><strong>Reference values</strong> distinguish local screening settings, biological context and statistical expectations.</li><li><strong>Gene profiles:</strong> chromatin TSS–TES profiles use CPM; RNA profiles use per-gene exonic coverage normalization, detailed in the RNA section. SEM is across contributing genes, not biological replicates.</li><li><strong>Genomic distribution:</strong> chromatin uses mapped alignments, Micro-C uses endpoint legs, and RNA region fractions use aligned reference bases excluding skipped introns and deletions.</li><li><strong>Duplicates:</strong> Pairtools marks sample-level coordinate duplicates as DD. Shallow sequencing does not estimate PCR duplication or library complexity. Native Micro-C contact deduplication is an internal validity filter.</li></ul></div></details><p><a href='../summary/index.html'>All samples</a> · <a href='../qc/adapters/{sample}.reads.tsv.gz'>Adapter table</a> · <a href='../qc/adapters/{sample}.index_assignments.tsv.gz'>Index assignments</a></p></footer>{REPORT_TOOLTIP_SCRIPT}</body></html>"""
    if rna_summary is not None:
        from rna_report import color_transcript_assets
        report = color_transcript_assets(report, sample, scopes, rna_summary, out/'qc/rna'/sample,
                                         out/'reports/assets', f'../qc/rna/{sample}', 'assets/')
    (out/"reports"/f"{sample}.html").write_text(add_pdf_export(report))
    return metric


def write_report_summary(metrics, cfg, out, prior_provenance, report_params):
    """Write cross-sample tables, the overview page and provenance."""
    import csv
    columns=["sample","assay","reference_assembly","mean_read_q","outer_complete_unit_rate","index_assignment_rate","index_ambiguous_fraction","index_unassigned_fraction","cis_ratio_ci95_low","cis_ratio_ci95_high","raw_reads","sequencing_yield_call","complete_read_rate","adapter_dimer_rate","adapter_dimer_call","index_architecture","primary_mapping_rate","spike_in_reads","spike_in_fraction","collision_bin_size","collision_oe","collision_fragments","collision_observed_pairs","collision_null_expected_pairs","usable_fragment_rate","both_legs_mapping_rate","mapq_pass_pairs","same_fragment_filter_pass_pairs","unique_valid_contacts","valid_contact_rate","cis_ratio","long_range_cis_min_bp","long_range_cis_contacts","long_range_cis_fraction","contact_evidence_status","metric_profile","metric_definition_version"]
    for m in metrics:
        m.update({k:m.get("metric_settings",{}).get(k) for k in ("metric_profile","metric_definition_version")})
    with open(out/"summary"/"qc_metrics.tsv","w") as handle:
        handle.write('\t'.join(columns)+'\n')
        for m in metrics: handle.write('\t'.join(fmt(m.get(k)) for k in columns)+'\n')
    # Include structure fields in a dedicated numeric table without rounding source counts.
    ligation_columns=["sample","scope_id","scope_name"]+[key for key in METRIC_META if key.startswith("ligation_")]
    with open(out/"summary"/"ligation_qc_metrics.tsv", "w") as handle:
        handle.write("\t".join(ligation_columns)+"\n")
        for m in metrics:
            if m["assay"] != "microc": continue
            for scope in [dict(m,scope_id="total",scope_name="Total library")]+m.get("sub_libraries",[]):
                row={"sample":m["sample"],**scope};handle.write("\t".join(fmt(row.get(k)) for k in ligation_columns)+"\n")
    mnase_rows=[]
    for path in sorted((out/"qc/ligation").glob("*.mnase.json")):
        value=json.loads(path.read_text())
        for sid,stats in [("total",value)]+list(value["scopes"].items()):
            mnase_rows.append(dict(sample=value["sample"],scope_id=sid,**{k:v for k,v in stats.items() if k.startswith("mnase_")}))
    with open(out/"summary/mnase_qc_metrics.tsv","w") as handle:
        writer=csv.DictWriter(handle,list(mnase_rows[0]) if mnase_rows else ["sample","scope_id"],delimiter="\t");writer.writeheader();writer.writerows(mnase_rows)
    scoped_columns = ["sample", "scope_id", "scope_name", "assay", "reference_assembly"] + list(METRIC_META)
    with open(out/"summary"/"sublibrary_qc_metrics.tsv", "w") as handle:
        handle.write("\t".join(scoped_columns) + "\n")
        for m in metrics:
            for sub in m.get("sub_libraries", []):
                row = {"sample": m["sample"], "assay": m["assay"], "reference_assembly": m["reference_assembly"], **sub}
                handle.write("\t".join(fmt(row.get(k)) for k in scoped_columns) + "\n")
    ranges = merged_reference_ranges(cfg)
    summary_plot=out/"plots"/"summary"/"all_samples_qc.svg"; summary_plot.write_text(summary_overview_svg(metrics,ranges,{m["assay"]: assay_reference_ranges(cfg,m["assay"]) for m in metrics}))
    rows=''.join(f"<tr><td><a href='../reports/{m['sample']}.html'>{html.escape(m['sample'])}</a></td><td>{html.escape(m['assay'])}</td><td>{html.escape(str(m.get('reference_assembly','')))}</td><td>{m['raw_reads']:,}</td><td>{display_value('mean_read_q',m.get('mean_read_q'),ranges)}</td><td>{display_value('complete_read_rate',m.get('complete_read_rate'),ranges)}</td><td>{display_value('raw_n50_bp',m.get('raw_n50_bp'),ranges)}</td><td>{display_value('insert_median_bp',m.get('insert_median_bp'),ranges)}</td><td>{display_value('primary_mapping_rate',m.get('primary_mapping_rate'),ranges)}</td><td>{display_value('adapter_dimer_rate',m.get('adapter_dimer_rate'),ranges)}</td></tr>" for m in metrics)
    index=f"""<!doctype html><html lang='en'><head><meta charset='utf-8'><meta name='viewport' content='width=device-width,initial-scale=1'><title>Library reports · ONT-LibQC</title>{favicon('../reports/assets/ont-libqc-mark.svg')}<style id="report-base-style">{REPORT_CSS}</style></head><body>{brand_bar('../reports/assets/ont-libqc-mark.svg', 'Library reports')}<header class='report-hero'><h1>Library reports</h1><p>Sequencing and library quality across samples.</p></header><div class='note'>Open a sample for its complete report, including Chromatin, Micro-C or RNA/CLIP beta metrics.</div><img class='summary-plot' src='../plots/summary/all_samples_qc.svg' alt='All samples technical QC summary'><div class='table-scroll' tabindex='0'><table><tr><th>Sample</th><th>Assay</th><th>Assembly</th><th>Raw reads</th><th>Mean Q</th><th>Complete</th><th>Raw N50</th><th>Insert median</th><th>Mapping</th><th>Adapter dimer</th></tr>{rows}</table></div><p><a href='qc_metrics.tsv'>Download total-library TSV</a> · <a href='sublibrary_qc_metrics.tsv'>Download sublibrary TSV</a></p>{REPORT_TOOLTIP_SCRIPT}</body></html>"""
    (out/"summary"/"index.html").write_text(index)
    json_dump({"created_utc": prior_provenance.get("created_utc") or datetime.datetime.now(datetime.timezone.utc).isoformat(),
               "report_generated_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
               "analyses": {m["sample"]: m["analysis_provenance"] for m in metrics},
               "config": cfg, "python_version": sys.version,
               "metric_schema_version": 3,
               "metric_definitions": {"native": metric_profile(), "matched_comparison": metric_profile("matched_comparison")},
               "library_size_profiles": {"main_weighting": "length", "count_based_default_collapsed": True,
                   "length_definition": "insert + sequence-derived adapters and physical index lengths",
                   "maximum_bp": float(report_params.get("size_trace_max_bp", 2000)),
                   "bandwidth": float(report_params.get("size_trace_bandwidth", .045)),
                   "coordinate": "log1p(bp/50)", "measured_fluorescence": False}}, out/"summary"/"provenance.json")



def cmd_report(args):
    """Coordinate report generation without owning rendering or metric formatting."""
    import yaml
    with open(args.config) as handle:
        cfg = yaml.safe_load(handle)
    out = Path(args.out)
    for directory in ('reports/assets', 'plots/summary', 'summary'):
        (out/directory).mkdir(parents=True, exist_ok=True)
    provenance_path = out/'summary/provenance.json'
    prior = json.loads(provenance_path.read_text()) if provenance_path.exists() else {}
    copy_mark(out/'reports/assets')
    params = configured_params(cfg)
    annotations, backgrounds = {}, {}
    metrics = [render_sample_report(sample, cfg, out, params, prior, annotations, backgrounds, args)
               for sample in cfg['samples']]
    write_report_summary(metrics, cfg, out, prior, params)
