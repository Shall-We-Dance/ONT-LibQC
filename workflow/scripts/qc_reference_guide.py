"""Report-only reference context. Advisory expectations never create QC cutoffs."""
import html

SOURCES={
    'cutrun':('Skene & Henikoff, 2017','https://elifesciences.org/articles/21856'),
    'microc':('Mapping Mammalian 3D Genome Interactions with Micro-C-XL','https://www.jove.com/t/64579/mapping-mammalian-3d-genome-interactions-with-micro-c-xl'),
}

# text, basis, explanation; these are annotation, not statistical reference intervals.
GUIDE={
 'mean_read_q':('Chemistry-specific','Measurement context','No validated ONT mean-read Q cutoff is assigned. Phred Q20 corresponds to an estimated 1% error probability, but chemistry, basecaller calibration and the averaging method matter.'),
 'raw_median_bp':('Library-dependent','Structure context','Raw ONT molecule length includes recovered library and adapter structure. Concatemers and incomplete molecules change this distribution; it is not the NGS insert length.'),
 'raw_n50_bp':('Library-dependent','Structure context','N50 measures the raw molecule-size distribution. A larger value is not inherently better for a finished short-read library.'),
 'insert_median_bp':('Target / protocol-specific','Biological context','CUT&RUN can recover both short protein-associated fragments and nucleosome-sized fragments. The expected mixture depends on target and protocol; no universal median-size interval is assigned.'),
 'outer_complete_unit_rate':('Both P5 and P7 expected','Library design','The expected P5/P7 library architecture includes both outer handles. ONT read truncation, errors and incomplete sequence recovery can reduce observed handle detection, so 100% is a design expectation rather than a calibrated ONT pass threshold.'),
 'index_assignment_rate':('Prefer high assignment','Library identity','High all-unit assignment is desirable, but no validated universal percentage is available for inferred ONT indexes. Assess against the known index sheet and examine ambiguous/unassigned units.'),
 'index_ambiguous_fraction':('Prefer 0%','Library identity','Ambiguous assignments should be minimized. Nearby indexes, ONT errors and missing sequence can increase the fraction; no new color threshold is implied.'),
 'index_unassigned_fraction':('Prefer 0%','Library identity','Unassigned units reduce confidence in sublibrary attribution. The achievable rate depends on index structure and sequence recovery; no new color threshold is implied.'),
 'scope_unit_fraction':('Expected pooling ratio','Experimental design','Compare with the intended library mixture when known. Equal shares are not assumed; barcode assignment and ONT recovery may distort proportions.'),
 'concatemer_rate':('Preparation-dependent','Structure context','Multiple adapter-bounded units within one ONT molecule are not the same as multiple Micro-C genomic fragments ligated within one unit. No universal desirable concatemer fraction is established.'),
 'units_per_complete_read':('≥1 by definition','Metric definition','At least one unit defines a reconstructed read. Values above one describe ONT molecules containing multiple units, not a biological success threshold.'),
 'adapter_dimer_call':('≥20 units; see burden','Local screening rule','The configured raw-molecule dimer-burden thresholds and a minimum of 20 reconstructed units determine the call. This is an ONT screening rule, not a calibrated NGS loading assessment.'),
 'adapter_dimer_unit_rate':('Prefer 0%','Library purity objective','This denominator is reconstructed units. The configured adapter-dimer call instead uses raw molecules containing at least one reconstructed dimer-like unit; these percentages are not interchangeable.'),
 'raw_structure_dimer_rate':('Prefer 0%','Library purity objective','Minimize evidence of zero-insert library structures. Raw size alone does not establish an adapter dimer, and unresolved template candidates are not included in confirmed burden.'),
 'sequencing_yield_call':('Local yield screen','Local screening rule','Uses minimum raw-read and reconstructed-unit requirements from the effective pipeline configuration. It is not a biological success/failure label.'),
 'mycoplasma_status':('NOT_DETECTED','Investigator reporting policy','Zero qualifying Mycoplasma signal is green; any supporting signal is red. A negative shallow screen does not establish that the library is contamination-free. Zero screened inserts is unavailable.'),
 'mycoplasma_qualifying_reads':('0','Investigator reporting policy','Desired supporting-read count is zero. Any qualifying target read triggers a red indicator; no-screen data are not a measured zero.'),
 'mycoplasma_qualifying_loci':('0','Investigator reporting policy','Desired supported target-locus count is zero. Any positive signal triggers a red indicator, while its specificity and depth remain available in screen details.'),
 'mycoplasma_qualifying_fraction':('0%','Investigator reporting policy','Zero qualifying signal is green and any nonzero signal is red. The denominator is screened adapter-trimmed inserts, not original culture cell abundance.'),
 'mycoplasma_screened_reads':('Depth-dependent','Sampling context','More screened inserts improve sensitivity. No read-count target can guarantee absence of Mycoplasma; the negative detection limit depends on the screening design and recovery bias.'),
 'spike_in_status':('Experimental record','Experimental design','Reference candidates are assessed against declared spike-in metadata. Sequence evidence alone does not identify deliberate addition.'),
 'spike_in_fraction':('Experimental record','Experimental design','There is no universal spike-in percentage. Added amounts, species, extraction and alignment determine recovery. Candidate matches do not establish intentional spike-in.'),
 'cis_ratio':('No biological cutoff','Descriptive contact composition','Cis / all unique valid contacts. The interval assumes independent contacts; count sufficiency does not establish biological success.'),
 'collision_contributing_molecules':('Evidence count','Sampling context','Distinct original ONT molecules contributing eligible chromatin fragments. Several units can share one parent molecule.'),
 'collision_oe':('1 = uniform-null baseline','Statistical expectation','The uniform-base null has O/E = 1 with terminal-bin-width correction. O/E >1 indicates excess concentration relative to that null; mappability, repeats, copy number and target biology still matter. This is not an enrichment pass cutoff.'),
 'collision_bin_size':('Fixed analysis setting','Metric definition','The active analysis uses 10 kb bins. Changing bin size changes the spatial scale and may change O/E; this is not a sample-quality criterion.'),
 'collision_observed_pairs':('Depth / concentration','Sampling context','The count depends on the number and genomic distribution of eligible alignments. Sparse or zero observations do not establish failed enrichment.'),
 'collision_null_expected_pairs':('N choose 2 × p₀','Statistical expectation','The expected same-bin pair count is the number of eligible alignment pairs times the terminal-bin-corrected uniform-base collision probability.'),
 'long_range_cis_min_bp':('Analysis definition','Metric definition','This is the configured minimum separation defining long-range cis, currently 10 kb, rather than a sample-quality interval.'),
 'cis_ratio_ci95_low':('Lower uncertainty bound','Statistical interval','Interpret together with the cis estimate, upper bound and unique-contact denominator. The Wilson interval assumes independent contacts and does not include recovery bias.'),
 'cis_ratio_ci95_high':('Upper uncertainty bound','Statistical interval','Interpret together with the cis estimate, lower bound and unique-contact denominator. No independent pass threshold is assigned to this endpoint.'),
 'long_range_cis_ci95_low':('Lower uncertainty bound','Statistical interval','Uncertainty for the long-range cis fraction. No universal bound is validated for this shallow ONT endpoint assay.'),
 'long_range_cis_ci95_high':('Upper uncertainty bound','Statistical interval','Uncertainty for the long-range cis fraction. The interval depends on the number of qualifying contacts.'),
 'ligation_resolved_fraction':('No validated cutoff','Alignment context','Resolution depends on block length, uniqueness, query coverage and alignment quality. A high resolved fraction does not alone prove successful ligation.'),
 'ligation_multifragment_fraction':('Protocol-dependent','Structure context','Evidence for multiple genomic blocks is conditional on library recovery and alignment resolution. No calibrated minimum ligation-efficiency fraction is established.'),
 'ligation_multifragment_fraction_resolved':('Protocol-dependent','Structure context','This fraction excludes unresolved units from its denominator. It is not absolute ligation efficiency.'),
 'same_fragment_filter_pass_rate':('Protocol-dependent','Contact filtering context','The short inward filter can remove real nearby contacts. Retention is descriptive and no universal minimum is imposed.'),
 'mapq_pass_pair_rate':('Library / reference-specific','Alignment context','Retention depends on read length, repeats, reference choice and mapper confidence. No calibrated universal pair-MAPQ retention rate is available.'),
 'valid_contact_rate':('Protocol-dependent','Contact filtering context','Uses unique endpoint contacts per mappable unit. Thresholds from deep NGS Hi-C or different fragment filtering rules cannot be transferred without validation.'),
 'contact_deduplication_retention':('Depth-dependent','Sampling context','Coordinate coincidence depends on sampling depth and local enrichment. At shallow depth this is not a reliable estimate of library complexity.'),
 'long_range_cis_fraction':('Protocol / cell-specific','Biological context','Depends on folding, cell type, fragment recovery and distance definitions. Compare with matched controls; no universal percentage is validated for this ONT assay.'),
 'long_range_within_cis_fraction':('Protocol / cell-specific','Biological context','Uses cis contacts only as denominator. Compare with matched controls under the same distance cutoff and filtering rules.'),
}

COUNT_KEYS={'raw_reads','library_units','mappable_library_units','outer_complete_units','usable_fragments','collision_fragments','both_legs_mapped','mapq_pass_pairs','same_fragment_filter_pass_pairs','unique_valid_contacts','long_range_cis_contacts','cis_lt1kb','cis_1_10kb','cis_10_100kb','cis_ge100kb','spike_in_screened_reads','spike_in_reads'}


def context(key,metric=None):
    if key=='insert_median_bp' and (metric or {}).get('assay')=='microc':
        return ('Protocol-specific','Biological context','Two approximately 150-bp nucleosomal fragments can yield an approximately 300-bp insert before sequencing adapters. Final size selection and higher-order ligation change the distribution; this is context, not a fixed median-size pass range.')
    if key in GUIDE:return GUIDE[key]
    if key in COUNT_KEYS or key.startswith('ligation_'):
        return ('Depth-dependent','Sampling context','Absolute counts depend on recovered sequencing depth and the eligibility filters. More observations improve estimation, but no universal sample-quality count is assigned here.')
    if key in {'canonical_adapter_dimer_units','probable_adapter_dimer_units','unreconstructed_probable_dimer_reads'}:
        return ('Prefer 0; assess fraction','Library purity objective','Absolute dimer/candidate counts depend on depth. Use raw-molecule dimer burden for screening; unresolved candidates remain separate.')
    if key=='short_unmapped_units':return ('Inspect size distribution','Alignment eligibility','Micro-C units shorter than 50 bp cannot supply two nonoverlapping 25-bp legs. This count is descriptive.')
    if key=='same_fragment_pairs':return ('Protocol-dependent','Contact filtering context','A short inward orientation is a heuristic exclusion, not proof that both ends came from the same original fragment.')
    return ('No universal cutoff','Descriptive context','No validated universal reference interval is available for this metric in shallow ONT library QC.')


def reference_notes_html():
    link=lambda k:f"<a href='{SOURCES[k][1]}' target='_blank' rel='noopener'>{html.escape(SOURCES[k][0])}</a>"
    return ("<details class='reference-notes'><summary>Reference intervals: basis and sources</summary>"
            "<p>Hover a reference value for its basis. Numeric colored ranges are local screening settings unless stated otherwise. Biological expectations, statistical baselines and design objectives provide context and do not introduce new pass/fail thresholds.</p>"
            "<p><strong>Adapter dimers.</strong> ONT raw-molecule dimer burden is a count-based screening measurement. It is not interchangeable with electrophoresis mass fraction or platform-specific loading limits. The configured thresholds are local screening settings.</p>"
            "<p><strong>Library handles.</strong> The expected library architecture includes both P5 and P7. Observed ONT handle recovery also depends on sequencing and reconstruction.</p>"
            "<p><strong>Fragment sizes.</strong> CUT&amp;RUN has target-dependent short and nucleosomal components. Micro-C MNase titration evaluates approximately 150-bp nucleosomal DNA and longer products before ligation. Neither provides a universal acceptable insert median or a calibrated threshold for the inferred fragments in this report. "+link('cutrun')+"; "+link('microc')+".</p></details>")
