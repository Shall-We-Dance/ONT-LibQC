# Sublibrary similarity

Multiplexed chromatin reports with at least two resolved sublibraries include a
**Sublibrary similarity** panel in Size and genomic profiles, before the collapsed
individual profiles. It compares genomic alignment distributions within the same
sample and reference. The pooled total is excluded.

The panel provides Spearman (default) and Pearson correlation, 10 kb / 100 kb
(default) / 1 Mb windows, an average-linkage dendrogram, and downloadable correlation
matrices and sparse bin counts. It works offline in the report ZIP; PDF export
prints the currently selected heatmap.

## Counting and correlation

- Use the existing eligible primary alignments and configured MAPQ threshold
  (normally 20); retain only confidently assigned, non-dimer sublibrary units.
  Secondary, supplementary, duplicate-flagged, QC-failed and alignment-ineligible
  records are excluded by the reporting pipeline.
- Use the same consecutive, non-overlapping primary nuclear genome bins for all
  sublibraries, with chromosome-terminal bins clipped to the reference length.
  Mitochondrial and non-primary sequences are excluded.
- Each alignment counts once in each bin overlapped by its aligned M/=/X blocks.
  Soft clips, deletions and skipped regions do not contribute coverage. Multiple
  blocks in the same bin still contribute one count.
- All genomic bins are retained by default, including shared-zero bins. An optional
  filter excludes bins with zero counts across *all* sublibraries. This can induce
  negative correlation in sparse data, so the filtered view carries a short notice.
- Spearman uses average ranks for ties. Pearson uses raw counts; dividing each
  library's vector by its own total to obtain CPM gives the same coefficient.
  No outlier removal or downsampling is applied.
- Cluster using average linkage on distance `1 - correlation`, with deterministic
  ties. Constant vectors or fewer than two comparison bins yield unavailable
  coefficients, displayed as a dash rather than zero. Undefined libraries are
  placed last and omitted from the dendrogram.

Counts, correlation, reference chromosome selection, eligible alignment counts,
window sizes, filtering mode and clustering are retained in
`assets/<sample>.sublibraries.similarity.json` and the accompanying SVG/TSV files.
The count TSVs are sparse: missing genome windows have all-zero counts. They can
be expanded using the matching reference chromosome sizes. Implicit shared-zero
rows are included exactly in default correlations and tie ranks.

The feature follows the bin-count/correlation workflow of deepTools
[multiBamSummary](https://deeptools.readthedocs.io/en/latest/content/tools/multiBamSummary.html)
and [plotCorrelation](https://deeptools.readthedocs.io/en/latest/content/tools/plotCorrelation.html).
It uses the pipeline's recovered-insert eligibility rules, so it is not a claim of
bit-for-bit equivalence to running deepTools with its default BAM filters.

At shallow depth, sampling and resolution can dominate the result. Correlation is
descriptive and does not establish replicate identity, antibody specificity or
experimental success. This genomic-distribution panel currently applies to
chromatin libraries; it does not substitute for RNA expression correlation or
Micro-C contact-map comparison.
