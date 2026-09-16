<p><img src="assets/ont-libqc-mark.svg" width="40" height="40" alt="ONT-LibQC"></p>

# ONT-LibQC · Example reports

Open `index.html` directly in a browser. Keep this entire folder together.

- [HL240](samples/HL240.html): chromatin library, hg38, 2,936 raw ONT molecules.
- [HL120](samples/HL120.html): Micro-C library, dm6, 4,799 raw ONT molecules.

These are real sample reports from the ONT-LibQC analysis. Sample measurements and figures are retained. The entry page links directly to these two reports; downloadable cohort tables contain only these two samples.

HL240 supports local analysis of hg38 BED regions: choose a BED, select **Plot mode** (`Auto`, `Reference point` or `Region (scaled)`), confirm hg38 coordinates and click **Draw metaplot**. The result shows the selected filename and provides SVG, data and provenance downloads. Selecting another BED clears the previous browser result and requires coordinate confirmation again. These interactive results stay in your browser; export them before leaving the page.

HL120 includes Pairtools category fractions after duplicate marking, contact, ligation-structure and MNase-size evidence. Its promoter/gene/distal pie has been removed. The Pairtools analysis uses overlapping 150-bp insert ends, fastp and BWA-MEM2; every parsed pair stays in the denominator, including DD, WW and Other. Counts, fractions, settings and per-unit types are downloadable. All linked figures and report data are included. Execution paths in metadata use portable placeholders.
