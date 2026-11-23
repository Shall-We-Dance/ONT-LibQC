# Shallow single-cell RNA QC (Beta)

The local scRNA workflow takes a shallow ONT FASTQ from a finished sequencing library. It identifies supported assay layouts from their barcode positions and adjacent sequences, then measures the yield of gene-assigned molecules. One run represents one FASTQ. Expected cell recovery and paired short-read data are optional experimental context.

## Run a sample

Use the project's RNA environment with Python, Matplotlib, STARlong, samtools and minimap2. Prepare the barcode files once. The downloader verifies the SHA-256 of each file. Subsequent analyses read them locally.

```bash
python workflow/scripts/ont_qc.py prepare-scrna-barcodes --out refs/sc-barcodes

python workflow/scripts/ont_qc.py scrna \
  --fastq sample.fastq.gz --sample sample \
  --assembly danRer11 --gtf refs/danRer11/genes.gtf \
  --star-index refs/danRer11/star --barcodes refs/sc-barcodes \
  --threads 3 --out results/sample
```

Use a STAR index prepared from the same assembly and GTF. The RNA reference builder records both input checksums. The scRNA workflow checks its manifest when available. The GTF digest, chromosome dictionary digest, barcode resource hashes, input digest and STAR version are recorded with the results. See the [RNA reference setup](RNA_BETA.md) and [UCSC zebrafish reference](https://hgdownload.soe.ucsc.edu/goldenPath/danRer11/bigZips/).

## Reconstruction and assay recognition

The standard adapter reconstruction recovers both ONT orientations, splits concatemers, excludes adapter dimers and assigns physical i5/i7 libraries. Unclaimed raw reads are also scanned for ordered MGI Read 1 and Read 2 primer ends. This fallback requires both boundaries and produces a separate unindexed scope. On each unit, the physical Read 1 orientation gives the coordinate system for testing cell barcodes.

| Candidate | Barcode evidence at its expected end | Independent structural evidence |
| --- | --- | --- |
| 10x Chromium 3′ v2, v3/v3.1 or v4 | Read 1 starts with a 16-base barcode from the corresponding on-list | A 10- or 12-base UMI followed by poly(T) |
| 10x Chromium 5′ v1/v2 or v3 | Read 1 starts with the corresponding 16-base barcode | A 10- or 12-base UMI followed by a template-switch sequence |
| MGI DNBelab C4 | Two 10-base barcode rounds, each in the C4 on-list | Linked layout with three intervening or flanking motifs, or the contiguous two-barcode layout, followed by a 10-base UMI and poly(T) |
| Parse Evercode WT, 58-base layout | Read 2 carries a 10-base UMI and three 8-base on-list barcodes at offsets 10, 30 and 50 | Two recurrent 12-base linkers between the rounds |

10x 3′ and 5′ chemistries with a shared barcode list are distinguished by the sequence after the UMI. For Parse, the detector learns the two linkers only from units matching all three barcode lists exactly. Each linker must recur in at least 60% of those independent reads. Both random-primed and poly(T)-primed inserts are eligible. The report names only the layout supported by the evidence. Additional kits can be introduced with a sourced geometry, barcode list and positive and negative controls. [10x library structures](https://kb.10xgenomics.com/s/article/360000939852-What-is-the-difference-between-Single-Cell-3-and-5-Gene-Expression-libraries) · [MGI C4 barcode definitions](https://github.com/MGI-tech-bioinformatics/DNBelab_C_Series_scRNA-analysis-software/blob/2b5c226c189ae24c1dc3ebab242e8cdd65a112f9/config/DNBelabC4_scRNA_readStructure.json) · [MGI read structures](https://mgi-tech-bioinformatics.github.io/DNBelab_C_Series_HT_scRNA-analysis-software/Document/site/en/doc/parameter/scRNA.html) · [Parse Read 2 layout and on-list](https://caltech-bioinformatics-resource-center.github.io/Tutorials/pages/technologies/parse.html).

Each physical library is evaluated independently. A call requires complete barcode and anchor evidence in at least 10 independent raw reads, at least 3 different barcode combinations and at least 5% of the raw reads represented by eligible units in that scope. Exact barcode matches at the expected end must occur at least three times as often as at the opposite end. The report keeps the exact and structure-supported counts for every candidate. Multiple qualifying layouts yield a mixed or ambiguous call. Insufficient or unresolved evidence leaves the scope unassigned. These thresholds govern assay recognition, not cell calling.

## Read and molecule accounting

1. Resolve cell barcodes. For Parse and C4, accept at most one unambiguous, low-quality substitution (base Q ≤20) across all rounds. Parse also permits a one-base indel within each independently supported linker when the best decoding is unique. The 10x branch uses exact on-list barcodes.
2. Accept UMIs with A/C/G/T bases, mean base Q ≥15 and at least two distinct bases. Trim the barcode cassette, anchored poly(T) and supported terminal template-switch handles. An uninterrupted homopolymer at a primer boundary is trimmed as a single run. Retain cDNA of at least 30 bases. Truncated tags with an uncertain boundary are excluded.
3. Align trimmed cDNA to the selected assembly with the shared splice-aware STARlong workflow. Count a gene when there is one primary, uniquely placed MAPQ-pass alignment (default MAPQ 20) and one annotated gene on the expected strand. The report distinguishes multiple-locus, chimeric and unmapped cDNA.
4. Group observations by physical library, barcode combination, annotated gene and exact UMI. Both exonic and intronic evidence contribute to the inclusive count. A key with any exonic evidence also enters the exonic-only count. This gives both views when the cells-versus-nuclei preparation is unknown.
5. Report observed genes, UMIs, barcode ranks, mitochondrial and splice-junction evidence, regional coverage, and rarefaction through the observed sequencing depth. Rarefaction samples the countable gene-assigned units without replacement. Repeated keys measure resampling at the observed depth.

Index-unassigned units can contribute to chemistry and alignment QC. The default expression counts use assigned physical libraries. When the FASTQ is independently known to come from one capture, pass `--single-capture` to pool its matching barcode–gene–UMI keys across index scopes. Conflicting assigned physical libraries or detected layouts are rejected in this mode.

For Parse, one well may have different first-round barcodes for poly(T) and random priming. If a verified, kit-specific mapping is available, `--parse-well-map wells.tsv` takes a complete table with columns `round1_barcode` and `well`. Otherwise the output uses the original three-round combinations. A singleton-heavy rank distribution supports a library QC reading. Estimating recovered cells, doublets, cell types or differential expression requires additional evidence.

## Outputs

**Library identity** lists the reference assembly, assay, library mode and physical sample indexes. The top-pair percentage is the share of index-informative units supporting the most frequent pair. **Sublibrary composition** uses all reconstructed library units, including dimers, and shows the share assigned to each library. Ambiguous and unassigned indexes appear separately when present. Native unindexed units appear as their own group.

Use **Edit names** beside Sublibrary composition or Physical indexes to label an assigned library. The report shows each name alongside its stable library ID. The v1 preview saves names with the shared report. Offline reports save names in the current browser. Expand **Unresolved indexes** to see units with an unassigned or ambiguous index. These units still count toward the whole-library denominator.

**Chemistry detection** reports the fraction of non-dimer units in each sample-index group that match a barcode layout. For Parse, all three cell barcodes must match their lists and both intervening linkers must match at the expected Read 2 positions. An insert can meet those criteria even when its outer sample index is unresolved. The structure diagram uses the same P5-to-P7 orientation and read-start arrows as Chromatin. In the enlarged insert, Read 2 enters from the right and encounters the UMI first. Expand the evidence for the candidate kits and observed Parse linker sequences.

The **Support by library preparation** chart compares Parse Evercode WT, 10x Chromium and MGI DNBelab C4 across all reconstructed non-dimer units. Each unit counts once. Matches to several versions of the same family count toward that family. A unit matching different families goes into **Multiple families**. **No tested layout matched** includes units with incomplete or unrecognized barcode structure. The separate library-level chemistry call also requires support from distinct ONT reads and barcode combinations.

**Parse structure in unmatched units** takes that last group as its denominator. In each index group with a learned Parse linker pair, it checks BC3, linker 1, BC2, linker 2 and BC1 at their expected Read 2 positions, then searches nearby for shifted landmarks. Each unit gets one category. Intact surrounding features can localize a mismatch to a barcode or linker. Two ordered linkers can reveal a shift in position or spacing. Other combinations of landmarks provide partial support when the affected site is unclear. The same search at the opposite physical end provides a comparison.

Expand **Barcode, linker and end evidence** for the per-site counts. A near barcode is one substitution from a unique listed sequence. The nearby linker search allows one substitution, insertion or deletion per linker. Possible truncations include short inserts with surviving landmarks and inserts with barcode-supported landmarks unusually close to the reconstructed Read 2 boundary. A missing cassette prefix, an indel or a shifted reconstruction boundary can produce similar patterns. Where original read coordinates are available, the report also checks how much sequence continues beyond that boundary. These diagnostics do not change barcode decoding or gene and UMI counts.

**Species composition** compares the trimmed, barcode- and UMI-qualified cDNA against a local panel of Bowtie 2 genomes. Supply `--species-panel panel.json` or place `species-panel.json` beside the STAR index directory. Its `indexes` object maps each assembly name to a Bowtie 2 index prefix. Paths are relative to the panel file or absolute. If the host and panel assemblies have different names, map the host name to the corresponding same-species panel assembly in `assembly_aliases`. At least two genomes are required. `--bowtie2` selects the executable. The report includes all qualified cDNA, each physical index scope and a comparison excluding identified rRNA and mt-rRNA. Genomic local alignment can miss fragments spanning splice junctions.

**Transcript regions** uses uniquely aligned cDNA passing the configured MAPQ threshold. The figure divides the bases assigned to 5′ UTR, CDS or 3′ UTR among those three categories. The expanded table includes other exonic, intronic and intergenic bases and gives percentages of all aligned bases. The calculation uses matched CIGAR blocks, transcript strand and explicit or inferred UTRs. When transcripts overlap, each aligned base is counted once.

The **transcript metaplot** follows the spliced representative transcript from 5′ to 3′. The 5′ UTR, CDS and 3′ UTR each have 30 bins. Coverage in each bin is divided by its exonic width and by that gene’s mean coverage across its representative transcript. Genes are then averaged with equal weight. The band shows the standard error across genes. An annotated region with zero coverage contributes zeros. A gene missing the region is omitted from that region’s average. The descriptive curve includes genes supported by a single cDNA fragment.

The RNA family plot compares all barcode- and UMI-qualified cDNA with the multimapped subset. Each fragment contributes once per group, including unclassified fragments. The RNA-seq family classifier checks all reported genomic alignments. A qualifying match to the independent rRNA panel takes precedence. The gene and UMI matrix uses the original unique-alignment criteria.

Use `--rna-family-gtf` for a same-assembly annotation with explicit gene biotypes, and `--rna-source-panel` for an extended RNA-seq panel containing `panel.fa` and its checksum manifest `panel.json`. By default the workflow uses the STAR GTF and the bundled RNA-seq panel. A same-assembly UCSC table supplied with `--repeatmasker rmsk.txt.gz` adds repeat evidence. Each candidate locus is evaluated over its aligned bases. At least 50% coverage by rRNA repeats at every reported locus supports a repeat-based rRNA assignment. The evidence table shows partial and non-rRNA repeat overlaps separately.

| File | Contents |
| --- | --- |
| `report.html` | Self-contained report with aggregate chemistry evidence, stage denominators and figures |
| `summary.json` | Aggregate metrics, rank counts, rarefaction points and reference provenance |
| `rna_families.svg` | RNA family percentages for runs with qualified cDNA |
| `transcript_regions.svg` | Proportions and base counts for 5′ UTR, CDS and 3′ UTR when at least one base is assigned |
| `transcript_metaplot.svg` | Observed relative coverage along the representative transcript when coverage is available |
| `barcode_gene_counts.tsv.gz` | Sparse long-format matrix with physical library, barcode, gene ID, inclusive UMIs and exonic UMIs |

`--keep-work` also retains the private reconstruction tables, trimmed cDNA FASTQ and BAM. The reporting directory is published only after all steps succeed. Use a fresh output path for each run.
