# RNA-seq and CLIP-seq QC · beta

These modes analyze ONT reads of **completed NGS libraries**. The workflow reconstructs adapter-bounded library inserts before splice-aware alignment. An ONT molecule can contain several library units; observed units, genomic alignments and inferred original RNA molecules are different counting units. This release does not implement direct-RNA basecalling or full-length transcript discovery.

Select **RNA-seq · beta** or **CLIP-seq · beta**, provide FASTQ and a prepared reference, and use automatic library-structure detection unless an optional override is needed. Chemistry, selection method, target protein and strandedness are not required inputs. Results describe the evidence available in the reads; a low rRNA fraction, terminal A runs or a 3′-biased profile cannot uniquely identify the original preparation protocol.

## What the report measures

- RNA composition: rRNA, mitochondrial rRNA/mRNA, tRNA, snRNA, snoRNA, protein-coding and other annotated RNA families. Ambiguous family assignments and unassigned reads remain visible. Reference completeness limits these estimates.
- Transcript regions: exon, intron, intergenic, CDS, 5′UTR and 3′UTR coverage. Counts use aligned blocks; skipped introns and deletion gaps do not contribute covered bases or inflate insert length.
- Physical strand: sense, antisense, missing orientation and conflicting annotation, using the recovered P5-to-P7 library direction. Informative denominators and a confidence interval accompany the estimate.
- Splicing: spliced inserts, junction observations and agreement with annotated junctions. Annotation is available to the aligner, so annotated junction support is not independent experimental validation.
- Coverage: exon-concatenated gene-body profiles, a continuous 5′UTR–CDS–3′UTR metaplot, and supported 5′/3′ coverage imbalance. Insufficient coverage is reported explicitly.
- Composition concentration: observed gene counts, dominant genes and depth-dependent detection curves. These are library-unit summaries, not quantitative transcript abundance or deduplicated molecule counts.
- UMI evidence from supported automatic structure detection or a custom configuration: recovery status, exact coordinate/UMI groups and observed collisions. No UMI error correction, deduplication or complexity extrapolation is performed.

The same analysis is available for resolved index-defined sublibraries. Reports retain uncertain or low-depth outcomes instead of assigning a pass/fail result to every metric.

## RNA source reference evidence

A bundled public panel supplies an additional rRNA check even when reads map poorly to the selected host genome. It contains 66 representative mature rRNA sequences from human, mouse, *Escherichia coli*, *Bacillus subtilis*, *Pseudomonas aeruginosa* and *Staphylococcus aureus*. Human and mouse transcripts come from NCBI RefSeq; bacterial sequences are rRNA features extracted using NCBI assembly annotations. The [panel manifest](../workflow/resources/rna_sources/panel.json) records accessions, source URLs, coordinates, sequence hashes and the complete FASTA hash.

Each prepared insert competes across the panel using minimap2. An assignment requires at least 80% query coverage, 85% alignment identity and 30 matching bases. Reference organisms scoring within 95% of the best alignment are retained as ambiguous candidates. Reports preserve the full analysis-input denominator, thresholds and unassigned fraction for each sublibrary. A candidate organism is a reference-sequence match, not a definitive species or contamination diagnosis. This deliberately small panel cannot detect every RNA source, and panel percentages overlap host-genome annotation percentages; they must not be added together. Missing panel data or mapper produce an explicit unavailable result.

## Automatic internal-tag evidence

Within each resolved index-defined sublibrary, the preparer tests the published GRIP/eCLIP X1 structure in physical P5-to-P7 direction: a five-base random tag, fixed `CCTATAT` or `TGCTATT` barcode, RNA-derived insert, and ten-base P7 random tag. Detection requires enough units, repeated exact barcode support and diverse tags at both ends. Only a supported sublibrary is eligible for automatic trimming; individual uncertain units retain their sequence. Successful units record the P5 random tag and fixed barcode separately, recover the P7 N10 UMI with its quality, and remove the supported technical bases before alignment. Barcode evidence in one index does not activate trimming in another. The measured sequence structure is stronger evidence than a protocol name and does not establish the original chemistry by itself.

The physical structure is documented in [Sun et al., Nature Chemistry (2022)](https://doi.org/10.1038/s41557-022-01038-4) and its supplementary oligo table; P7 random-mer UMI handling follows the [ENCODE eCLIP protocol](https://www.encodeproject.org/documents/842f7424-5396-424a-a1a3-3f18707c3222/@@download/attachment/eCLIP_SOP_v1.P_110915.pdf). No random bases are inferred merely from low alignment or end composition.

## CLIP sequence preferences

CLIP motif discovery uses **HOMER 5.1** (`homer2 denovo`) on technically trimmed, unambiguously gene-assigned inserts oriented to the annotated RNA strand. rRNA and terminal soft-clipped query bases are excluded; hard-clipped records are omitted from motif discovery. Identical oriented sequences are collapsed for motif analysis only; the remaining QC counts are unchanged. At least 100 distinct sequences are required. Discovery searches widths 5, 6 and 8 on the forward RNA strand with eight dinucleotide-preserving shuffles per sequence, one CPU, a maximum of 2,000 sequences / 2 million bases, and a 90-second limit per width. An unavailable or failed engine is reported explicitly.

The report shows the leading motif logos, observed sequence support, shuffled-background support and enrichment. Exact consensus matches of at least six bases to a standard adapter (either strand) are flagged as adapter sequence matches, not interpreted as RBP specificity. Native HOMER log10 P values are retained as discovery scores in supporting details and `motifs.tsv`; they are not independent validation or motif-level FDR. Temporary target/background FASTA files are removed after discovery and never enter the downloadable report. FIMO scans supplied motifs and is not used as a de novo discovery replacement. The background does not reproduce IP/input biology or transcript abundance; peak calling, replicate IDR and single-base crosslink calling remain outside this release.

RNA/CLIP reports include an always-visible **Contamination screening** section for Mycoplasma and foreign-sequence screening, including total and sublibrary denominators. Previously disabled screens remain explicitly not run. Species-level foreign-sequence evidence remains in the candidate-screen section.

## STARlong and memory use

The runner builds **STARlong 2.7.11b** explicitly. The regular STAR executable has a different compiled read-length limit; replacing STARlong with ordinary STAR is not equivalent. Variable-length FASTQ sequences and qualities are retained. The RNA wrapper uses local spliced alignment, reduced indel penalties, bounded seeding/multimapping, and explicit minimum aligned support. Its full command and alignment counts are included in the report.

References are prepared with suffix-array sparsity 4, a preindex of at most 12 bases, an adaptive chromosome bin size, and embedded annotation with `sjdbOverhang=100`. These options exchange some speed for reduced memory. The genome-format version recorded by this STAR release is **2.7.4a**, distinct from the executable version.

Online mapping does not retain shared-memory indexes, perform two-pass junction insertion, or rebuild references. The workflow caps STAR at four threads. BAM sorting runs after alignment with a separate 256 MiB sorting buffer. The worker reserves the greater of its configured base memory, operator reservation and index-array sizes plus 3 GiB, rounded upward. Existing total-memory and concurrency limits still apply.

For the validated human reference, the three loaded index arrays occupy about **9.64 GiB**, and the worker reserves **13 GiB**. A synthetic human test with variable 50–10,000 bp queries ran within a 13 GiB container and used about 9.81 GiB peak mapper/sorter RSS. Actual requirements depend on the assembly, annotation, thread count and input; reservation is computed from each prepared index rather than assuming every human reference has the same size.

Source definitions: [STAR indexing parameters](https://github.com/alexdobin/STAR/blob/2.7.11b/source/parametersDefault), [compiled long-read limits](https://github.com/alexdobin/STAR/blob/2.7.11b/source/IncludeDefine.h), and [STAR build targets](https://github.com/alexdobin/STAR/blob/2.7.11b/source/Makefile).

## Prepare references as an operator

Use a catalog whose FASTA, GTF and chromosome sizes describe the same assembly. The script builds indexes locally, checks source hashes and chromosome sizes, writes a preparation manifest, and creates a separate administrator catalog. It does not download references or modify its input catalog.

```bash
python portal/scripts/prepare_star_references.py \
  --catalog /srv/ont-qc/references/catalog.admin.json \
  --output /srv/ont-qc/references/star-2.7.11b-sparse4 \
  --star /path/to/STARlong \
  --threads 8 --build-memory-gib 48
```

Use `--only hg38 mm39` to select references. Build memory is a one-time preparation limit, separate from online mapping memory. Build on a machine with enough available RAM and disk space. Completed indexes are reused only when their input hashes and preparation parameters agree; unfinished indexes remain in a staging directory and are not advertised as ready. Readiness requires the expected executable/index-format versions, sparsity, manifest and complete index files.

After reviewing the generated catalog, merge its `star_index` and `star_alignment_memory_bytes` fields into the deployment catalog and regenerate the public reference catalog. Keep source sequences and indexes outside the web image and mount them read-only for analysis jobs. Rebuild the runner from [Dockerfile.runner](../portal/deploy/Dockerfile.runner) when upgrading a deployment without STARlong.

## Standalone configuration

Add these fields to an existing valid configuration; retain the usual matching FASTA, annotation, chromosome sizes and library-structure settings:

```yaml
references:
  hg38:
    star_index: /path/to/prepared/star/hg38
    star_alignment_memory_bytes: 13958643712  # example; use preparation output
alignment:
  rna_mapper: starlong
samples:
  rna_library:
    type: rnaseq
    reference: hg38
    fastq: [/path/to/rna_library.fastq.gz]
  clip_library:
    type: clip
    reference: hg38
    fastq: [/path/to/clip_library.fastq.gz]
tools:
  starlong: /path/to/STARlong
```

When technical bases inside the extracted insert are known, optional per-sample `rna_trim_p5_bp` and `rna_trim_p7_bp` values remove fixed lengths at the **physical** P5/P7 ends before RNA alignment. Defaults are zero; the trim sides are exchanged for reverse-oriented units and the report records the operation. Nonzero manual trimming or custom `library_structure.umis` extraction disables automatic internal-tag trimming to prevent removing the same bases twice. Do not infer arbitrary linker or UMI lengths from a biological end bias. Custom UMI definitions remain optional.

Detailed evidence is supplied as JSON and downloadable tables for units, gene counts, junctions, profiles, rarefaction, motifs and cohort similarity. Use these denominators and evidence-status fields when comparing shallow libraries.
