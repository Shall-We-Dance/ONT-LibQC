<p><img src="portal/static/mark.svg" width="40" height="40" alt=""></p>

# ONT-LibQC

[Web server](https://ontlibqc.org) · [Example reports](https://ontlibqc.org/example/index.html) · [Run locally](#run-your-data) · [Self-hosting](portal/README.md)

ONT-LibQC checks completed **Chromatin and Micro-C libraries** using shallow ONT sequencing. It reconstructs adapter-bounded library units, resolves physical i5/i7 index pairs, and reports library size, mapping and assay-specific metrics. Run it in the browser or locally with Snakemake.

## Try it online

Upload FASTQ files at **[ontlibqc.org](https://ontlibqc.org)**, then choose a reference and assay for each library. Email is optional; save the private task link to return to your results. The server accepts up to 16 files per batch, 50 MB each, and keeps reports for 7 days. Do not submit patient data or other sensitive personal data.

Library structure is detected automatically. Use **Customize** for fixed-length UMIs or custom adapters. [Help](https://ontlibqc.org/help) covers settings, input limits and result interpretation.

| Example | Reference | Report contents |
| --- | --- | --- |
| [HL240 - Chromatin](https://ontlibqc.org/example/samples/HL240.html) | hg38 | Library size, gene coverage and custom BED profiles |
| [HL120 - Micro-C](https://ontlibqc.org/example/samples/HL120.html) | dm6 | Pairtools categories, contacts, ligation structure and MNase size evidence |

Use **Load example** on the website to try either FASTQ. The [saved reports](example_report/README.md) also work offline: clone or download this repository and open `example_report/index.html` in a browser.

## Run your data

From the directory containing `Snakefile`:

```bash
conda env create -f environment.yaml
conda activate ont-library-qc
cp config.example.yaml config.yaml
# Edit config.yaml with your FASTQ paths, assays and reference files.
snakemake --cores 8 --rerun-incomplete
```

The default target is `all`. Start with [config.example.yaml](config.example.yaml):

| Setting | What to provide |
| --- | --- |
| `samples.<name>.fastq` | One FASTQ path or a list; plain text or gzip |
| `samples.<name>.type` | `chromatin` or `microc` |
| `references` | FASTA, chromosome sizes, GTF, mapper indexes and primary chromosome names |
| `samples.<name>.library_structure` | Optional custom adapters and fixed-length UMIs; [structure guide](portal/README.md#library-structure-and-fixed-length-umis) |
| `samples.<name>.custom_regions` | Optional BED, assembly, label and plot mode: `auto`, `reference_point` or `region` |
| `spike_in` / `mycoplasma_screen` | Optional competitive screens with prepared reference indexes |
| `params` / `qc_reference_ranges` | Analysis settings and local cutoffs; [metric definitions](workflow/scripts/metric_definitions.json) |
| `tools` | Optional executable overrides; otherwise tools are found on `PATH` |

**BWA-MEM2 is the default Micro-C aligner.** To use classic BWA-MEM, set `mapper: bwa-mem` for the sample and provide a separately built `bwa_index` (`bwa index`). The chosen variant is used for endpoint mapping, full-insert ligation analysis and Pairtools. A portal operator can override the default with a reference's `microc_mapper` setting.

<details>
<summary>Index reuse and memory settings</summary>

- `minimap2_index` reuses a saved `.mmi` built with `-k 9 -w 5`.
- `bwa_batch_bases` bounds classic BWA batches; its default is 1,000,000 bases.
- `samtools_sort_memory` controls memory per sorting thread. The portal uses `128M`.
- Reference preparation and server resource limits are described in the [deployment guide](portal/README.md#reference-and-index-preparation).

</details>

## Read the report

Open `report/index.html` after the workflow finishes. Share the entire `report/` folder, or download **Complete report (ZIP)** from an online result. Extract the ZIP and open `report/index.html`.

| Location | Contents |
| --- | --- |
| `report/index.html` | All-sample overview; online jobs open their sample directly |
| `report/samples/` | Sample reports, including resolved sublibraries |
| `report/assets/` | Figures and figure source tables |
| `report/data/` | QC metrics, adapter/index tables, screening results and analysis settings |
| `report/downloads.html` | Links to all report data |

Chromatin reports accept additional BED files in the browser. Choose a plot mode, confirm the reference coordinates, then select **Draw metaplot**. The submitted BED plot stays available for comparison. Export new plots with their own controls before leaving the page; they are not added to the saved report ZIP.

### Micro-C interpretation

Pairtools analyzes overlapping insert ends of at most 150 bp using fastp and the configured BWA aligner. It complements the native contact and full-insert analyses. All parsed pairs, including DD, WW and Other, remain in the denominator; fastp exclusions are listed separately.

DD marks coordinate duplicates at the observed depth. It does **not** estimate PCR duplication or library complexity. Counts, fractions, settings and per-unit pair types are in `report/data/pairtools/`. See the [Pairtools category definitions](https://pairtools.readthedocs.io/en/latest/parsing.html). Promoter/gene/distal pies are shown only for Chromatin libraries.

### Working files

Calculation intermediates go to `output_dir` (default `.work`). After the report and its links pass validation, the workflow removes this directory, including BAM/SAM, extracted FASTQ, coverage tracks and matrices. Failed runs retain working files for resumption. Original FASTQ and supplied references are preserved.

Use `report_dir` to change the report location. It must be separate from `output_dir`, and neither directory may contain input or reference files. Recomputing a completed analysis requires the original inputs.

## Self-hosting and support

The [portal guide](portal/README.md) covers deployment, reference preparation, resource limits and private report access. Limits are configurable by the operator.

For questions or problems, [open an issue](https://github.com/Shall-We-Dance/ONT-LibQC/issues) or email [hl2522@cornell.edu](mailto:hl2522@cornell.edu).

## License

[MIT License](LICENSE).
