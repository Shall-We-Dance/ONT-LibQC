# ONT library QC

Snakemake workflow for shallow ONT QC of completed Illumina chromatin and Micro-C libraries. It reconstructs adapter-bounded library units, identifies physical index pairs, and produces an interactive HTML report.

## Example report

**[Open the example report](example_report/index.html)**

Download or clone the repository, then open `example_report/index.html` in your browser. The report is already generated and works locally. GitHub's file view displays the HTML source.

The example uses synthetic DNA and includes a chromatin sample with two sublibraries and a Micro-C sample. Explore library size profiles, gene coverage, custom BED metaplots, contact metrics, ligation structure and MNase size evidence. These are demonstration data. The [example BED](example_report/assets/chromatin.custom_regions.bed) uses the synthetic reference coordinates and can be uploaded in the chromatin report.

## Run your data

```bash
conda env create -f environment.yaml
conda activate ont-library-qc
cp config.example.yaml config.yaml
# Edit FASTQ paths, assay types, references and optional screens.
snakemake --cores 8 --rerun-incomplete
```

Run commands from the directory containing `Snakefile`. The default target is `all`.

- `samples.<name>.fastq`: one FASTQ path or a list, plain or gzip compressed.
- `samples.<name>.type`: `chromatin` or `microc`.
- `references`: FASTA, chromosome sizes, GTF annotation, prebuilt mapper indexes and primary chromosome names.
- `samples.<name>.custom_regions`: optional BED, assembly and label. Chromatin reports also accept a BED in the browser.
- `spike_in` / `mycoplasma_screen`: optional [competitive screens](docs/screens.md).
- `params` / `qc_reference_ranges`: analysis settings and local screening cutoffs; see [metric definitions](workflow/scripts/metric_definitions.json).
- `tools`: optional executable overrides; otherwise commands are found on `PATH`.

## Output

Open `report/index.html` after the run. Copy the entire `report/` folder to share it. BED uploads and plot exports run locally in the browser.

| Location | Contents |
| --- | --- |
| `report/index.html` | All-sample overview |
| `report/samples/` | Sample reports and sublibrary results |
| `report/assets/` | Figures and figure source tables |
| `report/data/` | QC metrics, adapter/index tables, screens, ligation and MNase evidence |
| `report/downloads.html` | All downloadable report data |

All calculation intermediates go to `output_dir` (default `.work`). After the report is generated and its links are checked, the workflow removes BAM/BAI, SAM, extracted FASTQ, coverage tracks, matrices and the working directory. Failed runs retain working files for resumption. Original FASTQ and supplied references are preserved.

Set `report_dir` to change the report location. It must be separate from `output_dir`; both directories must be outside the input/reference paths. Recomputing a completed analysis requires the original inputs because calculation intermediates have been removed.