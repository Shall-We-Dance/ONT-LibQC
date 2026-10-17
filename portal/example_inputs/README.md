# Example FASTQ inputs

These research libraries are the complete small inputs used for the [public example reports](https://ontlibqc.org/example/index.html).

| Download | Assay | Reference | Reads | Source |
| --- | --- | --- | ---: | --- |
| [HL120.fastq.gz](https://ontlibqc.org/example-inputs/HL120.fastq.gz) | Micro-C | dm6 | 4,799 | Fruit fly research library |
| [HL240.fastq.gz](https://ontlibqc.org/example-inputs/HL240.fastq.gz) | Chromatin | hg38 | 2,936 | U2OS cell line |

No patient samples are included. Read identifiers, sequences and qualities are retained; acquisition metadata after identifiers and optional plus-line comments are removed. The files have not been downsampled.

**Load example** selects the reference and assay, uses Auto library detection and the server's default mapper, and leaves submission to the user. Updated references or analysis versions may produce results different from the saved reports.

[manifest.json](manifest.json) records download checksums, sizes and settings. Only the named example files are exposed for public download; private uploads are not served from this directory.
