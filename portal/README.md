[ONT-LibQC](../README.md) / Web portal

# Web portal

A self-hosted interface for uploading FASTQ files, running ONT-LibQC and sharing private reports. Each file has its own settings and task link. Email notifications are optional.

**[Public server](https://ontlibqc.org)** · [Help](https://ontlibqc.org/help) · [Example reports](https://ontlibqc.org/example/index.html)

## Submit a library

1. Add `.fastq`, `.fq`, `.fastq.gz` or `.fq.gz` files. Basecall POD5/FAST5 files first.
2. Choose a **Reference genome** and **Library assay** for each file. Neither is preselected. With several files, **Apply to all** copies these two settings; other settings and sample names stay independent.
3. Leave the library structure on **Auto**, or customize adapters and fixed-length UMIs below.
4. Add an email address if you want notification, accept the data declarations and submit. Save the private task links.

Mycoplasma and foreign-sequence screens start enabled when the selected reference has a prepared panel. Users can switch them off; unavailable panels stay disabled. Declaring an intentional spike-in is a separate choice.

**Do not upload human patient data, including de-identified patient data, or other sensitive personal data.** This applies to reads, BED files, filenames and sample names. The upload form requires separate permission and data-policy acknowledgments; these are user declarations, not automatic data screening.

## Library structure and fixed-length UMIs

The diagram follows one physical strand:

```text
5′ P5 - i5 - adapter - insert - adapter - i7 - P7 3′
```

**Auto** detects P5/P7 handles, inner adapters and physical i5/i7 indexes. **Customize** supports:

| Control | Behavior |
| --- | --- |
| Add UMI | Click a gap, then set a fixed length of 1–64 bp |
| Move or remove UMI | Use each row's controls, or drag a block/row handle to a gap; dropping on another UMI swaps positions |
| Adapter preset | Select TruSeq / Nextera, or keep automatic detection |
| Custom sequences | Set P5/P7 handles and either inner adapter independently; each accepts 16–80 A/C/G/T bases |
| Index lengths | Set expected lengths from 0–24 bp, with a two-base ONT indel tolerance |
| Apply structure to all | Copy the complete structure settings to other files in the submission |

Every configured UMI remains visible in its own row. Supported positions are `before_i5`, `after_i5`, `before_insert`, `after_insert`, `before_i7` and `after_i7`. Returning to Auto retains the editor draft but submits no custom overrides.

### Sequence orientation

Read 1/2 are directions, not physical sequence blocks. Short arrows start at the insert end, or at an adjacent UMI between the insert and adapter. Primers bind the adapters; their binding sequences are not read at the start of Read 1/2.

Custom sequences follow the displayed **5′ P5-to-P7 strand**. The standard P7 end is `ATCTCGTATGCCGTCTTCTGCTTG`, the reverse complement of the usual P7 oligo. Likewise, `adapter_i7_sequence` is not the opposite-strand primer sequence.

Set inner adapters through **Adapter & index settings**, or click an Adapter block. A filled `adapter_i5_sequence` or `adapter_i7_sequence` overrides that side of the preset. A blank side uses the preset, or both supported presets when Auto is selected. These sequences affect reconstruction, dimer templates, full-library length modeling and report diagrams.

### UMI interpretation

UMIs are removed from indexes and mapped genomic inserts, recorded individually and together in the per-unit table, and included in modeled full-library lengths. Reports show the layout and extraction counts for the total library and resolved sublibraries.

- Missing outer anchors leave the corresponding UMI unknown. Gaps shorter than the configured layout are rejected.
- Fixed-offset extraction does not correct ONT indels within random UMIs.
- Assay counting and deduplication are unchanged; no UMI correction or UMI-based deduplication is performed.
- Random-UMI layouts cannot use the unreconstructed-read zero-insert template screen. Reconstructed insert-length dimer calls still apply.

The report diagram shows i5/i7 sequences, adapter types and UMI lengths. Expand **P5/P7 & adapter sequences** for full motifs and evidence. Matched motifs are distinguished from configured but undetected custom sequences.

## Read and save results

Open a report through its complete private task link, including the part after `#`. Viewing and downloads require no cookies or browser storage. Short-lived report/figure links are renewed from the task page; if one expires, reopen the original private task link. That link works until deletion or retention expiry.

Online jobs open their sample directly. **Download complete report (ZIP)** includes saved HTML, figures, tables and settings. Extract it and open `report/index.html`; keep the folder together.

### Custom BED profiles

Chromatin reports keep the BED plot submitted with the analysis and accept additional BEDs in the browser. Choose **Plot mode** (`Auto`, `Reference point` or `Region (scaled)`), confirm the assembly coordinates and click **Draw metaplot**.

Browser profiles accept up to 100,000 retained regions and a 10 MiB BED. Directory and unreadable selections are rejected. Selecting another file resets coordinate confirmation and clears the previous browser result. New reports retain the original BED filename; older reports may lack it.

Export new browser-generated plots with their own controls before leaving. They are not added to the server's saved ZIP.

### Progress

Upload progress tracks bytes sent. Analysis progress tracks completed workflow steps, which have different runtimes; the percentage does not estimate remaining time. Updates survive page refreshes. The UI respects **Reduce Motion** and uses no looping decorative animation.

The worker samples progress about every 3 seconds. The browser polls every 10 seconds (30 seconds when hidden), with jitter and rate-limit backoff. In split deployments, bridge transfer time and its 5-second cycle can delay updates. Computation continues through SSH interruptions; reports appear only after the complete bundle is received and validated.

## Public examples and help

**Load example** selects HL120 (dm6 / Micro-C) or HL240 (hg38 / Chromatin), fills the reference and assay, and keeps Auto structure detection. It does not submit a job or accept the data declarations.

The [example FASTQs](example_inputs/README.md) retain read IDs, sequences and qualities, with acquisition metadata removed. Only these named downloads are exposed at `/example-inputs/`.

The public `/help` page covers formats, settings, interpretation, limits and common errors. `/license` and `/LICENSE` provide the full MIT License. These pages and the example reports require no authentication.

## Host your own server

| Guide | Contents |
| --- | --- |
| [Deployment](docs/deployment.md) | Installation, HTTPS, renewal, migration and separate compute hosts |
| [Reference and index preparation](docs/references.md) | Reuse or build references, memory needs and optional panels |
| [Security and limits](docs/security.md) | Isolation, private access, request pools, quotas and retention |
| [Email notifications](docs/email.md) | Optional SMTP delivery, batch summaries, retries and sender authentication |

Operators control reference paths, executable names, downloads and resource limits. Users cannot supply arbitrary server paths or commands.

## Reference and index preparation

Prepare or adopt references before starting the server. The [reference guide](docs/references.md) lists commands, supported assemblies and index requirements. BWA-MEM2 remains the default Micro-C aligner; classic BWA-MEM is optional when separately configured.

## Email notification

An address is optional. If email is enabled and an address is supplied, one summary is sent after all files in the submission finish. Otherwise, users access results through saved private task links. See [SMTP setup and delivery behavior](docs/email.md).

## Isolation and data handling

The public web service is separate from isolated analysis containers. Private task credentials authorize access; reads and working files are cleaned after processing, and reports expire after the configured retention period. See [security, limits and retention](docs/security.md).
