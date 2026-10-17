[Portal guide](../README.md) / References

# Reference and index preparation

## Reuse prepared references

Adopt an existing workflow configuration without copying or rebuilding its references:

```bash
/srv/ont-qc/venv/bin/python portal/scripts/prepare_references.py adopt \
  --config /path/to/existing/config.yaml \
  --output /srv/ont-qc/references/common/catalog.admin.json
```

This reads reference settings only. It does not copy samples, raw reads, tools or reports from the source configuration. It checks that the FASTA, annotation, chromosome sizes, and complete Bowtie2 and BWA-MEM2 index sets are present. Optional Mycoplasma and spike-in panels are exposed only when their indexes are ready. Symlinked reference FASTA files retain the correct BWA index prefix.

## Build new references

From the repository root:

```bash
/srv/ont-qc/venv/bin/python portal/scripts/prepare_references.py build \
  --output /srv/ont-qc/references/common --threads 2
```

The operator script downloads the reference FASTA and matching GTF over HTTPS. Docker containers build chromosome sizes and Bowtie2, BWA-MEM2 and Minimap2 indexes with networking disabled. Prepared files live in the output directory outside the Docker images, so they persist across container upgrades and can be mounted read-only by analysis jobs. Starting the portal itself does not download references; initialize it with the generated `catalog.admin.json` after preparation completes.

The [default source list](../deploy/references.common.yaml) pins an annotated RefSeq assembly for each species:

| Species | Assembly | NCBI reference |
| --- | --- | --- |
| Human | GRCh38.p14 | [GCF_000001405.40](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000001405.40/) |
| Mouse | GRCm39 | [GCF_000001635.27](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000001635.27/) |
| Rat | GRCr8 | [GCF_036323735.1](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_036323735.1/) |
| Fruit fly | Release 6 plus ISO1 MT | [GCF_000001215.4](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000001215.4/) |
| C. elegans | WBcel235 | [GCF_000002985.6](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000002985.6/) |
| Zebrafish | GRCz12ab | [GCF_052040795.1](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_052040795.1/) |
| Budding yeast | R64 | [GCF_000146045.2](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000146045.2/) |
| Fission yeast | ASM294v3 | [GCF_000002945.2](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000002945.2/) |
| Arabidopsis | TAIR10.1 | [GCF_000001735.4](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000001735.4/) |
| Rice | AGIS1.0 | [GCF_034140825.1](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_034140825.1/) |

Use the assembly and annotation versions pinned in the manifest. A replacement assembly needs matching annotation to support the full report. Update these sources deliberately; historical alternatives such as hg19 and mm10 are not in the default list.

Use `--only hg38 dm6` to prepare a subset, or `--sources /path/to/sources.yaml` with [the custom template](../deploy/references.example.yaml). FASTA and GTF chromosome names from the common list are normalized together using NCBI assembly reports; generated `chrom.sizes` and the catalog show the names accepted by BED uploads.

## Memory and interrupted builds

Reference preparation is sequential, starting with smaller genomes. It defaults to 2 CPUs, with a memory cap estimated separately for each tool and genome. Bowtie2 uses larger blocks to reduce repeated genome scans while retaining that memory cap. `--memory 96g` is the maximum permitted build allocation, not the allocation for every reference. BWA-MEM2 indexing needs approximately 28 bytes per genome base plus overhead, so a new human-sized reference requires roughly 80–90 GiB during this one-time step. The script checks available memory with 16 GiB reserved by default and fails before launching an undersized build. The portal's analysis limits remain controlled separately by `limits.yaml`.

Completed references are reused. Interrupted builds keep `.building-<reference>` and `build.log`; rerun the same command to resume completed download/index steps. Changed source settings require a new output directory. Each completed reference is recorded in the generated catalog; unfinished references are never marked ready. Do not replace indexes used by active analysis jobs.

To reuse your own already prepared, assembly-compatible references, add `--reuse-catalog /path/to/catalog.admin.json`. Their original annotation and screening settings are retained. Source URLs, file paths and preparation commands are controlled exclusively by the operator. Uploaded files cannot request downloads or builds.

## Minimap2 indexes

Minimap2 indexes use the workflow's existing short-unit settings (`-k 9 -w 5`) and contain the full reference in one index. Online jobs load this prepared index instead of rebuilding it for every submission. Its header is checked for compatible parameters; runtime memory is based on index size plus reserved overhead and remains subject to the 20 GiB online cap. Index preparation and online resource caps are independent. See the [Minimap2 indexing options](https://lh3.github.io/minimap2/minimap2.html) for how saved indexes retain their settings.

The default `build` command prepares all three aligners. To add Minimap2 indexes to an existing administrator catalog:

```bash
/srv/ont-qc/venv/bin/python portal/scripts/prepare_references.py minimap2 \
  --catalog /srv/ont-qc/references/common/catalog.admin.json \
  --output /srv/ont-qc/references/minimap2 --threads 2
```

Use the generated `catalog.admin.json` when initializing or updating the portal. An adopted reference without a matching Minimap2 index still supports its prepared Bowtie2 and BWA-MEM2 modes. The optional `minimap2_index` field in a workflow configuration points to the saved `.mmi` file; standalone CLI runs retain the FASTA fallback when it is omitted.

## Optional screening panels

Host-competitive Mycoplasma and spike-in panels remain operator-configured in the custom source template. A newly prepared genome has these options disabled until a matching panel is provided; the script does not substitute another species' screen.
