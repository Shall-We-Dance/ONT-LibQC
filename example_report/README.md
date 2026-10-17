# Example reports

[View online](https://ontlibqc.org/example/index.html), or open `index.html` in a browser after downloading this folder. Keep the folder together so figures and data links work.

| Sample | Assay | Reference | Raw ONT reads |
| --- | --- | --- | ---: |
| [HL240](samples/HL240.html) | Chromatin | hg38 | 2,936 |
| [HL120](samples/HL120.html) | Micro-C | dm6 | 4,799 |

These reports contain real sample measurements, figures and downloadable data. Cohort tables include only these two samples; execution paths in metadata use portable placeholders. [Example FASTQs](../portal/example_inputs/README.md) are available for running your own analysis.

## Chromatin: custom BED profiles

In HL240, choose an hg38 BED file, select **Plot mode**, confirm the coordinates and click **Draw metaplot**. The report offers SVG, data and provenance downloads for the new result.

Selecting another BED clears the previous browser result and requires coordinate confirmation again. Export new results before leaving; they stay in your browser and are not added to the saved report.

## Micro-C: contacts and structure

HL120 includes Pairtools categories, contact metrics, ligation structure and MNase size evidence. Pairtools uses overlapping 150-bp insert ends, fastp and BWA-MEM2. Every parsed pair remains in the denominator, including DD, WW and Other.

Counts, fractions, settings and per-unit pair types are included. DD is a coordinate-duplicate category, not an estimate of PCR duplication or library complexity. Promoter/gene/distal pies apply only to Chromatin reports.
