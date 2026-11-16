# Refresh retained report presentation

`workflow/scripts/report_refresh.py` updates retained standalone report bundles
without FASTQ, BAM, reference genomes or mapping tools. It regenerates the size
and genomic profiles card from retained numerical data, updates the main report
styles and tooltips, and removes the withdrawn RNA cohort panel. Other sections,
analysis dates, versions, metrics and provenance remain intact.

## Inventory and preview

Keep the inventory outside the repository. It contains private deployment paths.
Use a bundle root containing `manifest.json`, `samples/`, `assets/` and `data/`:

```json
{
  "schema_version": 1,
  "reports": [
    {"root": "/path/to/retained-report"}
  ]
}
```

An optional `manifest_sha256` on each entry rejects an unexpected bundle version.

```sh
python workflow/scripts/report_refresh.py --manifest /path/to/inventory.json
```

The default is a dry run. It renders into a temporary sibling directory, checks
all report links and sample counts, verifies that scientific records remain
unchanged, and reports the changed filenames. Original report files stay intact.
A sibling lock file prevents simultaneous refreshes by this command.

## Apply and rollback

Run as the report service account, or with permission to preserve its ownership.
Use an idle publication window: the refresh lock does not stop workers, retention
cleanup or another deployment. Each report is a separate transaction; earlier
successful reports remain updated if a later inventory entry fails.

```sh
python workflow/scripts/report_refresh.py --manifest /path/to/inventory.json --apply
python workflow/scripts/report_refresh.py --rollback /path/to/retained-report /path/to/returned-backup
```

After validation, the old directory becomes a uniquely named sibling backup and
the staged directory takes its place. Publication failure restores the original.
The rename pair can briefly make the report unavailable; schedule accordingly.
The command detects changes made to the original while preparing the update.
It preserves file ownership and rebuilds an existing `report.zip` with the same
`report/` layout, excluding the ZIP itself. Rollback restores the entire previous
bundle and ZIP, and retains the displaced version as another backup. Backups are
private operational files: retain them outside published routes and delete them
according to the original report's retention policy.

## Retained evidence and compatibility

New profiles use schema-versioned `*.profile.json` files containing means and SEM
bounds computed directly from the deepTools matrix. SEM follows deepTools:
population standard deviation across non-missing values divided by the square
root of the total matrix row count. Missing bins remain missing. The reader checks
axes and means against the retained TSV and rejects malformed or inconsistent
sidecars. It never substitutes fabricated uncertainty.

New genomic-region figures retain `*.regions.json`. Historical bundles can still
use their profile TSV plus validated SVG bands and embedded region metrics. If
required numeric inputs or identifiable scope sections are absent, refresh fails
before publication. It does not reanalyse missing data, reinterpret old metrics,
or silently adopt a different genome/background.
