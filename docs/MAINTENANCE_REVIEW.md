# Maintainability review and remediation

All seven issues from the maintenance review now have bounded implementations.
Scientific thresholds, mapping parameters, backgrounds and metric denominators
are unchanged. This is not an exhaustive security audit or biological validation.

| Reviewed issue | Resolution |
| --- | --- |
| CLI and report monolith | `ont_qc.py` is a compatibility facade and CLI. `qc_analysis.py` owns analysis helpers; `report_metrics.py` owns formatting; `report_render.py` separates per-sample rendering, summary/provenance writing and command orchestration. Existing helper imports remain compatible. |
| Analysis imports report entry point; duplicate adapters | `library_sequences.py` owns adapter sequences, reverse complementation, FASTQ iteration and edit distance. Alignment, motif and structure consumers import it directly without loading report presentation. |
| Plot data recovered from SVG | New reports retain versioned numerical profile/SEM and genomic-region JSON alongside TSV/SVG assets. Bundling preserves them through intermediate cleanup. Numeric readers validate compatibility; SVG readers remain only for historical reports. |
| Accumulated CSS overrides | Report base styles and tooltip JS are external resources. Stylesheets use readable declarations and preserve cascade order/media conditions. 100 superseded declarations and inactive RNA-comparison portal styles were removed. |
| Portal backend/frontend monoliths | The application factory registers job, link, report and public route groups with shared access services. Browser task cards and public server statistics have explicit controller interfaces in separate files. |
| Incomplete similarity withdrawal | Formal report generation no longer calls RNA cohort comparison. Its renderer is isolated as an inactive prototype; the chromatin prototype and scientific tests remain available for future work. |
| Ad hoc report refresh scripts | `report_refresh.py` provides an external inventory, dry run, validation, immutable scientific records, ownership preservation, ZIP rebuild, per-report backups, failure recovery and rollback. See `REPORT_REFRESH.md`. |

Workflow inputs track the extracted analysis/report resources. CI checks every
browser JavaScript file. The release remains v1.1.0; this maintenance change does
not claim a new release or alter historical analysis timestamps.

## Validation

Local checks passed: 205 pipeline/report tests (9 optional-tool skips) and 154
portal tests, plus
numeric SEM agreement with deepTools' masked-array formula, sidecar retention
through bundle cleanup, transactional refresh/rollback, browser uploads, public
counts, short/legacy private links, profile highlighting and offline/print paths.
Report and portal computed styles were compared before/after consolidation at
1440, 900, 390 and 320 px and in print media; all computed properties matched on
the tested pages. The existing public example bundle also passed refresh dry run.

Full genome mapping and production deployment are separate from these maintenance
checks. Local browser validation uses Chromium; CI additionally runs WebKit and
Python 3.11/3.13. Optional-tool skips remain explicit in the test results.

## Compatibility deliberately retained

- Long private links, report grants and historical provenance remain supported.
- Existing validated SVG readers remain for bundles without numeric sidecars.
- Inactive similarity mathematics and tests are retained, outside live reports.
- Worker cleanup/transaction exception handling remains in place.
- Custom BED calculation remains in Python and browser JavaScript so downloaded
  reports continue to work offline; parity tests cover both implementations.

Further changes should follow these boundaries, rather than adding report CSS,
route handlers or sequence constants back into the CLI/application factory.
