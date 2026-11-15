# Maintainability review

This review covers the pipeline CLI, report generation, portal entry points,
workflow dependencies and regression setup. It is a maintenance assessment, not
an exhaustive security audit or a revalidation of the biological methods.

## Changes in this pass

- Move the size-profile CSS and JavaScript out of Python string literals into
  `workflow/scripts/report_profiles.css` and `report_profiles.js`. The renderer
  still embeds both resources, so offline reports and the report CSP continue
  to work. Resource bytes and every renderer function are unchanged.
- Declare the profile renderer, its resources, profile-band reader,
  mitochondrial renderer and spike-in matrix renderer as inputs to the reports
  rule. Previously changes to these helpers were not tracked as report inputs.
  They now participate in Snakemake's normal input-change detection; timestamp
  changes still depend on the configured rerun triggers.
- Syntax-check every browser JavaScript file in CI, including the profile
  interaction and custom-region scripts, instead of maintaining a partial list.

No thresholds, mapping parameters, metric denominators or numerical algorithms
are modified by this pass. Existing published reports do not need regeneration
for this resource extraction.

## Remaining priorities

| Priority | Evidence | Suggested next change |
| --- | --- | --- |
| High | `ont_qc.py` is 2,832 lines. `cmd_report` alone is 274 lines and mixes loading BAM-derived evidence, writing plot assets, rendering HTML and recording provenance. | Extract pure metric formatting and report assembly first, with before/after output comparisons. Keep analysis functions and their CLI unchanged during that step. |
| High | `chromatin_alignment.py` and `rna_motifs.py` import helpers or adapter presets from `ont_qc.py`, which also loads report styles at import time. P5 is defined in both `ont_qc.py` and `library_structure.py`. | Introduce one dependency-free library-sequence module and retain compatibility exports from the CLI. Check adapter orientation and foldback fixtures before migrating consumers. |
| High | `report_profile_band.py` recovers SEM bounds from deepTools SVG geometry; `report_profiles.read_region_metrics` reads metrics embedded in SVG. | Save versioned numeric profile/band and region data alongside new plots. Use those data for new reports; retain the existing validated SVG reader only for historical reports. Do not invent SEM values from means. |
| Medium | `report_modern.css` overrides base CSS embedded in `ont_qc.py`; portal and report selectors recur across later overrides and media queries. | Consolidate by component and breakpoint with desktop, mobile and print comparisons. Repeated selectors inside different media queries are not automatically redundant. |
| Medium | `portal/ont_portal/app.py:create_app` contains 415 lines; `portal/static/app.js` contains 798 lines covering uploads, private links, jobs and status UI. | Separate route groups and browser concerns incrementally. Preserve authentication, retention, legacy links and rate-limit tests at each boundary. File length alone is not a defect. |
| Medium | The chromatin similarity prototype is inactive, but `cmd_report` still calls `rna_report.render_cohort_qc` whenever RNA summaries exist; that renderer includes a cohort-similarity panel and runs `cohort_analysis`. | Reconcile the intended withdrawal across both report paths. Do not assume that retiring the sublibrary panel also retires RNA cohort comparisons. |
| Medium | There is a bundle builder and validator, but no single checked-in command for presentation-only updates of retained private reports and their ZIPs. | Add a manifest-driven refresh command with dry run, validation, backups and rollback, rather than relying on per-change HTML substitutions. Keep deployment inventory and credentials outside the repository. |

## Code that should not be deleted blindly

- `report_similarity.py` and its tests are an explicitly retained, inactive
  prototype (see `SUBLIBRARY_SIMILARITY.md`), not a live pipeline cost.
- Legacy private-link handling preserves access to existing results.
- Historical provenance and SVG compatibility readers preserve old reports.
- Broad exception handlers around worker teardown and transaction rollback
  protect cleanup; assess each handler rather than deleting them by pattern.
- Custom-region computation exists in Python and browser JavaScript because
  uploaded BEDs can be explored offline. Keep parity tests; replacing both with
  a new runtime would be a separate architectural decision.

## Validation and limits

The local maintenance checks passed: 192 pipeline/report tests (9 skipped),
154 portal tests, the Chromium profile interaction check (desktop, keyboard,
touch, offline and print), JavaScript syntax checks, and Snakemake rule parsing.
The initial portal run could not launch an absent Playwright browser build;
using the installed Chromium executable resolved that environment issue.

The skipped tests require optional tools or fixtures. This pass did not run a
full reference-genome mapping workflow, benchmark production throughput or
modify production services. CI also checks Python 3.11/3.13 and Chromium/WebKit.
