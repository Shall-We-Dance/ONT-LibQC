include: "workflow/rules/common.smk"
rule all:
    default_target: True
    input:
        REPORT
include: "workflow/rules/prepare.smk"
include: "workflow/rules/alignment.smk"
include: "workflow/rules/spikein.smk"
include: "workflow/rules/mycoplasma.smk"
include: "workflow/rules/ligation.smk"
include: "workflow/rules/metrics.smk"
include: "workflow/rules/tracks.smk"
include: "workflow/rules/reports.smk"


# Cleanup runs only after a validated bundle from this run has been delivered.
onstart:
    from report_bundle import initialize_work
    initialize_work(OUT, REPORT, config)

onsuccess:
    from report_bundle import cleanup_work
    cleanup_work(OUT, REPORT)
