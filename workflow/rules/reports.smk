rule resolved_report_config:
    output:
        temp(f"{OUT}/qc/resolved_config.json")
    params:
        serialized=json.dumps(dict(config, params=PARAMS), sort_keys=True, indent=2),
    run:
        from pathlib import Path
        Path(str(output[0])).parent.mkdir(parents=True, exist_ok=True)
        Path(str(output[0])).write_text(params.serialized + "\n")

rule reports:
    input:
        script="workflow/scripts/ont_qc.py",
        extension_script="workflow/scripts/spikein_auto.py",
        definitions="workflow/scripts/metric_definitions.json",
        definition_script="workflow/scripts/qc_definitions.py",
        custom_script="workflow/scripts/custom_regions.py",
        custom_ui="workflow/scripts/custom_regions_ui.py",
        custom_browser="workflow/scripts/custom_regions.js",
        modern_css="workflow/scripts/report_modern.css",
        branding="workflow/scripts/report_brand.py",
        brand_mark="workflow/scripts/report_mark.svg",
        reference_guide="workflow/scripts/qc_reference_guide.py",
        mnase=expand(f"{OUT}/qc/ligation/{{sample}}.mnase.json", sample=MICROC),
        custom_beds=sorted({sc["custom_regions"]["bed"] if isinstance(sc["custom_regions"], dict) else sc["custom_regions"] for sc in config["samples"].values() if sc.get("custom_regions")}),
        ligation=expand(f"{OUT}/qc/ligation/{{sample}}.json", sample=MICROC),
        ligation_plots=expand(f"{OUT}/plots/ligation/{{sample}}.total.svg", sample=MICROC),
        ligation_examples=expand(f"{OUT}/plots/ligation/{{sample}}.examples.svg", sample=MICROC),
        mnase_plots=expand(f"{OUT}/plots/ligation/{{sample}}.total.mnase.svg", sample=MICROC),
        plot_metrics="workflow/scripts/qc_plot_metrics.py",
        config=f"{OUT}/qc/resolved_config.json",
        karyotype_script="workflow/scripts/karyotype_plot.R",
        metrics=expand(f"{OUT}/qc/metrics/{{sample}}.json", sample=SAMPLES),
        bam=expand(f"{OUT}/bam/{{sample}}.bam", sample=SAMPLES),
        annotations=sorted({reference_for(sample)["annotation"] for sample in SAMPLES}),
        chrom_sizes=sorted({reference_for(sample)["chrom_sizes"] for sample in SAMPLES}),
        reference_gene_beds=expand(f"{OUT}/reference/{{reference}}.genes.bed", reference=sorted({config["samples"][sample].get("reference", DEFAULT_REFERENCE) for sample in SAMPLES})),
        assignments=expand(f"{OUT}/qc/adapters/{{sample}}.index_assignments.tsv.gz", sample=SAMPLES),
        gene_profile_matrix=expand(f"{OUT}/qc/gene_profile/{{sample}}.matrix.gz", sample=CHROMATIN),
        gene_profile_data=expand(f"{OUT}/qc/gene_profile/{{sample}}.profile.tsv", sample=CHROMATIN),
        gene_profile_plot=expand(f"{OUT}/plots/gene_profile/{{sample}}.TSS_TES.svg", sample=CHROMATIN),
        reads=expand(f"{OUT}/qc/adapters/{{sample}}.reads.tsv.gz", sample=SAMPLES),
        adapter=expand(f"{OUT}/qc/adapters/{{sample}}.summary.json", sample=SAMPLES),
        raw=expand(f"{OUT}/qc/lengths/{{sample}}.raw.tsv", sample=SAMPLES),
        insert_lengths=expand(f"{OUT}/qc/lengths/{{sample}}.insert.tsv", sample=SAMPLES),
    output:
        index=temp(f"{OUT}/summary/index.html"),
        tsv=temp(f"{OUT}/summary/qc_metrics.tsv"),
        sublibrary_tsv=temp(f"{OUT}/summary/sublibrary_qc_metrics.tsv"),
        ligation_tsv=temp(f"{OUT}/summary/ligation_qc_metrics.tsv"),
        mnase_tsv=temp(f"{OUT}/summary/mnase_qc_metrics.tsv"),
        summary_plot=temp(f"{OUT}/plots/summary/all_samples_qc.svg"),
        brand_mark=temp(f"{OUT}/reports/assets/ont-libqc-mark.svg"),
        reports=temp(expand(f"{OUT}/reports/{{sample}}.html", sample=SAMPLES)),
        full_library_weighted=temp(expand(f"{OUT}/reports/assets/{{sample}}.total.full_library_length_weighted.svg", sample=SAMPLES)),
        full_library_count=temp(expand(f"{OUT}/reports/assets/{{sample}}.total.full_library_count_based.svg", sample=SAMPLES)),
        full_library_profiles=temp(expand(f"{OUT}/reports/assets/{{sample}}.total.full_library_profiles.tsv", sample=SAMPLES)),
        full_library_units=temp(expand(f"{OUT}/reports/assets/{{sample}}.total.full_library_units.tsv", sample=SAMPLES)),
        full_library_models=temp(expand(f"{OUT}/reports/assets/{{sample}}.total.full_library_model.json", sample=SAMPLES)),
        tss_assets=temp(expand(f"{OUT}/reports/assets/{{sample}}.TSS_TES.svg", sample=CHROMATIN)),
        karyotype_assets=temp(expand(f"{OUT}/reports/assets/{{sample}}.total.karyotype.svg", sample=SAMPLES)),
    log:
        f"{OUT}/logs/reports/report.log"
    params:
        python=PYTHON,
        out=lambda wc, output: os.path.dirname(os.path.dirname(str(output.index))),
    conda:
        "../envs/ont-qc.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {log})
        {params.python} workflow/scripts/ont_qc.py report --config {input.config} --out {params.out} > {log} 2>&1
        """


rule report_bundle:
    input:
        rendered=rules.reports.output,
        config=rules.resolved_report_config.output,
        script="workflow/scripts/report_bundle.py",
        branding="workflow/scripts/report_brand.py",
        brand_mark="workflow/scripts/report_mark.svg",
        modern_css="workflow/scripts/report_modern.css",
        metrics=expand(f"{OUT}/qc/metrics/{{sample}}.json", sample=SAMPLES),
        reads=expand(f"{OUT}/qc/adapters/{{sample}}.reads.tsv.gz", sample=SAMPLES),
        assignments=expand(f"{OUT}/qc/adapters/{{sample}}.index_assignments.tsv.gz", sample=SAMPLES),
        adapter=expand(f"{OUT}/qc/adapters/{{sample}}.summary.json", sample=SAMPLES),
        mycoplasma=expand(f"{OUT}/qc/mycoplasma/{{sample}}.json", sample=SAMPLES),
        spikein=expand(f"{OUT}/qc/spike_in/{{sample}}.json", sample=SAMPLES),
        spikein_assignments=expand(f"{OUT}/qc/spike_in/{{sample}}.assignments.tsv.gz", sample=SAMPLES),
        gene_profiles=expand(f"{OUT}/qc/gene_profile/{{sample}}.profile.tsv", sample=CHROMATIN),
        ligation=expand(f"{OUT}/qc/ligation/{{sample}}.json", sample=MICROC),
        units=expand(f"{OUT}/qc/ligation/{{sample}}.units.tsv.gz", sample=MICROC),
        segments=expand(f"{OUT}/qc/ligation/{{sample}}.segments.tsv.gz", sample=MICROC),
        mnase=expand(f"{OUT}/qc/ligation/{{sample}}.mnase.json", sample=MICROC),
        mnase_lengths=expand(f"{OUT}/qc/ligation/{{sample}}.mnase_fragments.tsv.gz", sample=MICROC),
    output:
        report=directory(REPORT),
    params:
        python=PYTHON,
        work=OUT,
    conda:
        "../envs/ont-qc.yaml"
    shell:
        "{params.python:q} {input.script:q} --work {params.work:q} --report {output.report:q} --config {input.config:q}"
