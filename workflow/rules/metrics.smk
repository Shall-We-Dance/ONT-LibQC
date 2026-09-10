rule sample_metrics:
    input:
        ligation=lambda wc: f"{OUT}/qc/ligation/{wc.sample}.json" if assay(wc.sample) == "microc" else [],
        ligation_script="workflow/scripts/microc_ligation.py",
        spikein_script="workflow/scripts/spikein_auto.py",
        bam=f"{OUT}/bam/{{sample}}.bam",
        bai=f"{OUT}/bam/{{sample}}.bam.bai",
        adapter=f"{OUT}/qc/adapters/{{sample}}.summary.json",
        reads=f"{OUT}/qc/adapters/{{sample}}.reads.tsv.gz",
        index_assignments=f"{OUT}/qc/adapters/{{sample}}.index_assignments.tsv.gz",
        mycoplasma=f"{OUT}/qc/mycoplasma/{{sample}}.json",
        spikein=f"{OUT}/qc/spike_in/{{sample}}.json",
        raw=f"{OUT}/qc/lengths/{{sample}}.raw.tsv",
        insert_lengths=f"{OUT}/qc/lengths/{{sample}}.insert.tsv",
        chrom_sizes=reference_input("chrom_sizes"),
        script="workflow/scripts/ont_qc.py",
        plot_metrics="workflow/scripts/qc_plot_metrics.py",
        definitions="workflow/scripts/metric_definitions.json",
        definition_script="workflow/scripts/qc_definitions.py",
    output:
        metrics=temp(f"{OUT}/qc/metrics/{{sample}}.json"),
        fragments=temp(f"{OUT}/fragments/{{sample}}.bed.gz"),
        pairs=temp(f"{OUT}/pairs/{{sample}}.pairs.gz"),
    log:
        f"{OUT}/logs/metrics/{{sample}}.log"
    params:
        ligation=lambda wc: f"{OUT}/qc/ligation/{wc.sample}.json" if assay(wc.sample) == "microc" else "",
        samtools=SAMTOOLS,
        python=PYTHON,
        assay=lambda wc: assay(wc.sample),
        primary=lambda wc: chromosome_setting(wc.sample, "primary_chromosomes"),
        mito=lambda wc: chromosome_setting(wc.sample, "mitochondrial_chromosomes"),
        assembly=lambda wc: str(reference_for(wc.sample).get("assembly", "unspecified")),
        mapper=lambda wc: sample_mapper(wc.sample),
        microc_long_range=int(PARAMS.get("microc_long_range_cis_min_bp", 10000)),
        mapq=MAPQ,
        min_contacts=int(PARAMS.get("microc_minimum_contacts", 200)),
        short_inward=int(PARAMS.get("microc_short_inward_bp", 200)),
        query_coverage=float(PARAMS.get("microc_minimum_query_coverage", .8)),
        dimer_pass=float((config.get("qc_reference_ranges", {}).get("adapter_dimer_rate") or {}).get("max", .05)),
        dimer_warning=float((config.get("qc_reference_ranges", {}).get("adapter_dimer_rate") or {}).get("warning_max", .10)),
        minimum_raw=int(PARAMS.get("minimum_raw_reads", 500)),
        minimum_units=int(PARAMS.get("minimum_complete_units", 200)),
        collision_bin=int(PARAMS.get("collision_qc_bin_size", 10000)),
    conda:
        "../envs/ont-qc.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {log})
        {params.python} workflow/scripts/ont_qc.py metrics \
          --sample {wildcards.sample} --assay {params.assay} --bam {input.bam} --ligation '{params.ligation}'  \
          --samtools {params.samtools} --adapter-summary {input.adapter} \
          --read-table {input.reads} \
          --index-assignments {input.index_assignments} --assembly '{params.assembly}' --mycoplasma {input.mycoplasma} \
          --spike-in {input.spikein} --chrom-sizes {input.chrom_sizes} --primary '{params.primary}' --mito '{params.mito}' \
          --aligner {params.mapper} \
          --mapq {params.mapq} --collision-bin-size {params.collision_bin} \
          --microc-long-range-cis-min-bp {params.microc_long_range} \
          --microc-minimum-contacts {params.min_contacts} --microc-short-inward-bp {params.short_inward} \
          --microc-minimum-query-coverage {params.query_coverage} \
          --dimer-pass-max {params.dimer_pass} --dimer-warning-max {params.dimer_warning} \
          --minimum-raw {params.minimum_raw} --minimum-units {params.minimum_units} \
          --metrics {output.metrics} \
          --fragments {output.fragments} --pairs {output.pairs} \
          > {log} 2>&1
        """

