rule microc_full_insert_ligation:
    input:
        fastq=f"{OUT}/prepared/{{sample}}.inserts.fastq.gz",
        assignments=f"{OUT}/qc/adapters/{{sample}}.index_assignments.tsv.gz",
        script="workflow/scripts/microc_ligation.py",
        mnase_script="workflow/scripts/microc_mnase.py",
        trace_script="workflow/scripts/qc_plot_metrics.py",
        reference=reference_input("fasta"),
    output:
        sam=temp(f"{OUT}/qc/ligation/{{sample}}.full_insert.sam"),
        json=temp(f"{OUT}/qc/ligation/{{sample}}.json"),
        units=temp(f"{OUT}/qc/ligation/{{sample}}.units.tsv.gz"),
        segments=temp(f"{OUT}/qc/ligation/{{sample}}.segments.tsv.gz"),
        plot=temp(f"{OUT}/plots/ligation/{{sample}}.total.svg"),
        examples=temp(f"{OUT}/plots/ligation/{{sample}}.examples.svg"),
        mnase=temp(f"{OUT}/qc/ligation/{{sample}}.mnase.json"),
        mnase_lengths=temp(f"{OUT}/qc/ligation/{{sample}}.mnase_fragments.tsv.gz"),
        mnase_plot=temp(f"{OUT}/plots/ligation/{{sample}}.total.mnase.svg"),
    log:
        f"{OUT}/logs/ligation/{{sample}}.log"
    threads: 4
    wildcard_constraints:
        sample="|".join(MICROC)
    params:
        bwa=lambda wc: bwa_tool(wc.sample),
        aligner=lambda wc: microc_bwa_mapper(wc.sample),
        bwa_options=lambda wc: bwa_batch_options(wc.sample),
        python=PYTHON,
        index=lambda wc: bwa_index(wc.sample),
        primary=lambda wc: chromosome_setting(wc.sample,"primary_chromosomes"),
        plot_prefix=lambda wc: f"{OUT}/plots/ligation/{wc.sample}",
        mapq=MAPQ,
        segment=int(PARAMS.get("ligation_minimum_segment_bp",40)),
        coverage=float(PARAMS.get("ligation_minimum_query_coverage",.8)),
        gap=int(PARAMS.get("ligation_maximum_query_gap",30)),
    conda:
        "../envs/ont-qc.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.sam}) $(dirname {log})
        {params.bwa} mem {params.bwa_options} -t {threads} -Y -k 15 -T 20 {params.index} {input.fastq} > {output.sam} 2> {log}
        {params.python} {input.script} --aligner {params.aligner} --sample {wildcards.sample} --fastq {input.fastq} --sam {output.sam} \
          --assignments {input.assignments} --primary '{params.primary}' --output {output.json} \
          --minimum-segment-bp {params.segment} --minimum-query-coverage {params.coverage} \
          --maximum-query-gap {params.gap} --mapq {params.mapq} --plot-prefix {params.plot_prefix} >> {log} 2>&1
        """
