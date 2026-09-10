rule spikein_screen:
    threads: 4
    input:
        fastq=f"{OUT}/prepared/{{sample}}.inserts.fastq.gz",
        script="workflow/scripts/ont_qc.py",
        auto_script="workflow/scripts/spikein_auto.py",
        database=lambda wc: [spikein_config(wc.sample).get("database_manifest", "resources/spikein_v1/manifest.json")] if spikein_config(wc.sample).get("enabled",False) and spikein_config(wc.sample).get("mode")=="auto" else [],
    output:
        json=temp(f"{OUT}/qc/spike_in/{{sample}}.json"),
        assignments=temp(f"{OUT}/qc/spike_in/{{sample}}.assignments.tsv.gz"),
    log:
        f"{OUT}/logs/spike_in/{{sample}}.log"
    params:
        mode=lambda wc: spikein_config(wc.sample).get("mode", "configured") if spikein_config(wc.sample).get("enabled",False) else "disabled",
        database=lambda wc: spikein_config(wc.sample).get("database_manifest", "resources/spikein_v1/manifest.json"),
        host=lambda wc: reference_name(wc.sample),
        declared=lambda wc: "--declared-spike-in" if spikein_config(wc.sample).get("declared_spike_in",False) else "",
        bowtie2=Bowtie2,
        python=PYTHON,
        enabled=lambda wc: "--enabled" if spikein_config(wc.sample).get("enabled", False) else "",
        index=lambda wc: str(spikein_config(wc.sample).get("combined_bowtie2_index", "")),
        regex=lambda wc: str(spikein_config(wc.sample).get("target_contig_regex", "")),
        label=lambda wc: str(spikein_config(wc.sample).get("label", "spike-in")),
        mapq=lambda wc: int(spikein_config(wc.sample).get("mapq", MAPQ)),
        minimum_aligned=lambda wc: int(spikein_config(wc.sample).get("minimum_aligned_bp", 80)),
        minimum_coverage=lambda wc: float(spikein_config(wc.sample).get("minimum_query_coverage", 0.80)),
        minimum_identity=lambda wc: float(spikein_config(wc.sample).get("minimum_identity", 0.85)),
        minimum_screened=lambda wc: int(spikein_config(wc.sample).get("minimum_screened_reads", 500)),
    conda:
        "../envs/ont-qc.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {log})
        if [ "{params.mode}" = "auto" ]; then
          {params.python} {input.auto_script} --fastq {input.fastq} --output {output.json} --database {params.database} \
            --host {params.host} --bowtie2 {params.bowtie2} --threads {threads} {params.declared} --mapq {params.mapq} \
            --minimum-aligned-bp {params.minimum_aligned} --minimum-query-coverage {params.minimum_coverage} \
            --minimum-identity {params.minimum_identity} > {log} 2>&1
        else
        {params.python} workflow/scripts/ont_qc.py spikein --fastq {input.fastq} --output {output.json}           {params.enabled} --index '{params.index}' --target-regex '{params.regex}' --label '{params.label}'           --bowtie2 {params.bowtie2} --mapq {params.mapq} --minimum-aligned-bp {params.minimum_aligned}           --minimum-query-coverage {params.minimum_coverage} --minimum-identity {params.minimum_identity}           --minimum-screened-reads {params.minimum_screened} > {log} 2>&1
        fi
        if [ ! -f {output.assignments:q} ]; then gzip -c < /dev/null > {output.assignments:q}; fi
        """

