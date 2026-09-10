rule mycoplasma_screen:
    input:
        fastq=f"{OUT}/prepared/{{sample}}.inserts.fastq.gz",
        script="workflow/scripts/ont_qc.py",
    output:
        json=temp(f"{OUT}/qc/mycoplasma/{{sample}}.json"),
    log:
        f"{OUT}/logs/mycoplasma/{{sample}}.log"
    params:
        bowtie2=Bowtie2,
        python=PYTHON,
        enabled=lambda wc: "--enabled" if mycoplasma_config(wc.sample).get("enabled", False) else "",
        index=lambda wc: str(mycoplasma_config(wc.sample).get("combined_bowtie2_index", "")),
        regex=lambda wc: str(mycoplasma_config(wc.sample).get("target_contig_regex", "mycoplasma_target|mycoplasmopsis|mycoplasmoides|metamycoplasma|mesomycoplasma|acholeplasma|spiroplasma|ureaplasma|mollicutes")),
        minimum_screened=lambda wc: int(mycoplasma_config(wc.sample).get("minimum_screened_reads", 500)),
    conda:
        "../envs/ont-qc.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {log})
        {params.python} workflow/scripts/ont_qc.py mycoplasma --fastq {input.fastq} --output {output.json} \
          {params.enabled} --index '{params.index}' --target-regex '{params.regex}' --bowtie2 {params.bowtie2} \
          --minimum-screened-reads {params.minimum_screened} > {log} 2>&1
        """

