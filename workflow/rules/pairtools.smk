rule microc_pairtools:
    input:
        script="workflow/scripts/microc_pairtools.py",
        inserts=f"{OUT}/prepared/{{sample}}.inserts.fastq.gz",
        reads=f"{OUT}/qc/adapters/{{sample}}.reads.tsv.gz",
        adapter=f"{OUT}/qc/adapters/{{sample}}.summary.json",
        chrom_sizes=reference_input("chrom_sizes"),
    output:
        json=temp(f"{OUT}/qc/pairtools/{{sample}}.json"),
        summary=temp(f"{OUT}/qc/pairtools/{{sample}}.summary.tsv"),
        units=temp(f"{OUT}/qc/pairtools/{{sample}}.units.tsv.gz"),
    log:
        f"{OUT}/logs/pairtools/{{sample}}.log"
    threads: 4
    params:
        python=PYTHON,
        bwa=BWA_MEM2,
        pairtools=TOOLS.get("pairtools", "pairtools"),
        fastp=TOOLS.get("fastp", "fastp"),
        index=lambda wc: reference_for(wc.sample)["bwa_mem2_index"],
        prefix=lambda wc: f"{OUT}/qc/pairtools/{wc.sample}",
    conda:
        "../envs/ont-qc.yaml"
    shell:
        r"""
        mkdir -p $(dirname {log:q})
        {params.python:q} {input.script:q} --sample {wildcards.sample:q} \
          --inserts {input.inserts:q} --read-table {input.reads:q} --adapter-summary {input.adapter:q} \
          --chrom-sizes {input.chrom_sizes:q} --index {params.index:q} \
          --bwa-mem2 {params.bwa:q} --pairtools {params.pairtools:q} --fastp {params.fastp:q} \
          --threads {threads} --output-prefix {params.prefix:q} > {log:q} 2>&1
        """
