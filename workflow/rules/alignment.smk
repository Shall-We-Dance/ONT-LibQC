rule align:
    input:
        query1=query1,
        query2=query2,
        reference=reference_input("fasta"),
        legs=lambda wc: f"{OUT}/prepared/{wc.sample}.legs.fastq.gz" if assay(wc.sample) == "microc" else [],
        index=index_dependency,
        rna_script="workflow/scripts/run_star_rna.py",
    output:
        bam=temp(f"{OUT}/bam/{{sample}}.bam"),
        bai=temp(f"{OUT}/bam/{{sample}}.bam.bai"),
    log:
        f"{OUT}/logs/alignment/{{sample}}.log"
    threads:
        lambda wc: min(4, int(PARAMS.get("alignment_threads", 4))) if assay(wc.sample) in {"rnaseq", "clip"} else int(PARAMS.get("alignment_threads", 4))
    params:
        python=PYTHON,
        starlong=STARLONG,
        rna_summary=lambda wc: f"{OUT}/qc/rna_alignment/{wc.sample}.json",
        samtools=SAMTOOLS,
        minimap2=MINIMAP2,
        bwa=lambda wc: bwa_tool(wc.sample),
        bwa_options=lambda wc: bwa_batch_options(wc.sample),
        sort_memory=str(ALIGNMENT.get('samtools_sort_memory','768M')),
        bowtie2=Bowtie2,
        assay=lambda wc: assay(wc.sample),
        mapper=lambda wc: sample_mapper(wc.sample),
        index=lambda wc: index_prefix(wc),
        minimap_query=lambda wc: f"{OUT}/prepared/{wc.sample}.legs.fastq.gz" if assay(wc.sample) == "microc" else query1(wc),
    conda:
        "../envs/ont-qc.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {log})
        if [ "{params.mapper}" = "starlong" ]; then
          {params.python:q} {input.rna_script:q} --fastq {input.query1:q} --index {params.index:q} \
            --output {output.bam:q} --summary {params.rna_summary:q} --star {params.starlong:q} \
            --samtools {params.samtools:q} --threads {threads} > {log:q} 2>&1
        else
        if [ "{params.mapper}" = "bowtie2" ]; then
          {params.bowtie2} --end-to-end --very-sensitive --no-unal -p {threads} -x {params.index} \
            -U {input.query1} 2> {log}
        elif [ "{params.mapper}" = "bwa-mem2" ] || [ "{params.mapper}" = "bwa-mem" ]; then
          if [ "{params.assay}" = "microc" ]; then
            {params.bwa} mem -SP {params.bwa_options} -t {threads} {params.index} {input.query1} {input.query2} 2> {log}
          else
            {params.bwa} mem {params.bwa_options} -t {threads} {params.index} {input.query1} 2> {log}
          fi
        else
          {params.minimap2} -ax sr -k 9 -w 5 -N 3 --secondary=yes -t {threads} \
            {params.index} {params.minimap_query} 2> {log}
        fi | {params.samtools} sort -@ {threads} -m {params.sort_memory:q} -o {output.bam} -
        fi
        {params.samtools} index -@ {threads} {output.bam}
        """
