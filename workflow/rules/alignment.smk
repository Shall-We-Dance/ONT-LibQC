rule align:
    input:
        query1=query1,
        query2=query2,
        reference=reference_input("fasta"),
        legs=lambda wc: f"{OUT}/prepared/{wc.sample}.legs.fastq.gz" if assay(wc.sample) == "microc" else [],
        index=index_dependency,
    output:
        bam=temp(f"{OUT}/bam/{{sample}}.bam"),
        bai=temp(f"{OUT}/bam/{{sample}}.bam.bai"),
    log:
        f"{OUT}/logs/alignment/{{sample}}.log"
    threads:
        int(PARAMS.get("alignment_threads", 4))
    params:
        samtools=SAMTOOLS,
        minimap2=MINIMAP2,
        bwa_mem2=BWA_MEM2,
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
        if [ "{params.mapper}" = "bowtie2" ]; then
          {params.bowtie2} --end-to-end --very-sensitive --no-unal -p {threads} -x {params.index} \
            -U {input.query1} 2> {log}
        elif [ "{params.mapper}" = "bwa-mem2" ]; then
          if [ "{params.assay}" = "microc" ]; then
            {params.bwa_mem2} mem -SP -t {threads} {params.index} {input.query1} {input.query2} 2> {log}
          else
            {params.bwa_mem2} mem -t {threads} {params.index} {input.query1} 2> {log}
          fi
        else
          {params.minimap2} -ax sr -k 9 -w 5 -N 3 --secondary=yes -t {threads} \
            {params.index} {params.minimap_query} 2> {log}
        fi | {params.samtools} sort -@ {threads} -o {output.bam} -
        {params.samtools} index -@ {threads} {output.bam}
        """

