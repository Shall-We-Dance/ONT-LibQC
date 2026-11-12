rule rna_prepare:
    input:
        script="workflow/scripts/rna_prepare.py",
        inserts=f"{OUT}/prepared/{{sample}}.inserts.fastq.gz",
        reads=f"{OUT}/qc/adapters/{{sample}}.reads.tsv.gz",
    output:
        fastq=temp(f"{OUT}/prepared/{{sample}}.rna.fastq.gz"),
        summary=temp(f"{OUT}/qc/rna_preparation/{{sample}}.json"),
        units=temp(f"{OUT}/qc/rna_preparation/{{sample}}.units.tsv.gz"),
    wildcard_constraints:
        sample="|".join(RNA) if RNA else "(?!)"
    params:
        python=PYTHON,
        p5=lambda wc: int(config["samples"][wc.sample].get("rna_trim_p5_bp", 0)),
        p7=lambda wc: int(config["samples"][wc.sample].get("rna_trim_p7_bp", 0)),
        structure=library_structure_config,
    conda:
        "../envs/ont-qc.yaml"
    shell:
        "{params.python:q} {input.script:q} --fastq {input.inserts:q} --read-table {input.reads:q} --output {output.fastq:q} --summary {output.summary:q} --units-tsv {output.units:q} --structure-json {params.structure:q} --p5 {params.p5} --p7 {params.p7}"


rule rna_sources:
    input:
        script="workflow/scripts/rna_sources.py",
        inserts=f"{OUT}/prepared/{{sample}}.rna.fastq.gz",
        reference="workflow/resources/rna_sources/panel.fa",
        manifest="workflow/resources/rna_sources/panel.json",
    output:
        units=temp(f"{OUT}/qc/rna_sources/{{sample}}.units.tsv.gz"),
        summary=temp(f"{OUT}/qc/rna_sources/{{sample}}.json"),
    log:
        f"{OUT}/logs/rna_sources/{{sample}}.log"
    wildcard_constraints:
        sample="|".join(RNA) if RNA else "(?!)"
    params:
        python=PYTHON,
        minimap2=MINIMAP2,
    conda:
        "../envs/ont-qc.yaml"
    shell:
        r"""
        mkdir -p $(dirname {log:q})
        {params.python:q} {input.script:q} --input {input.inserts:q} \
          --reference-fasta {input.reference:q} --reference-manifest {input.manifest:q} \
          --output {output.units:q} --summary {output.summary:q} --minimap2 {params.minimap2:q} > {log:q} 2>&1
        """


rule rna_qc:
    input:
        script="workflow/scripts/rna_qc.py",
        motifs="workflow/scripts/rna_motifs.py",
        bam=f"{OUT}/bam/{{sample}}.bam",
        bai=f"{OUT}/bam/{{sample}}.bam.bai",
        annotation=annotation_input,
        reads=f"{OUT}/qc/adapters/{{sample}}.reads.tsv.gz",
        assignments=f"{OUT}/qc/adapters/{{sample}}.index_assignments.tsv.gz",
        inserts=f"{OUT}/prepared/{{sample}}.rna.fastq.gz",
        preparation=f"{OUT}/qc/rna_preparation/{{sample}}.json",
        preparation_units=f"{OUT}/qc/rna_preparation/{{sample}}.units.tsv.gz",
        source_units=f"{OUT}/qc/rna_sources/{{sample}}.units.tsv.gz",
        source_summary=f"{OUT}/qc/rna_sources/{{sample}}.json",
    output:
        evidence=temp(directory(f"{OUT}/qc/rna/{{sample}}")),
    log:
        f"{OUT}/logs/rna/{{sample}}.log"
    wildcard_constraints:
        sample="|".join(RNA) if RNA else "(?!)"
    params:
        python=PYTHON,
        samtools=SAMTOOLS,
        assay=lambda wc: assay(wc.sample),
        mapq=MAPQ,
        assembly=lambda wc: str(reference_for(wc.sample).get("assembly", reference_name(wc.sample))),
    conda:
        "../envs/ont-qc.yaml"
    shell:
        r"""
        mkdir -p $(dirname {log:q})
        {params.python:q} {input.script:q} --bam {input.bam:q} --gtf {input.annotation:q} \
          --read-table {input.reads:q} --fastq {input.inserts:q} \
          --index-assignments {input.assignments:q} --preparation-json {input.preparation:q} --assembly {params.assembly:q} \
          --preparation-tsv {input.preparation_units:q} --source-assignments {input.source_units:q} --source-summary {input.source_summary:q} \
          --sample {wildcards.sample:q} --assay {params.assay:q} --out-dir {output.evidence:q} \
          --samtools {params.samtools:q} --mapq {params.mapq} > {log:q} 2>&1
        """
