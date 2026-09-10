rule chromatin_bigwig:
    input:
        bam=f"{OUT}/bam/{{sample}}.bam",
        bai=f"{OUT}/bam/{{sample}}.bam.bai",
        chrom_sizes=reference_input("chrom_sizes"),
        script="workflow/scripts/ont_qc.py",
    output:
        temp(f"{OUT}/tracks/{{sample}}.CPM.bw")
    log:
        f"{OUT}/logs/tracks/{{sample}}.log"
    params:
        python=PYTHON,
        samtools=SAMTOOLS,
        converter=BEDGRAPH2BW,
        mapq=MAPQ,
        primary=lambda wc: chromosome_setting(wc.sample, "primary_chromosomes"),
    wildcard_constraints:
        sample="|".join(CHROMATIN)
    conda:
        "../envs/ont-qc.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {log})
        {params.python} workflow/scripts/ont_qc.py track --bam {input.bam} --output {output} \
          --samtools {params.samtools} --converter {params.converter} --mapq {params.mapq} \
          --chrom-sizes {input.chrom_sizes} --primary '{params.primary}' > {log} 2>&1
        """



rule reference_gene_bed:
    input:
        annotation=annotation_for_reference,
        script="workflow/scripts/gtf_gene_bed.awk",
    output:
        temp(f"{OUT}/reference/{{reference}}.genes.bed")
    params:
        primary=primary_for_reference,
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output})
        awk -v primary='{params.primary}' -f {input.script} {input.annotation} | sort -k1,1 -k2,2n > {output}
        """

rule chromatin_gene_profile:
    input:
        bigwig=f"{OUT}/tracks/{{sample}}.CPM.bw",
        genes=lambda wc: f"{OUT}/reference/{reference_name(wc.sample)}.genes.bed",
    output:
        matrix=temp(f"{OUT}/qc/gene_profile/{{sample}}.matrix.gz"),
        data=temp(f"{OUT}/qc/gene_profile/{{sample}}.profile.tsv"),
        plot=temp(f"{OUT}/plots/gene_profile/{{sample}}.TSS_TES.svg"),
    log:
        f"{OUT}/logs/gene_profile/{{sample}}.log"
    threads:
        int(PARAMS.get("gene_profile_threads", 2))
    params:
        compute_matrix=COMPUTEMATRIX,
        plot_profile=PLOTPROFILE,
        upstream=int(PARAMS.get("gene_profile_upstream_bp", 5000)),
        downstream=int(PARAMS.get("gene_profile_downstream_bp", 5000)),
        body=int(PARAMS.get("gene_profile_body_bp", 5000)),
        bin=int(PARAMS.get("gene_profile_bin_size", 100)),
    wildcard_constraints:
        sample="|".join(CHROMATIN)
    conda:
        "../envs/ont-qc.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.matrix}) $(dirname {output.plot}) $(dirname {log})
        {params.compute_matrix} scale-regions \
          -S {input.bigwig} -R {input.genes} \
          --beforeRegionStartLength {params.upstream} \
          --afterRegionStartLength {params.downstream} \
          --regionBodyLength {params.body} --binSize {params.bin} \
          --missingDataAsZero --sortRegions keep \
          --startLabel TSS --endLabel TES -p {threads} \
          -o {output.matrix} > {log} 2>&1
        {params.plot_profile} -m {output.matrix} -o {output.plot} \
          --outFileNameData {output.data} --plotType se \
          --startLabel TSS --endLabel TES --yAxisLabel "CPM-normalized signal" \
          >> {log} 2>&1
        """
