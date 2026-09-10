rule prepare_library_units:
    input:
        fastqs=fastqs,
        script="workflow/scripts/ont_qc.py",
    output:
        inserts=temp(f"{OUT}/prepared/{{sample}}.inserts.fastq.gz"),
        legs=temp(f"{OUT}/prepared/{{sample}}.legs.fastq.gz"),
        leg1=temp(f"{OUT}/prepared/{{sample}}.leg1.fastq.gz"),
        leg2=temp(f"{OUT}/prepared/{{sample}}.leg2.fastq.gz"),
        reads=temp(f"{OUT}/qc/adapters/{{sample}}.reads.tsv.gz"),
        index_assignments=temp(f"{OUT}/qc/adapters/{{sample}}.index_assignments.tsv.gz"),
        summary=temp(f"{OUT}/qc/adapters/{{sample}}.summary.json"),
        raw_lengths=temp(f"{OUT}/qc/lengths/{{sample}}.raw.tsv"),
        insert_lengths=temp(f"{OUT}/qc/lengths/{{sample}}.insert.tsv"),
    log:
        f"{OUT}/logs/prepare/{{sample}}.log"
    params:
        python=PYTHON,
        assay=lambda wc: assay(wc.sample),
        leg=int(PARAMS.get("microc_leg_length", 100)),
        max_unit=int(PARAMS.get("maximum_library_unit_bp", 5000)),
        adapter_dimer_max=int(PARAMS.get("adapter_dimer_max_bp", 10)),
        ultrashort_max=int(PARAMS.get("ultrashort_insert_max_bp", 50)),
        dimer_raw_min=int(PARAMS.get("adapter_dimer_raw_min_bp", 120)),
        dimer_raw_max=int(PARAMS.get("adapter_dimer_raw_max_bp", 170)),
        dimer_identity=float(PARAMS.get("adapter_dimer_template_identity", 0.80)),
        single_pair_fraction=float(PARAMS.get("single_index_pair_fraction", 0.90)),
        min_index_info=int(PARAMS.get("minimum_index_informative_units", 20)),
        min_index_units=int(PARAMS.get("minimum_index_cluster_units", 20)),
        min_index_fraction=float(PARAMS.get("minimum_index_cluster_fraction", 0.01)),
        index_variant_max_edit=int(PARAMS.get("index_variant_max_edit", 2)),
        index_variant_dominance_ratio=float(PARAMS.get("index_variant_dominance_ratio", 5.0)),
        expected_pairs=expected_index_pairs,
    conda:
        "../envs/ont-qc.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {log})
        {params.python} workflow/scripts/ont_qc.py prepare \
          --sample {wildcards.sample} --assay {params.assay} \
          --fastq {input.fastqs} --inserts {output.inserts} --legs {output.legs} \
          --leg1 {output.leg1} --leg2 {output.leg2} \
          --read-table {output.reads} --index-assignments {output.index_assignments} --summary {output.summary} \
          --raw-lengths {output.raw_lengths} --insert-lengths {output.insert_lengths} \
          --leg-length {params.leg} --max-unit {params.max_unit} \
          --adapter-dimer-max-bp {params.adapter_dimer_max} \
          --ultrashort-insert-max-bp {params.ultrashort_max} \
          --adapter-dimer-raw-min-bp {params.dimer_raw_min} --adapter-dimer-raw-max-bp {params.dimer_raw_max} \
          --adapter-dimer-template-identity {params.dimer_identity} \
          --single-index-pair-fraction {params.single_pair_fraction} \
          --minimum-index-informative-units {params.min_index_info} \
          --minimum-index-cluster-units {params.min_index_units} \
          --minimum-index-cluster-fraction {params.min_index_fraction} \
          --index-variant-max-edit {params.index_variant_max_edit} \
          --index-variant-dominance-ratio {params.index_variant_dominance_ratio} \
          --expected-index-pairs '{params.expected_pairs}' \
          > {log} 2>&1
        """

