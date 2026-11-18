def species_indexes(sample):
    ref = config["samples"][sample].get("reference", DEFAULT_REFERENCE)
    panel = dict(config["species_composition"]["indexes"])
    panel[ref] = reference_for(sample)["bowtie2_index"]
    return panel


rule species_composition:
    input:
        raw=fastqs,
        table=rules.prepare_library_units.output.reads,
        summary=rules.prepare_library_units.output.summary,
        script="workflow/scripts/species_composition.py",
        sequences="workflow/scripts/library_sequences.py",
    output:
        temp(f"{OUT}/qc/species_composition/{{sample}}.json"),
    threads: 4
    params:
        python=PYTHON,
        bowtie2=Bowtie2,
        indexes=lambda wc: json.dumps(species_indexes(wc.sample)),
        reference=lambda wc: config["samples"][wc.sample].get("reference", DEFAULT_REFERENCE),
    log:
        f"{OUT}/logs/species_composition/{{sample}}.log"
    conda:
        "../envs/ont-qc.yaml"
    shell:
        "{params.python:q} {input.script:q} --raw {input.raw:q} --table {input.table:q} "
        "--summary {input.summary:q} --reference {params.reference:q} --indexes {params.indexes:q} "
        "--bowtie2 {params.bowtie2:q} --threads {threads} --output {output:q} > {log:q} 2>&1"
