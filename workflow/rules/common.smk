import os
import json
import sys
from pathlib import Path
sys.path.insert(0, str(Path(workflow.basedir)/"workflow/scripts"))
from qc_definitions import configured_params

configfile: "config.yaml"

SAMPLES = list(config["samples"])
OUT = config.get("output_dir", ".work")
REPORT = config.get("report_dir", "report")
from report_bundle import check_locations
check_locations(OUT, REPORT, config)
PARAMS = configured_params(config)
TOOLS = config.get("tools", {})
ALIGNMENT = config.get("alignment", {})
REFERENCES = config["references"]
DEFAULT_REFERENCE = config.get("default_reference", "dm6")
MYCOPLASMA = config.get("mycoplasma_screen", {})
SPIKE_IN = config.get("spike_in", {})
PYTHON = TOOLS.get("python", "python3")
MINIMAP2 = TOOLS.get("minimap2", "minimap2")
Bowtie2 = TOOLS.get("bowtie2", "bowtie2")
BWA_MEM2 = TOOLS.get("bwa_mem2", "bwa-mem2")
SAMTOOLS = TOOLS.get("samtools", "samtools")
BAMCOVERAGE = TOOLS.get("bamCoverage", "bamCoverage")
BEDGRAPH2BW = TOOLS.get("bedGraphToBigWig", "bedGraphToBigWig")
COMPUTEMATRIX = TOOLS.get("computeMatrix", "computeMatrix")
PLOTPROFILE = TOOLS.get("plotProfile", "plotProfile")
MAPQ = int(PARAMS.get("mapq", 20))

def reference_for(sample):
    key = config["samples"][sample].get("reference", DEFAULT_REFERENCE)
    return REFERENCES[key]

def reference_input(key):
    return lambda wc: reference_for(wc.sample)[key]

def annotation_input(wildcards):
    return reference_for(wildcards.sample)["annotation"]

def assay(sample):
    value = str(config["samples"][sample]["type"]).lower().replace("-", "").replace("_", "").replace(" ", "")
    if value in {"microc", "microclike"}:
        return "microc"
    if value in {"hic", "insituhic", "dilutionhic", "3d", "contact"}:
        raise ValueError(f"Sample {sample}: Hi-C must not be processed as Micro-C. See docs/HiC_extension_design.md and specify the restriction enzyme/kit before enabling the Hi-C branch.")
    if value != "chromatin":
        raise ValueError(f"Sample {sample}: unsupported assay type {value!r}")
    return "chromatin"

def fastqs(wildcards):
    files = config["samples"][wildcards.sample]["fastq"]
    if isinstance(files, str):
        files = [files]
    return [x if os.path.isabs(x) else os.path.join(workflow.basedir, x) for x in files]

def sample_mapper(sample):
    sample_override = config["samples"][sample].get("mapper")
    default = ALIGNMENT.get("microc_mapper", "bwa-mem2") if assay(sample) == "microc" else ALIGNMENT.get("chromatin_mapper", "bowtie2")
    mapper = str(sample_override or default).lower()
    if mapper not in {"bowtie2", "bwa-mem2", "minimap2"}:
        raise ValueError(f"Unsupported mapper for {sample}: {mapper}")
    if assay(sample) == "microc" and mapper == "bowtie2":
        raise ValueError("Micro-C requires bwa-mem2 or minimap2; Bowtie2 is only supported for chromatin")
    return mapper

def query1(wildcards):
    suffix = "leg1" if assay(wildcards.sample) == "microc" else "inserts"
    return f"{OUT}/prepared/{wildcards.sample}.{suffix}.fastq.gz"

def query2(wildcards):
    return f"{OUT}/prepared/{wildcards.sample}.leg2.fastq.gz" if assay(wildcards.sample) == "microc" else []

def index_dependency(wildcards):
    return []

def index_prefix(wildcards):
    mapper = sample_mapper(wildcards.sample)
    reference = reference_for(wildcards.sample)
    if mapper == "bowtie2":
        return reference["bowtie2_index"]
    if mapper == "bwa-mem2":
        return reference["bwa_mem2_index"]
    return reference["fasta"]

def expected_index_pairs(wildcards):
    pairs = config["samples"][wildcards.sample].get("index_pairs", [])
    if isinstance(pairs, dict):
        pairs = [{"name": name, **value} for name, value in pairs.items()]
    return ";".join(f"{x.get('name', 'library')}|{x.get('i5', '')}|{x.get('i7', '')}" for x in pairs)

def spikein_config(sample):
    settings = dict(SPIKE_IN)
    override = config["samples"][sample].get("spike_in")
    if isinstance(override, bool):
        settings["enabled"] = override
    elif isinstance(override, dict):
        settings.update(override)
    return settings

def mycoplasma_config(sample):
    settings = dict(MYCOPLASMA)
    settings["combined_bowtie2_index"] = reference_for(sample).get("mycoplasma_index", "")
    if not settings["combined_bowtie2_index"]:
        settings["enabled"] = False
    override = config["samples"][sample].get("mycoplasma_screen")
    if isinstance(override, bool):
        settings["enabled"] = override
    elif isinstance(override, dict):
        settings.update(override)
    return settings

def chromosome_setting(sample, key):
    reference = reference_for(sample)
    if key not in reference or reference[key] is None:
        return "auto"
    values = reference[key]
    if isinstance(values, str):
        return values
    return ",".join(values) if values else "none"

CHROMATIN = [s for s in SAMPLES if assay(s) == "chromatin"]

REFERENCE_NAMES = list(REFERENCES)

def reference_name(sample):
    return config["samples"][sample].get("reference", DEFAULT_REFERENCE)

def annotation_for_reference(wildcards):
    return REFERENCES[wildcards.reference]["annotation"]

def primary_for_reference(wildcards):
    values = REFERENCES[wildcards.reference].get("primary_chromosomes", [])
    return ",".join(values) if not isinstance(values, str) else values


MICROC = [s for s in SAMPLES if assay(s) == "microc"]
