"""Streaming reference normalization and bounded preparation resources."""
import math
import re
import shutil
from pathlib import Path

GIB = 1024**3


def memory_bytes(value):
    m = re.fullmatch(r"([1-9][0-9]*)([kmgt]?)", str(value).lower())
    if not m:
        raise ValueError("Memory sizes must be bytes or an integer followed by k, m, g or t")
    return int(m[1]) * 1024**("kmgt".index(m[2])+1 if m[2] else 0)


def tool_memory(program, bases, maximum):
    # BWA-MEM2 documents indexing memory of approximately 28 bytes per base.
    factor = 28 if program == "bwa-mem2" else 6 if program == "bowtie2-build" else 16 if program == "minimap2" else 0
    needed = max(2*GIB, math.ceil((factor*bases+2*GIB)/GIB)*GIB)
    cap = memory_bytes(maximum)
    if needed > cap:
        raise ValueError(f"{program} needs an estimated {needed//GIB} GiB; --memory caps builds at {cap/GIB:g} GiB. Increase the build cap on a suitable host.")
    return needed


def check_headroom(folder, memory, reserve):
    meminfo = Path("/proc/meminfo")
    if meminfo.is_file():
        available = next((int(s.split()[1])*1024 for s in meminfo.read_text().splitlines() if s.startswith("MemAvailable:")), 0)
        if available < memory+reserve:
            raise ValueError("Insufficient available RAM for this reference build and the reserved headroom. Rerun when resources are available.")
    if shutil.disk_usage(folder).free < 2*GIB:
        raise ValueError("Reference build needs at least 2 GiB of remaining disk space")


def assembly_aliases(report):
    aliases = {}
    for line in Path(report).read_text().splitlines():
        if not line.strip() or line.startswith("#"):
            continue
        fields = line.split("\t")
        if len(fields) < 9:
            raise ValueError("Malformed NCBI assembly report")
        name, role, molecule, location, genbank, relation, refseq, unit, length = fields[:9]
        ucsc = fields[9] if len(fields)>9 else "na"
        target = ucsc if ucsc != "na" else (("chrM" if molecule=="MT" else "chr"+molecule) if role=="assembled-molecule" else refseq)
        if target=="na":
            target=genbank
        if target=="na" or not re.fullmatch(r"[A-Za-z0-9_.-]+", target):
            raise ValueError("Unsupported chromosome identifier in assembly report")
        for alias in (name, genbank, refseq):
            if alias != "na":
                if alias in aliases and aliases[alias] != target:
                    raise ValueError("Conflicting chromosome aliases")
                aliases[alias]=target
        aliases[target]=target
    if not aliases:
        raise ValueError("Empty chromosome alias table")
    return aliases


def normalize_names(stage):
    """Rename both sources consistently without loading chromosomes into memory."""
    stage=Path(stage);aliases=assembly_aliases(stage/"assembly_report.txt")
    seen=set();features=0
    with (stage/"genome.fa").open() as source, (stage/"genome.renamed.fa").open("w") as out:
        for line in source:
            if line.startswith(">"):
                old=line[1:].split()[0]
                if old not in aliases:
                    raise ValueError("FASTA chromosome absent from assembly report: "+old)
                name=aliases[old]
                if name in seen:
                    raise ValueError("Duplicate normalized FASTA chromosome: "+name)
                seen.add(name);out.write(">"+name+"\n")
            else:
                if line.strip() and not seen:
                    raise ValueError("Malformed reference FASTA")
                out.write(line)
    with (stage/"genes.gtf").open() as source, (stage/"genes.renamed.gtf").open("w") as out:
        for line in source:
            if not line.strip() or line.startswith("#"):
                out.write(line);continue
            fields=line.rstrip("\n").split("\t")
            if len(fields)!=9 or fields[0] not in aliases or aliases[fields[0]] not in seen:
                raise ValueError("Annotation chromosome is absent from the matching FASTA")
            fields[0]=aliases[fields[0]]
            features+=fields[2] in {"gene","transcript","exon"}
            out.write("\t".join(fields)+"\n")
    if not seen or not features:
        raise ValueError("Reference needs FASTA sequences and gene annotation features")
    (stage/"genome.renamed.fa").replace(stage/"genome.fa")
    (stage/"genes.renamed.gtf").replace(stage/"genes.gtf")
