#!/usr/bin/env python3
"""Build memory-conscious STARlong indexes from an operator-owned catalog.

This never downloads data and never edits the input catalog. Indexes are built in
staging, verified, then atomically made available. Rerunning reuses only indexes
whose source hashes, version and parameters still match.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import re
import shutil
import subprocess
import time

VERSION = "2.7.11b"
CORE = ("Genome", "SA", "SAindex", "genomeParameters.txt", "chrNameLength.txt", "chrName.txt", "chrLength.txt", "chrStart.txt")
GIB = 1024**3


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024**2), b""):
            digest.update(block)
    return digest.hexdigest()


def read_sizes(path):
    sizes = {}
    for line in Path(path).read_text().splitlines():
        if line.strip():
            name, size = line.split()[:2]
            sizes[name] = int(size)
    if not sizes or any(size < 1 for size in sizes.values()):
        raise ValueError(f"Invalid chromosome sizes: {path}")
    return sizes


def atomic_json(path, value):
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(json.dumps(value, indent=2) + "\n")
    os.replace(tmp, path)


def build_one(reference_id, reference, args):
    if not re.fullmatch(r"[A-Za-z0-9_.-]+", reference_id) or reference_id in {".", ".."}:
        raise ValueError(f"Invalid reference identifier: {reference_id!r}")
    fasta = Path(reference["fasta"]).resolve(strict=True)
    annotation = Path(reference["annotation"]).resolve(strict=True)
    sizes = read_sizes(reference["chrom_sizes"])
    sequence_bases = sum(sizes.values())
    parameters = {
        "genomeSAsparseD": 4,
        "genomeSAindexNbases": min(12, max(1, int(math.log2(sequence_bases) / 2 - 1))),
        "genomeChrBinNbits": min(18, max(1, int(math.log2(max(sequence_bases / len(sizes), 1000))))),
        "sjdbOverhang": 100,
    }
    source_hashes = {"fasta": sha256(fasta), "annotation": sha256(annotation)}
    target = args.output / reference_id
    manifest_path = target / "ontlibqc-star-index.json"
    if target.exists():
        manifest = json.loads(manifest_path.read_text())
        if manifest.get("star_version") != VERSION or manifest.get("parameters") != parameters or manifest.get("source_sha256") != source_hashes:
            raise ValueError(f"Existing index differs from requested source/parameters: {target}")
        for name in CORE:
            size = manifest.get("index_files", {}).get(name)
            if not isinstance(size, int) or size < 1 or not (target / name).is_file() or (target / name).stat().st_size != size:
                raise ValueError(f"Incomplete existing STAR index: {target / name}")
        print(f"{reference_id}: verified existing {target}", flush=True)
    else:
        staging = args.output / f".building-{reference_id}"
        staging.mkdir(parents=True, exist_ok=True)
        command = [args.star, "--runMode", "genomeGenerate", "--runThreadN", str(args.threads),
                   "--genomeDir", str(staging), "--genomeFastaFiles", str(fasta),
                   "--sjdbGTFfile", str(annotation), "--limitGenomeGenerateRAM", str(args.build_memory_gib * GIB)]
        for key, value in parameters.items():
            command.extend(["--" + key, str(value)])
        print(f"{reference_id}: building sparse STAR index", flush=True)
        start = time.monotonic()
        with (staging / "build.log").open("w") as log:
            subprocess.run(command, check=True, stdout=log, stderr=subprocess.STDOUT, cwd=staging)
        for name in CORE:
            if not (staging / name).is_file() or not (staging / name).stat().st_size:
                raise ValueError(f"STAR did not produce a complete index: {name}")
        if read_sizes(staging / "chrNameLength.txt") != sizes:
            raise ValueError(f"STAR chromosome names or lengths differ from catalog: {reference_id}")
        index_files = {name: (staging / name).stat().st_size for name in CORE}
        loaded_bytes = sum(index_files[name] for name in ("Genome", "SA", "SAindex"))
        manifest = {
            "schema_version": 1, "reference_id": reference_id, "star_version": VERSION,
            "genome_format_version": "2.7.4a",
            "parameters": parameters, "source_sha256": source_hashes,
            "index_files": index_files,
            "star_alignment_memory_bytes": math.ceil((loaded_bytes + 3 * GIB) / GIB) * GIB,
            "build_seconds": round(time.monotonic() - start, 3),
            "build_threads": args.threads, "build_memory_limit_bytes": args.build_memory_gib * GIB,
            "sequence_bases": sequence_bases,
            "annotation_note": "Annotation is embedded; no runtime junction insertion is needed.",
        }
        atomic_json(staging / "ontlibqc-star-index.json", manifest)
        staging.rename(target)
        print(f"{reference_id}: complete; mapping reservation {manifest['star_alignment_memory_bytes'] / GIB:g} GiB", flush=True)
    result = dict(reference)
    result["star_index"] = str(target)
    result["star_alignment_memory_bytes"] = manifest["star_alignment_memory_bytes"]
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--catalog", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--star", default="STARlong")
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--build-memory-gib", type=int, default=48)
    parser.add_argument("--only", nargs="+")
    args = parser.parse_args()
    if args.threads < 1 or args.build_memory_gib < 1:
        parser.error("threads and build-memory-gib must be positive")
    args.star = str(Path(shutil.which(args.star) or args.star).resolve(strict=True))
    version = subprocess.check_output([args.star, "--version"], text=True).strip()
    if version != VERSION:
        parser.error(f"Expected STARlong {VERSION}; found {version!r}")
    args.output = args.output.resolve()
    args.output.mkdir(parents=True, exist_ok=True)
    catalog = json.loads(args.catalog.read_text())
    references = catalog["references"]
    selected = args.only or list(references)
    unknown = set(selected) - set(references)
    if unknown:
        parser.error(f"Unknown reference identifiers: {sorted(unknown)}")
    for reference_id in selected:
        references[reference_id] = build_one(reference_id, references[reference_id], args)
        atomic_json(args.output / "catalog.admin.json", catalog)


if __name__ == "__main__":
    main()
