"""Competitive, insert-level screening against a fixed panel of genome indexes."""
import argparse
import collections
import csv
import gzip
import html
import json
import re
import subprocess
import tempfile
from pathlib import Path

from library_sequences import fastq_records, microc_leg_span, rc
from report_palette import bar_ink

NAMES = {"ARS_UCD2": "Cattle", "hg38": "Human", "mm39": "Mouse", "dm6": "Fruit fly",
         "ecoli": "E. coli", "ce11": "C. elegans", "sacCer3": "Budding yeast",
         "TAIR10": "Arabidopsis", "GRCz12ab": "Zebrafish"}
RAW = "Reads without a reconstructed insert"
AMBIG = "Ambiguous species"
NONE = "No confident match"
TOTAL = "Total inserts"

# ggsci pal_bmj("default"), keyed by reference so subsets and reorderings keep their colors.
# https://nanx.me/ggsci/reference/pal_bmj.html
SPECIES_COLORS = {
    "ARS_UCD2": "#2a6ebb", "hg38": "#cd202c", "mm39": "#f0ab00", "dm6": "#00b2a9",
    "ecoli": "#e37222", "ce11": "#c50084", "GRCz12ab": "#7d5cc6", "TAIR10": "#69be28",
    "sacCer3": "#747678", AMBIG: "#a7aab8", NONE: "#e2e7ec",
}


def queries(raw, table, output, expected_reads, leg_length=None):
    by_read = collections.defaultdict(list)
    with gzip.open(table, "rt") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            by_read[row["read_id"]].append(row)
    scopes = []
    seen = set()
    def emit(scope, sequence):
        unit = len(scopes)
        if leg_length is None:
            output.write(f">{unit}\n{sequence}\n")
        elif scope == RAW:
            output.write(f">{2 * unit}\n{sequence}\n")
        else:
            length = microc_leg_span(len(sequence), leg_length)
            if length:
                output.write(f">{2 * unit}\n{sequence[:length]}\n")
                output.write(f">{2 * unit + 1}\n{rc(sequence[-length:])}\n")
        scopes.append(scope)
    for name, seq, _ in fastq_records(raw):
        if name in seen:
            raise ValueError("Duplicate raw read identifier")
        seen.add(name)
        rows = by_read.pop(name, [])
        if not rows:
            emit(RAW, seq)
        for row in rows:
            start, end = int(row["insert_start"]), int(row["insert_end"])
            if (len(seq) != int(row["raw_length"]) or not 0 <= start <= end <= len(seq)
                    or end - start != int(row["insert_length"])):
                raise ValueError("Original FASTQ does not match retained insert coordinates")
            if row["is_adapter_dimer"] == "1":
                continue
            if start == end:
                raise ValueError("Empty non-dimer insert")
            scope = row["scope_id"] if row["assignment_status"] == "assigned" else ""
            if not scope:
                scope = "Ambiguous index" if row["assignment_status"] == "ambiguous" else "Unassigned index"
            emit(scope, seq[start:end])
    if by_read or len(seen) != expected_reads:
        raise ValueError("Original FASTQ does not match the report's reads")
    return scopes


def parse_hit(line):
    fields = line.rstrip().split("\t")
    if int(fields[1]) & (4 | 256 | 2048):
        return None
    ops = collections.Counter()
    for length, op in re.findall(r"(\d+)([MIDNSHP=X])", fields[5]):
        ops[op] += int(length)
    aligned = sum(ops[k] for k in "M=XI")
    tags = {t[:2]: int(t[5:]) for t in fields[11:] if t.startswith(("AS:i:", "NM:i:"))}
    eligible = (aligned >= 30 and aligned >= 0.70 * len(fields[9])
                and tags["NM"] <= 0.15 * (aligned + ops["D"]))
    mq = int(fields[4])
    window = (fields[2], (int(fields[3]) - 1) // 10000) if 20 <= mq < 255 else None
    return int(fields[0]), tags["AS"], eligible, window


def species_call(candidates):
    ordered = sorted(candidates, key=lambda hit: hit[1], reverse=True)
    if not ordered or not ordered[0][2]:
        return NONE, None
    label, score, _, window = ordered[0]
    if len(ordered) > 1 and score - ordered[1][1] < 6:
        return AMBIG, None
    return label, window


def classify(scopes, hits, panel, reference, leg_length=None):
    counts = collections.defaultdict(collections.Counter)
    windows = collections.defaultdict(lambda: collections.defaultdict(set))
    end_evidence = collections.Counter()
    counts[TOTAL]
    counts[RAW]
    stride = 2 if leg_length is not None else 1
    if len(scopes) * stride != len(hits):
        raise ValueError("Missing query results")
    for i, scope in enumerate(scopes):
        label, window = species_call(hits[i * stride])
        loci = [window] if window is not None else []
        if stride == 2 and scope != RAW:
            other, other_window = species_call(hits[i * stride + 1])
            if AMBIG in (label, other):
                label, evidence = AMBIG, "ambiguous"
            elif label == other:
                evidence = "neither" if label == NONE else "concordant"
            elif NONE in (label, other):
                label = other if label == NONE else label
                evidence = "one_end"
            else:
                label, evidence = AMBIG, "discordant"
            end_evidence[evidence] += 1
            if other_window is not None:
                loci.append(other_window)
        for group in ([RAW] if scope == RAW else [TOTAL, scope]):
            counts[group][label] += 1
            if loci and label not in (AMBIG, NONE):
                windows[group][label].update(loci)
    rows = [{"scope": k, "n": sum(counts[k].values()), "counts": dict(counts[k]),
             "windows": {s: len(w) for s, w in windows[k].items()}}
            for k in [TOTAL, *sorted(counts.keys() - {TOTAL, RAW}), RAW]]
    result = {"method": "competitive-local-v1", "reference": reference, "panel": list(panel), "rows": rows}
    if leg_length is not None:
        result.update(method="competitive-microc-ends-v1", maximum_end_bp=leg_length, end_evidence=dict(end_evidence))
    return result


def align_hits(fasta, size, indexes, bowtie2, threads):
    hits = [[] for _ in range(size)]
    for species, prefix in indexes.items():
        if not size:
            break
        cmd = [bowtie2, "--very-sensitive-local", "--no-unal", "--no-hd", "-f", "-p",
               str(threads), "-x", str(prefix), "-U", str(fasta)]
        with subprocess.Popen(cmd, stdout=subprocess.PIPE, text=True) as process:
            for line in process.stdout:
                hit = parse_hit(line)
                if hit is not None:
                    idx, score, eligible, window = hit
                    hits[idx].append((species, score, eligible, window))
            if process.wait():
                raise RuntimeError(f"Species screen failed for {species}")
    return hits


def percent(n, total):
    value = 100 * n / total
    return "&lt;0.1%" if 0 < value < 0.1 else f"{value:.1f}%" if n else "0%"


def counts_table(title, scopes, species_keys, labels):
    if not scopes:
        return ""
    heads = "".join(f"<th scope='col'>{html.escape(r['scope'])}</th>" for r in scopes)
    totals = "".join(f"<td>{r['n']:,}</td>" for r in scopes)
    rows = [f"<tr class='species-screened'><th scope='row'>Screened</th>{totals}</tr>"]
    for key in species_keys:
        cells = []
        for scope in scopes:
            n = scope["counts"].get(key, 0)
            if n:
                cells.append(f"<td><span>{n:,}</span> <small>{percent(n, scope['n'])}</small></td>")
            else:
                cells.append("<td class='species-zero'>0</td>" if scope["n"] else "<td>—</td>")
        rows.append(f"<tr><th scope='row'>{labels[key]}</th>{''.join(cells)}</tr>")
    return (f"<div class='table-scroll' tabindex='0' aria-label='{title}'><table class='species-counts'>"
            f"<caption>{title}</caption><thead><tr><th scope='col'>Species</th>{heads}</tr></thead>"
            f"<tbody>{''.join(rows)}</tbody></table></div>")


def render(data):
    columns = [*data["panel"], AMBIG, NONE]
    colors = {k: SPECIES_COLORS.get(k, '#536b80') for k in columns}
    labels = {k: html.escape(NAMES.get(k, k) + (" (expected)" if k == data["reference"] else ""))
              for k in columns}

    rows = data["rows"]
    extras = {TOTAL, RAW, "Unassigned index", "Ambiguous index"}
    libraries = [r for r in rows if r["scope"] not in extras]
    whole = [r for r in rows if r["scope"] in extras]
    table_order = sorted(data["panel"], key=lambda key: key != data["reference"]) + [AMBIG, NONE]
    tables = (counts_table("Sublibraries", libraries, table_order, labels)
              + counts_table("Whole library", whole, table_order, labels))
    bars = []
    visible = set()
    for scope in libraries:
        if not scope["n"]:
            continue
        segments, descriptions = [], []
        for key in columns:
            n = scope["counts"].get(key, 0)
            if not n:
                continue
            visible.add(key)
            description = f"{labels[key]} {percent(n, scope['n'])}"
            descriptions.append(description)
            width = 100 * n / scope["n"]
            ink = bar_ink(colors[key])
            segments.append(f"<span class='species-segment' title='{description}' style='width:{width:.8f}%;"
                            f"background:{colors[key]};color:{ink}'><b>{percent(n, scope['n'])}</b></span>")
        label = html.escape(scope["scope"])
        description_list = ", ".join(descriptions)
        bars.append(f"<div class='species-row'><span class='mixture-label'>{label}</span>"
                    f"<div class='species-bar' role='img' aria-label='{label}, {description_list}'>"
                    + "".join(segments) + "</div></div>")
    legend = "".join(f"<li><i aria-hidden='true' style='background:{colors[k]}'></i>{labels[k]}</li>"
                     for k in columns if k in visible)
    chart = ("".join(bars) + f"<ul class='species-legend'>{legend}</ul>" if bars
             else "<p class='muted'>No inserts assigned to a sublibrary.</p>")
    microc_note = ""
    if data.get("method") == "competitive-microc-ends-v1":
        ends = data["end_evidence"]
        microc_note = (f"<p>Micro-C inserts are screened from both ends (up to {int(data['maximum_end_bp'])} bp per end). "
                       f"Both ends agreed in {ends.get('concordant', 0):,} inserts. One end assigned in "
                       f"{ends.get('one_end', 0):,}. Ends assigned to different species ({ends.get('discordant', 0):,}) "
                       "and ambiguous end matches count as ambiguous inserts. Each insert counts once.</p>")
    windows = " · ".join(html.escape(NAMES.get(k,k)) + f" {rows[0]['windows'].get(k,0):,}" for k in data["panel"])
    return ("<section class='report-section' id='species-composition'><h2>Species composition</h2>"
            "<style>.species-row{display:grid;grid-template-columns:minmax(5em,18%) minmax(0,1fr);"
            "align-items:center;gap:1em;margin:.65em 0}.species-row .mixture-label{overflow-wrap:anywhere}"
            ".species-bar{display:flex;height:1.2em;border-radius:4px;overflow:hidden}"
            ".species-segment{container:species-segment / inline-size;display:grid;place-items:center;overflow:hidden}"
            ".species-segment b{display:none;font:600 11px/1 system-ui;white-space:nowrap}"
            "@container species-segment (min-width:46px){.species-segment b{display:inline}}"
            ".species-legend{display:flex;flex-wrap:wrap;gap:.35em 1em;padding:0;list-style:none}"
            ".species-legend li{display:inline-flex;align-items:center;gap:.4em}"
            ".species-legend i{width:.85em;height:.85em;border-radius:2px}"
            ".species-bar,.species-legend i{print-color-adjust:exact;-webkit-print-color-adjust:exact}</style>"
            + chart + "<details><summary>Counts and whole-library totals</summary>"
            "<p>Counts and percentages use the screened total in each column. Total inserts excludes adapter dimers. "
            "Reads without a reconstructed insert are counted separately.</p>"
            + tables + microc_note + "<p>Distinct 10 kb windows with a "
            "confidently placed total-library insert (MAPQ ≥20): " + windows +
            ".</p><p>Bowtie 2 compares local alignments across the listed genomes. A match requires at least "
            "30 aligned bases, 70% query coverage, 85% identity and a 6-point lead over the next genome. "
            "Close matches are ambiguous. RNA splice junctions can reduce assignment. A species match can also reflect a spike-in. Raw-read matches alone cannot identify free genomic DNA.</p></details></section>")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--raw", nargs="+", required=True)
    for key in ("table", "summary", "reference", "indexes", "output"):
        parser.add_argument("--"+key, required=True)
    parser.add_argument("--threads", type=int, default=2)
    parser.add_argument("--bowtie2", default="bowtie2")
    parser.add_argument("--microc-leg-length", type=int)
    args = parser.parse_args()
    if args.microc_leg_length is not None and args.microc_leg_length <= 0:
        parser.error("Micro-C leg length must be positive")
    summary = json.loads(Path(args.summary).read_text())
    leg_length = args.microc_leg_length
    if leg_length is None and summary.get("assay") == "microc":
        leg_length = 100
    target = Path(args.output)
    target.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(dir=target.parent) as tmp:
        fasta = Path(tmp) / "queries.fa"
        with fasta.open("w") as out:
            scopes = queries(args.raw, args.table, out, summary["raw"]["count"], leg_length)
        indexes = json.loads(args.indexes)
        if args.reference not in indexes or len(indexes) < 2:
            raise ValueError("Species panel needs the selected reference and a competitor")
        size = len(scopes) * (2 if leg_length is not None else 1)
        hits = align_hits(fasta, size, indexes, args.bowtie2, args.threads)
        result = classify(scopes, hits, indexes, args.reference, leg_length)
        pending = Path(tmp) / "aggregate.json"
        pending.write_text(json.dumps(result, indent=2) + "\n")
        pending.replace(target)


if __name__ == "__main__":
    main()
