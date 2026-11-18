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

from library_sequences import fastq_records

NAMES = {"ARS_UCD2": "Cattle", "hg38": "Human", "mm39": "Mouse", "dm6": "Fruit fly",
         "ecoli": "E. coli", "ce11": "C. elegans", "sacCer3": "Budding yeast",
         "TAIR10": "Arabidopsis", "GRCz12ab": "Zebrafish"}
RAW = "Reads without a reconstructed insert"
AMBIG = "Ambiguous species"
NONE = "No confident match"
TOTAL = "Total inserts"


def queries(raw, table, output, expected_reads):
    by_read = collections.defaultdict(list)
    with gzip.open(table, "rt") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            by_read[row["read_id"]].append(row)
    scopes = []
    seen = set()
    def emit(scope, sequence):
        output.write(f">{len(scopes)}\n{sequence}\n")
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


def classify(scopes, hits, panel, reference):
    counts = collections.defaultdict(collections.Counter)
    windows = collections.defaultdict(lambda: collections.defaultdict(set))
    counts[TOTAL]
    counts[RAW]
    if len(scopes) != len(hits):
        raise ValueError("Missing query results")
    for scope, candidates in zip(scopes, hits):
        ordered = sorted(candidates, key=lambda hit: hit[1], reverse=True)
        label, window = NONE, None
        if ordered and ordered[0][2]:
            label, score, _, window = ordered[0]
            if len(ordered) > 1 and score - ordered[1][1] < 6:
                label, window = AMBIG, None
        for group in ([RAW] if scope == RAW else [TOTAL, scope]):
            counts[group][label] += 1
            if window is not None:
                windows[group][label].add(window)
    rows = [{"scope": k, "n": sum(counts[k].values()), "counts": dict(counts[k]),
             "windows": {s: len(w) for s, w in windows[k].items()}}
            for k in [TOTAL, *sorted(counts.keys() - {TOTAL, RAW}), RAW]]
    return {"method": "competitive-local-v1", "reference": reference, "panel": list(panel), "rows": rows}


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


def render(data):
    columns = [*data["panel"], AMBIG, NONE]
    colors = {k: color for k, color in zip(data["panel"], ["#296ab5", "#18857e", "#8862a9", "#bf6b2d", "#587b49", "#985679", "#417a92", "#715a49"])}
    colors.update({AMBIG: "#a7aab8", NONE: "#e2e7ec"})
    labels = {k: html.escape(NAMES.get(k, k) + (" (expected)" if k == data["reference"] else ""))
              for k in columns}

    def percent(n, total):
        value = 100 * n / total
        return "&lt;0.1%" if 0 < value < 0.1 else f"{value:.1f}%" if n else "0%"

    def row(r):
        cells = []
        for species in columns:
            n = r["counts"].get(species, 0)
            value = "—"
            if r["n"]:
                value = f"{n:,} ({percent(n, r['n'])})"
            cells.append(f"<td>{value}</td>")
        return "<tr><th scope='row'>" + html.escape(r["scope"]) + f"</th><td>{r['n']:,}</td>" + "".join(cells) + "</tr>"
    rows = data["rows"]
    bars = []
    visible = set()
    for scope in rows:
        if scope["scope"] in {TOTAL, RAW, "Unassigned index", "Ambiguous index"} or not scope["n"]:
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
            segments.append(f"<span title='{description}' style='width:{width:.8f}%;background:{colors.get(key, '#536b80')}'></span>")
        label = html.escape(scope["scope"])
        description_list = ", ".join(descriptions)
        bars.append(f"<div class='species-row'><span class='mixture-label'>{label}</span>"
                    f"<div class='species-bar' role='img' aria-label='{label}, {description_list}'>"
                    + "".join(segments) + "</div></div>")
    legend = "".join(f"<li><i aria-hidden='true' style='background:{colors.get(k, '#536b80')}'></i>{labels[k]}</li>"
                     for k in columns if k in visible)
    chart = ("".join(bars) + f"<ul class='species-legend'>{legend}</ul>" if bars
             else "<p class='muted'>No inserts assigned to a sublibrary.</p>")
    windows = " · ".join(html.escape(NAMES.get(k,k)) + f" {rows[0]['windows'].get(k,0):,}" for k in data["panel"])
    return ("<section class='report-section' id='species-composition'><h2>Species composition</h2>"
            "<style>.species-row{display:grid;grid-template-columns:minmax(5em,18%) minmax(0,1fr);"
            "align-items:center;gap:1em;margin:.65em 0}.species-row .mixture-label{overflow-wrap:anywhere}"
            ".species-bar{display:flex;height:1.2em;border-radius:4px;overflow:hidden}"
            ".species-legend{display:flex;flex-wrap:wrap;gap:.35em 1em;padding:0;list-style:none}"
            ".species-legend li{display:inline-flex;align-items:center;gap:.4em}"
            ".species-legend i{width:.85em;height:.85em;border-radius:2px}"
            ".species-bar,.species-legend i{print-color-adjust:exact;-webkit-print-color-adjust:exact}</style>"
            + chart + "<details><summary>Counts and whole-library totals</summary>"
            "<p>Percentages are within each row. Total inserts includes reconstructed, non-dimer inserts. "
            "Reads without a reconstructed insert are counted separately.</p>"
            "<div class='table-scroll'><table><thead><tr><th>Library</th><th>Screened</th>" +
            "".join(f"<th scope='col'>{labels[k]}</th>" for k in columns) + "</tr></thead><tbody>" +
            "".join(map(row, rows)) + "</tbody></table></div><p>Distinct 10 kb windows with a "
            "confidently placed total-library insert (MAPQ ≥20): " + windows +
            ".</p><p>Bowtie 2 compares local alignments across the listed genomes. A match requires at least "
            "30 aligned bases, 70% query coverage, 85% identity and a 6-point lead over the next genome. "
            "Close matches are ambiguous. RNA splice junctions and Micro-C ligation junctions can reduce assignment. A species match can also reflect a spike-in. Raw-read matches alone cannot identify free genomic DNA.</p></details></section>")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--raw", nargs="+", required=True)
    for key in ("table", "summary", "reference", "indexes", "output"):
        parser.add_argument("--"+key, required=True)
    parser.add_argument("--threads", type=int, default=2)
    parser.add_argument("--bowtie2", default="bowtie2")
    args = parser.parse_args()
    target = Path(args.output)
    target.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(dir=target.parent) as tmp:
        fasta = Path(tmp) / "queries.fa"
        with fasta.open("w") as out:
            scopes = queries(args.raw, args.table, out, json.loads(Path(args.summary).read_text())["raw"]["count"])
        indexes = json.loads(args.indexes)
        if args.reference not in indexes or len(indexes) < 2:
            raise ValueError("Species panel needs the selected reference and a competitor")
        hits = align_hits(fasta, len(scopes), indexes, args.bowtie2, args.threads)
        result = classify(scopes, hits, indexes, args.reference)
        pending = Path(tmp) / "aggregate.json"
        pending.write_text(json.dumps(result, indent=2) + "\n")
        pending.replace(target)


if __name__ == "__main__":
    main()
