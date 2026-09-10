#!/usr/bin/env python3
"""Dependency-light ONT QC for finished Illumina chromatin and Micro-C libraries."""

import argparse
import bisect
import collections
import gzip
import html
import json
import math
import datetime
import os
import statistics
import subprocess
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from qc_definitions import profile as metric_profile, configured_params, metric_settings, summarize_contacts, contact_evidence_status
from qc_plot_metrics import normalized_collision_stats, collision_metric_fields, svg_size_trace, size_trace_data
from spikein_auto import scoped_summary as auto_spikein_summary
from qc_reference_guide import context as reference_context, reference_notes_html


def rc(seq):
    return seq.translate(str.maketrans("ACGTNacgtn", "TGCANtgcan"))[::-1]


P5 = "AATGATACGGCGACCACCGAGATCTACAC"
P7 = "CAAGCAGAAGACGGCATACGAGAT"
TRU_R1 = "ACACTCTTTCCCTACACGACGCTCTTCCGATCT"
TRU_R2 = "GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT"
NEX_R1 = "TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG"
NEX_R2 = "GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG"

FAMILIES = {
    "truseq": (TRU_R1, TRU_R2),
    "nextera": (NEX_R1, NEX_R2),
}


def open_text(path, mode="rt"):
    return gzip.open(path, mode) if str(path).endswith(".gz") else open(path, mode)


def fastq_records(paths):
    for path in paths:
        with open_text(path) as handle:
            while True:
                name = handle.readline().rstrip()
                if not name:
                    break
                seq = handle.readline().rstrip()
                plus = handle.readline().rstrip()
                qual = handle.readline().rstrip()
                if not name.startswith("@") or not plus.startswith("+") or len(seq) != len(qual):
                    raise ValueError(f"Malformed FASTQ record in {path}: {name[:80]}")
                yield name[1:].split()[0], seq.upper(), qual


def edit_distance(a,b,limit=None):
    """Exact bit-vector Levenshtein; bounded calls saturate above the cutoff."""
    if limit is not None and abs(len(a)-len(b))>limit:
        return limit+1
    if not a or not b:
        value=max(len(a),len(b))
        return min(value,limit+1) if limit is not None else value
    masks={}
    for i,c in enumerate(a):
        masks[c]=masks.get(c,0)|(1<<i)
    positive,negative,score=~0,0,len(a)
    last=1<<(len(a)-1)
    for c in b:
        equal=masks.get(c,0)
        xv=equal|negative
        xh=(((equal&positive)+positive)^positive)|equal
        ph=negative|~(xh|positive)
        mh=positive&xh
        score+=bool(ph&last)-bool(mh&last)
        ph=(ph<<1)|1
        mh <<= 1
        positive=mh|~(xv|ph)
        negative=ph&xv
    return min(score,limit+1) if limit is not None else score

def approximate_matches(seq, pattern, error_fraction=0.20):
    """Return non-overlapping local matches using seed votes plus edit-distance correction."""
    m = len(pattern)
    max_edit = max(2, int(math.floor(m * error_fraction)))
    candidate_votes = collections.Counter()
    for k in (9, 7, 5):
        for offset in range(0, m - k + 1, max(3, k // 2)):
            seed = pattern[offset:offset + k]
            pos = seq.find(seed)
            while pos >= 0:
                candidate_votes[pos - offset] += 1
                pos = seq.find(seed, pos + 1)
        if candidate_votes:
            break
    if not candidate_votes:
        for pos in range(max(0, len(seq) - m + 1)):
            mismatches = sum(x != y for x, y in zip(pattern, seq[pos:pos + m]))
            if mismatches <= max_edit:
                candidate_votes[pos] += 1

    clusters = []
    for pos, votes in sorted(candidate_votes.items()):
        if clusters and pos - clusters[-1][0][-1] <= 4:
            clusters[-1][0].append(pos)
            clusters[-1][1] += votes
        else:
            clusters.append([[pos], votes])
    candidates = []
    for positions, votes in clusters:
        center = statistics.median(positions)
        chosen = max(positions, key=lambda x: (candidate_votes[x], -abs(x - center)))
        candidates.append((chosen, votes))

    scored = []
    for expected, votes in candidates:
        quick = []
        for pos in range(max(0, expected - 3), min(len(seq) - m + 1, expected + 4)):
            quick.append((sum(x != y for x, y in zip(pattern, seq[pos:pos + m])), pos))
        if not quick:
            continue
        hamming, quick_start = min(quick)
        if hamming > max_edit + 2:
            # Hamming distance is not an indel-safe exclusion criterion.
            # Retain the seed-voted start for bounded edit-distance refinement.
            quick_start = max(0, expected)
        best = None
        for pos in range(max(0, quick_start - 1), min(len(seq), quick_start + 2)):
            for length in range(max(12, m - 2), m + 3):
                stop = pos + length
                if stop > len(seq):
                    continue
                distance = edit_distance(pattern, seq[pos:stop], max_edit)
                key = (distance, abs(length - m), pos, stop)
                if best is None or key < best:
                    best = key
        if best and best[0] <= max_edit:
            scored.append({"start": best[2], "end": best[3], "edits": best[0],
                           "identity": 1.0 - best[0] / max(m, best[3] - best[2])})
    scored.sort(key=lambda x: (x["edits"], -x["identity"], x["start"]))
    kept = []
    for hit in scored:
        if any(not (hit["end"] <= old["start"] + 5 or hit["start"] >= old["end"] - 5) for old in kept):
            continue
        kept.append(hit)
    return sorted(kept, key=lambda x: x["start"])


def matches_in_window(seq, pattern, lo, hi):
    lo, hi = max(0, lo), min(len(seq), hi)
    hits = approximate_matches(seq[lo:hi], pattern)
    for hit in hits:
        hit["start"] += lo
        hit["end"] += lo
    return hits

def find_units(seq, max_unit):
    all_units = []
    matches_by_family = {}
    for family, (r1, r2) in FAMILIES.items():
        hits = {
            "r1f": approximate_matches(seq, r1),
            "r1r": approximate_matches(seq, rc(r1)),
            "r2f": approximate_matches(seq, r2),
            "r2r": approximate_matches(seq, rc(r2)),
        }
        matches_by_family[family] = hits
        for orientation, left_key, right_key in (("p5_to_p7", "r1f", "r2r"),
                                                  ("p7_to_p5", "r2f", "r1r")):
            used_right = set()
            for left in hits[left_key]:
                options = [(i, right) for i, right in enumerate(hits[right_key])
                           if i not in used_right and 0 <= right["start"] - left["end"] <= max_unit]
                if not options:
                    continue
                i, right = min(options, key=lambda x: x[1]["start"])
                used_right.add(i)
                all_units.append({
                    "family": family, "orientation": orientation,
                    "left": left, "right": right,
                    "start": left["end"], "end": right["start"],
                    "adapter_edits": left["edits"] + right["edits"],
                })
    # Prefer the best-scoring architecture when two families claim the same interval.
    all_units.sort(key=lambda u: (u["adapter_edits"], -(u["end"] - u["start"])))
    kept = []
    for unit in all_units:
        overlap = False
        for old in kept:
            intersect = max(0, min(unit["end"], old["end"]) - max(unit["start"], old["start"]))
            if max(unit["left"]["start"], old["left"]["start"]) < min(unit["right"]["end"], old["right"]["end"]):
                overlap = True
                break
        if not overlap:
            kept.append(unit)
    return sorted(kept, key=lambda u: u["start"]), matches_by_family


def nearest_before(hits, coordinate, max_gap=18):
    valid = [h for h in hits if 0 <= coordinate - h["end"] <= max_gap]
    return max(valid, key=lambda h: h["end"]) if valid else None


def nearest_after(hits, coordinate, max_gap=18):
    valid = [h for h in hits if 0 <= h["start"] - coordinate <= max_gap]
    return min(valid, key=lambda h: h["start"]) if valid else None


def add_indexes(seq, units):
    for unit in units:
        left, right = unit["left"], unit["right"]
        i5_obs = i7_obs = None
        p5_hit = p7_hit = None
        if unit["orientation"] == "p5_to_p7":
            p5_hit = nearest_before(matches_in_window(seq, P5, left["start"] - len(P5) - 20, left["start"]), left["start"])
            p7_hit = nearest_after(matches_in_window(seq, rc(P7), right["end"], right["end"] + len(P7) + 20), right["end"])
            if p5_hit:
                i5_obs = seq[p5_hit["end"]:left["start"]]
            if p7_hit:
                i7_obs = seq[right["end"]:p7_hit["start"]]
            i5 = i5_obs
            i7 = i7_obs
        else:
            p7_hit = nearest_before(matches_in_window(seq, P7, left["start"] - len(P7) - 20, left["start"]), left["start"])
            p5_hit = nearest_after(matches_in_window(seq, rc(P5), right["end"], right["end"] + len(P5) + 20), right["end"])
            if p7_hit:
                i7_obs = seq[p7_hit["end"]:left["start"]]
            if p5_hit:
                i5_obs = seq[right["end"]:p5_hit["start"]]
            i5 = rc(i5_obs) if i5_obs is not None else None
            i7 = rc(i7_obs) if i7_obs is not None else None
        unit.update({"p5_hit": p5_hit, "p7_hit": p7_hit,
                     "i5_observed": i5_obs, "i7_observed": i7_obs,
                     "i5": i5, "i7": i7})

def quantile(values, p):
    if not values:
        return None
    x = sorted(values)
    pos = (len(x) - 1) * p
    lo, hi = int(math.floor(pos)), int(math.ceil(pos))
    return x[lo] if lo == hi else x[lo] * (hi - pos) + x[hi] * (pos - lo)


def n50(values):
    if not values:
        return None
    target = sum(values) / 2
    running = 0
    for value in sorted(values, reverse=True):
        running += value
        if running >= target:
            return value


def length_summary(values):
    return {
        "count": len(values), "total_bases": sum(values),
        "min": min(values) if values else None,
        "p10": quantile(values, .10), "median": quantile(values, .50),
        "p90": quantile(values, .90), "p99": quantile(values, .99),
        "max": max(values) if values else None, "n50": n50(values),
    }


def index_consensus(values):
    values = [v for v in values if v is not None and len(v) <= 14]
    if not values:
        return {"sequence": None, "length": None, "support": 0, "purity": None, "status": "not_observed"}
    length_counts = collections.Counter(map(len, values))
    chosen_length, length_support = length_counts.most_common(1)[0]
    same = [v for v in values if len(v) == chosen_length]
    if chosen_length == 0:
        consensus = ""
    else:
        consensus = "".join(collections.Counter(v[i] for v in same).most_common(1)[0][0]
                            for i in range(chosen_length))
    exact = sum(v == consensus for v in same)
    allowed = chosen_length in {0, 6, 8, 10}
    return {"sequence": consensus, "length": chosen_length, "support": len(same),
            "all_flanked_observations": len(values), "purity": exact / len(same),
            "status": "consensus" if allowed else "custom_or_indel_uncertain"}



def parse_expected_index_pairs(text):
    pairs = []
    if not text or text in {".", "none", "None"}:
        return pairs
    for item in text.split(";"):
        fields = item.split("|")
        if len(fields) != 3:
            raise ValueError("Expected index pairs must use name|i5|i7;name|i5|i7")
        name, i5, i7 = fields
        pairs.append({"name": name, "i5": i5, "i7": i7})
    return pairs


def index_close(observed, expected, maximum_edit=1):
    if observed is None:
        return True
    if abs(len(observed) - len(expected)) > maximum_edit:
        return False
    return edit_distance(observed, expected, maximum_edit) <= maximum_edit


def index_pair_distance(left, right, maximum_edit=1):
    """Capped total edit distance across a physical i5/i7 pair.

    Auto-discovered seeds within the configured ONT indel-rescue window can be
    treated as variants of a much stronger seed. Explicit ``index_pairs`` remain
    authoritative when intentionally similar indexes are used.
    """
    total = 0
    for key in ("i5", "i7"):
        remaining = maximum_edit - total
        total += edit_distance(left[key], right[key], remaining)
        if total > maximum_edit:
            break
    return total


def cluster_index_pairs(units, single_fraction=.90, minimum_informative=20,
                        minimum_cluster_units=20, minimum_cluster_fraction=.01,
                        expected_text="", variant_max_edit=2,
                        variant_dominance_ratio=5.0):
    informative = [u for u in units if u["p5_detected"] and u["p7_detected"]]
    exact = collections.Counter((u["i5"], u["i7"]) for u in informative
                                if u["i5"] is not None and u["i7"] is not None
                                and len(u["i5"]) <= 14 and len(u["i7"]) <= 14)
    expected = parse_expected_index_pairs(expected_text)
    seeds = []
    suppressed_variants = []
    if expected:
        seeds = [{"name": x["name"], "i5": x["i5"], "i7": x["i7"], "seed_support": 0}
                 for x in expected]
    else:
        support_floor = max(minimum_cluster_units,
                            math.ceil(max(1, len(informative)) * minimum_cluster_fraction))
        for (i5, i7), count in sorted(exact.items(), key=lambda x: (-x[1], x[0])):
            if count < support_floor:
                continue
            candidate = {"name": None, "i5": i5, "i7": i7, "seed_support": count}
            # ONT errors can create short deletions or longer insertions. A
            # low-count exact pair is a variant of an established pair only
            # when it is close in edit distance and the established pair has
            # substantially greater support. The support-ratio gate prevents
            # collapsing two genuine, similarly abundant physical indexes.
            dominated_by = None
            for seed in seeds:
                distance = index_pair_distance(candidate, seed, variant_max_edit)
                if (distance <= variant_max_edit
                        and seed.get("seed_support", 0) >= variant_dominance_ratio * count):
                    dominated_by = (seed, distance)
                    break
            if dominated_by is not None:
                seed, distance = dominated_by
                suppressed_variants.append({
                    "i5": i5, "i7": i7, "support": count,
                    "dominant_i5": seed["i5"], "dominant_i7": seed["i7"],
                    "dominant_support": seed.get("seed_support", 0),
                    "pair_edit_distance": distance,
                })
                continue
            seeds.append(candidate)
        if not seeds and exact:
            (i5, i7), count = exact.most_common(1)[0]
            seeds.append({"name": None, "i5": i5, "i7": i7, "seed_support": count})
    seeds.sort(key=lambda x: (-x["seed_support"], x["i5"], x["i7"]))
    for number, seed in enumerate(seeds, 1):
        seed["scope_id"] = f"lib{number:02d}"
        seed["name"] = seed["name"] or seed["scope_id"]
        seed["unit_ids"] = []
        seed["informative_support"] = 0
        seed["mappable_units"] = 0
        seed["adapter_dimer_units"] = 0
        seed["raw_read_ids"] = set()
    ambiguous = unassigned = 0
    assignments = []
    for unit in units:
        candidates = []
        for seed in seeds:
            if (not unit["p5_detected"] or index_close(unit["i5"], seed["i5"], variant_max_edit)) and \
               (not unit["p7_detected"] or index_close(unit["i7"], seed["i7"], variant_max_edit)):
                d5 = 0 if not unit["p5_detected"] else edit_distance(unit["i5"], seed["i5"], variant_max_edit)
                d7 = 0 if not unit["p7_detected"] else edit_distance(unit["i7"], seed["i7"], variant_max_edit)
                candidates.append((d5 + d7, seed))
        status, assigned = "unassigned", None
        if candidates:
            best_distance = min(x[0] for x in candidates)
            best = [seed for distance, seed in candidates if distance == best_distance]
            if len(best) == 1:
                assigned, status = best[0], "assigned"
            else:
                status = "ambiguous"
        if assigned is None:
            ambiguous += status == "ambiguous"
            unassigned += status == "unassigned"
        else:
            assigned["unit_ids"].append(unit["unit_id"])
            assigned["informative_support"] += int(unit["p5_detected"] and unit["p7_detected"])
            assigned["mappable_units"] += int(not unit["is_adapter_dimer"])
            assigned["adapter_dimer_units"] += int(unit["is_adapter_dimer"])
            assigned["raw_read_ids"].add(unit["read_id"])
        assignments.append({"unit_id": unit["unit_id"], "read_id": unit["read_id"],
                            "scope_id": assigned["scope_id"] if assigned else "",
                            "assignment_status": status, "i5": unit["i5"], "i7": unit["i7"],
                            "p5_detected": unit["p5_detected"], "p7_detected": unit["p7_detected"],
                            "is_adapter_dimer": unit["is_adapter_dimer"]})
    informative_n = len(informative)
    assigned_informative = sum(seed["informative_support"] for seed in seeds)
    top_fraction = max([seed["informative_support"] for seed in seeds] + [0]) / informative_n if informative_n else None
    supported = [seed for seed in seeds if seed["informative_support"] >= max(
        minimum_cluster_units, math.ceil(max(1, informative_n) * minimum_cluster_fraction))]
    if informative_n < minimum_informative:
        mode = "inconclusive_low_index_information"
    elif top_fraction is not None and top_fraction >= single_fraction and seeds:
        mode = "single_library"
    elif len(seeds) == 1 and len(supported) == 1:
        mode = "single_library"
    elif len(supported) >= 2:
        mode = "multiplexed_library"
    else:
        mode = "mixed_or_unresolved"
    public_pairs = []
    for seed in seeds:
        public_pairs.append({
            "scope_id": seed["scope_id"], "name": seed["name"], "i5": seed["i5"], "i7": seed["i7"],
            "support_units": len(seed["unit_ids"]), "informative_support": seed["informative_support"],
            "support_fraction": seed["informative_support"] / informative_n if informative_n else None,
            "mappable_units": seed["mappable_units"], "adapter_dimer_units": seed["adapter_dimer_units"],
            "raw_reads": len(seed["raw_read_ids"]), "unit_ids": seed["unit_ids"],
        })
    result = {"mode": mode, "single_library_fraction": single_fraction,
              "minimum_informative_units": minimum_informative,
              "informative_units": informative_n, "assigned_informative_units": assigned_informative,
              "assignment_rate": assigned_informative / informative_n if informative_n else None,
              "all_unit_assignment_rate": sum(bool(a["scope_id"]) for a in assignments) / len(units) if units else None,
              "ambiguous_fraction": ambiguous / len(units) if units else None,
              "unassigned_fraction": unassigned / len(units) if units else None,
              "top_pair_fraction": top_fraction, "ambiguous_units": ambiguous,
              "unassigned_units": unassigned, "pairs": public_pairs,
              "index_variant_max_edit": variant_max_edit,
              "index_variant_dominance_ratio": variant_dominance_ratio,
              "suppressed_variant_pairs": suppressed_variants,
              "sequence_orientation": "physical P5-side i5 / P7-side i7, normalized P5-to-P7"}
    return result, assignments

def zero_insert_template(family, i5, i7):
    r1, r2 = FAMILIES[family]
    return P5 + (i5 or "") + r1 + rc(r2) + (i7 or "") + rc(P7)


def full_length_template_identity(sequence, template):
    best = None
    for target in (template, rc(template)):
        limit = max(3, math.ceil(max(len(sequence), len(target)) * .20))
        distance = edit_distance(sequence, target, limit)
        identity = 1 - distance / max(len(sequence), len(target))
        best = identity if best is None else max(best, identity)
    return best


def json_dump(data, path):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as handle:
        json.dump(data, handle, indent=2, sort_keys=True)


def mean_error_q(quality):
    """Phred transform of mean base error probability (all bases, no trimming)."""
    return -10 * math.log10(sum(10 ** (-(ord(x) - 33) / 10) for x in quality) / len(quality)) if quality else None


def cmd_prepare(args):
    for path in (args.inserts, args.legs, args.leg1, args.leg2, args.read_table, args.index_assignments,
                 args.summary, args.raw_lengths, args.insert_lengths):
        Path(path).parent.mkdir(parents=True, exist_ok=True)
    raw_lengths, insert_lengths, mean_q = [], [], []
    families, orientations = collections.Counter(), collections.Counter()
    i5_values, i7_values, units_all, raw_all = [], [], [], []
    read_count = complete_reads = concatemer_reads = unit_count = 0
    seen_read_ids = set()
    for read_id, seq, qual in fastq_records(args.fastq):
        if read_id in seen_read_ids:
            raise ValueError(f"Duplicate FASTQ read ID across input files: {read_id}")
        seen_read_ids.add(read_id)
        read_count += 1
        raw_lengths.append(len(seq))
        q = mean_error_q(qual) if qual else 0
        mean_q.append(q)
        raw_all.append({"read_id": read_id, "seq": seq, "qual": qual, "mean_q": q})
        units, _ = find_units(seq, args.max_unit)
        add_indexes(seq, units)
        complete_reads += int(bool(units)); concatemer_reads += int(len(units) > 1)
        for number, unit in enumerate(units, 1):
            unit_count += 1
            unit_id = f"{read_id}|u{number}"
            ins_seq = seq[unit["start"]:unit["end"]]
            ins_qual = qual[unit["start"]:unit["end"]]
            insert_lengths.append(len(ins_seq))
            canonical_dimer = len(ins_seq) <= args.adapter_dimer_max_bp
            families[unit["family"]] += 1; orientations[unit["orientation"]] += 1
            i5_values.append(unit["i5"]); i7_values.append(unit["i7"])
            units_all.append({
                "read_id": read_id, "unit_id": unit_id, "raw_length": len(seq), "raw_seq": seq,
                "mean_q": q, "family": unit["family"], "orientation": unit["orientation"],
                "insert_start": unit["start"], "insert_end": unit["end"], "insert_length": len(ins_seq),
                "insert_seq": ins_seq, "insert_qual": ins_qual, "adapter_edits": unit["adapter_edits"],
                "i5_observed": unit["i5_observed"], "i5": unit["i5"],
                "i7_observed": unit["i7_observed"], "i7": unit["i7"],
                "p5_detected": unit["p5_hit"] is not None, "p7_detected": unit["p7_hit"] is not None,
                "is_adapter_dimer": canonical_dimer, "dimer_class": "canonical_adapter_dimer" if canonical_dimer else "not_dimer",
                "template_identity": None,
            })
    pairing, assignments = cluster_index_pairs(
        units_all, args.single_index_pair_fraction, args.minimum_index_informative_units,
        args.minimum_index_cluster_units, args.minimum_index_cluster_fraction, args.expected_index_pairs,
        args.index_variant_max_edit, args.index_variant_dominance_ratio)
    unit_scope = {row["unit_id"]: row["scope_id"] for row in assignments if row["scope_id"]}
    pair_by_scope = {pair["scope_id"]: pair for pair in pairing["pairs"]}
    raw_dimer_reads = set()
    for unit in units_all:
        scope_id = unit_scope.get(unit["unit_id"])
        pair = pair_by_scope.get(scope_id)
        if pair and args.adapter_dimer_raw_min_bp <= unit["raw_length"] <= args.adapter_dimer_raw_max_bp:
            template = zero_insert_template(unit["family"], pair["i5"], pair["i7"])
            unit["template_identity"] = full_length_template_identity(unit["raw_seq"], template)
            # A template can match a true 11-50 bp insert at >=80% identity.
            # Once both inner anchors resolve a nonzero insert, retain that evidence.
        if not unit["is_adapter_dimer"] and unit["insert_length"] <= args.ultrashort_insert_max_bp:
            unit["dimer_class"] = "ultrashort_insert"
        if unit["is_adapter_dimer"]:
            raw_dimer_reads.add(unit["read_id"])
    # Update pair counts after structure-aware dimer classification.
    unit_by_id = {u["unit_id"]: u for u in units_all}
    for pair in pairing["pairs"]:
        pair_units = [unit_by_id[x] for x in pair["unit_ids"]]
        pair["adapter_dimer_units"] = sum(u["is_adapter_dimer"] for u in pair_units)
        pair["mappable_units"] = sum(not u["is_adapter_dimer"] and (args.assay != "microc" or u["insert_length"] >= 50) for u in pair_units)
        pair["adapter_dimer_raw_reads"] = len({u["read_id"] for u in pair_units if u["is_adapter_dimer"]})
    # Search short unreconstructed reads against supported zero-insert templates.
    templates = []
    family = max(families, key=families.get) if families else "truseq"
    for pair in pairing["pairs"]:
        templates.append((pair["scope_id"], zero_insert_template(family, pair["i5"], pair["i7"])))
    unreconstructed_probable = 0
    reconstructed_reads = {u["read_id"] for u in units_all}
    for raw in raw_all:
        if raw["read_id"] in reconstructed_reads or not (args.adapter_dimer_raw_min_bp <= len(raw["seq"]) <= args.adapter_dimer_raw_max_bp):
            continue
        if templates and max(full_length_template_identity(raw["seq"], template) for _, template in templates) >= args.adapter_dimer_template_identity:
            unreconstructed_probable += 1
    adapter_dimer_units = sum(u["is_adapter_dimer"] for u in units_all)
    canonical_dimers = sum(u["dimer_class"] == "canonical_adapter_dimer" for u in units_all)
    probable_dimers = sum(u["dimer_class"] == "probable_adapter_dimer" for u in units_all)
    ultrashort_units = sum(u["dimer_class"] == "ultrashort_insert" for u in units_all)
    assignment_by_unit = {a["unit_id"]: a for a in assignments}
    table_header = ["read_id", "raw_length", "mean_q", "unit_id", "scope_id", "assignment_status",
                    "family", "orientation",
                    "insert_start", "insert_end", "insert_length", "adapter_edits", "is_adapter_dimer",
                    "dimer_class", "dimer_template_identity", "p5_detected", "p7_detected",
                    "i5_observed", "i5_normalized", "i7_observed", "i7_normalized"]
    with gzip.open(args.inserts, "wt") as insert_out, gzip.open(args.legs, "wt") as leg_out, \
            gzip.open(args.leg1, "wt") as leg1_out, gzip.open(args.leg2, "wt") as leg2_out, \
            gzip.open(args.read_table, "wt") as table:
        table.write("\t".join(table_header) + "\n")
        for unit in units_all:
            assignment = assignment_by_unit[unit["unit_id"]]
            if not unit["is_adapter_dimer"]:
                insert_out.write(f"@{unit['unit_id']}\n{unit['insert_seq']}\n+\n{unit['insert_qual']}\n")
            if args.assay == "microc" and not unit["is_adapter_dimer"] and unit["insert_length"] >= 50:
                leg_len = min(args.leg_length, unit["insert_length"] // 2)
                left_seq, left_q = unit["insert_seq"][:leg_len], unit["insert_qual"][:leg_len]
                right_seq, right_q = rc(unit["insert_seq"][-leg_len:]), unit["insert_qual"][-leg_len:][::-1]
                leg_out.write(f"@{unit['unit_id']}/1\n{left_seq}\n+\n{left_q}\n@{unit['unit_id']}/2\n{right_seq}\n+\n{right_q}\n")
                leg1_out.write(f"@{unit['unit_id']}\n{left_seq}\n+\n{left_q}\n")
                leg2_out.write(f"@{unit['unit_id']}\n{right_seq}\n+\n{right_q}\n")
            row = [unit["read_id"], unit["raw_length"], f"{unit['mean_q']:.3f}", unit["unit_id"],
                   assignment["scope_id"], assignment["assignment_status"],
                   unit["family"], unit["orientation"], unit["insert_start"], unit["insert_end"],
                   unit["insert_length"], unit["adapter_edits"], int(unit["is_adapter_dimer"]),
                   unit["dimer_class"], unit["template_identity"], int(unit["p5_detected"]), int(unit["p7_detected"]),
                   unit["i5_observed"], unit["i5"], unit["i7_observed"], unit["i7"]]
            table.write("\t".join("" if x is None else str(x) for x in row) + "\n")
    assignment_header = ["unit_id", "read_id", "scope_id", "assignment_status", "i5", "i7",
                         "p5_detected", "p7_detected", "is_adapter_dimer", "dimer_class"]
    with gzip.open(args.index_assignments, "wt") as out:
        out.write("\t".join(assignment_header) + "\n")
        for unit in units_all:
            row = dict(assignment_by_unit[unit["unit_id"]], is_adapter_dimer=unit["is_adapter_dimer"],
                       dimer_class=unit["dimer_class"])
            out.write("\t".join("" if row[k] is None else str(int(row[k])) if isinstance(row[k], bool)
                                else str(row[k]) for k in assignment_header) + "\n")
    with open(args.raw_lengths, "w") as out:
        out.writelines(f"{x}\n" for x in raw_lengths)
    with open(args.insert_lengths, "w") as out:
        out.writelines(f"{x}\n" for x in insert_lengths)
    i5, i7 = index_consensus(i5_values), index_consensus(i7_values)
    if i5["length"] == 0 and i7["length"] in {6, 8, 10}:
        architecture = "classic_single_index"
    elif i5["length"] in {6, 8, 10} and i7["length"] in {6, 8, 10}:
        architecture = "physical_dual_index"
    else:
        architecture = "partial_or_custom"
    raw_window_count = sum(args.adapter_dimer_raw_min_bp <= x <= args.adapter_dimer_raw_max_bp for x in raw_lengths)
    summary = {
        "sample": args.sample, "assay": args.assay, "raw": length_summary(raw_lengths),
        "inserts": length_summary(insert_lengths), "mean_read_q": statistics.mean(mean_q) if mean_q else None,
        "complete_reads": complete_reads, "complete_read_rate": complete_reads / read_count if read_count else 0,
        "library_units": unit_count, "units_per_raw_read": unit_count / read_count if read_count else 0,
        "adapter_dimer_max_bp": args.adapter_dimer_max_bp, "adapter_dimer_units": adapter_dimer_units,
        "canonical_adapter_dimer_units": canonical_dimers, "probable_adapter_dimer_units": probable_dimers,
        "unreconstructed_probable_dimer_reads": unreconstructed_probable,
        "ultrashort_insert_units": ultrashort_units, "ultrashort_insert_max_bp": args.ultrashort_insert_max_bp,
        "adapter_dimer_unit_rate": adapter_dimer_units / unit_count if unit_count else 0,
        "adapter_dimer_rate": len(raw_dimer_reads) / read_count if read_count else 0,
        "adapter_dimer_raw_window_bp": [args.adapter_dimer_raw_min_bp, args.adapter_dimer_raw_max_bp],
        "raw_reads_in_dimer_size_window": raw_window_count,
        "raw_dimer_size_window_rate": raw_window_count / read_count if read_count else 0,
        "raw_structure_dimer_reads": len(raw_dimer_reads),
        "raw_structure_dimer_rate": len(raw_dimer_reads) / read_count if read_count else 0,
        "mappable_library_units": sum(not u["is_adapter_dimer"] and (args.assay != "microc" or u["insert_length"] >= 50) for u in units_all),
        "short_unmapped_units": sum(not u["is_adapter_dimer"] and args.assay == "microc" and u["insert_length"] < 50 for u in units_all),
        "outer_complete_units": sum(u["p5_detected"] and u["p7_detected"] for u in units_all),
        "outer_complete_unit_rate": sum(u["p5_detected"] and u["p7_detected"] for u in units_all) / unit_count if unit_count else None,
        "quality_method": "mean across reads of -10log10(mean base error probability); all bases",
        "concatemer_reads": concatemer_reads, "concatemer_rate": concatemer_reads / read_count if read_count else 0,
        "adapter_families": dict(families), "orientations": dict(orientations),
        "i5": i5, "i7": i7, "index_architecture": architecture, "index_pairing": pairing,
        "notes": [
            "Index pairs are jointly clustered after P5-to-P7 normalization; missing outer anchors are not empty indexes.",
            "Canonical dimer-like units have a central gap at or below the cutoff. Unreconstructed template matches are candidates only and excluded from dimer burden.",
            "Ultrashort 11-50-bp inserts are reported separately and retained because genuine CUT&RUN fragments may be short.",
            "The raw dimer-size window supports, but does not alone define, adapter-dimer classification."
        ]}
    json_dump(summary, args.summary)
    print(json.dumps(summary, indent=2))


def cigar_ref_span(cigar):
    number, span = "", 0
    for char in cigar:
        if char.isdigit():
            number += char
        else:
            if char in "MDN=X" and number:
                span += int(number)
            number = ""
    return span


def sam_primary_records(bam, samtools):
    process = subprocess.Popen([samtools, "view", bam], text=True, stdout=subprocess.PIPE)
    best = {}
    for line in process.stdout:
        f = line.rstrip().split("\t")
        flag = int(f[1])
        if flag & (4 | 256 | 512 | 1024 | 2048):
            continue
        mate_suffix = "/1" if flag & 64 else ("/2" if flag & 128 else "")
        key = f[0] + mate_suffix
        record = {"qname": key, "flag": flag, "chrom": f[2], "start": int(f[3]) - 1,
                  "mapq": int(f[4]), "cigar": f[5], "strand": "-" if flag & 16 else "+"}
        record["end"] = record["start"] + cigar_ref_span(record["cigar"])
        import re
        ops = [(int(n), op) for n, op in re.findall(r"(\d+)([MIDNSHP=X])", record["cigar"])]
        aligned = sum(n for n, op in ops if op in "MI=X")
        query_span = sum(n for n, op in ops if op in "MISH=X")
        record["query_coverage"] = aligned / query_span if query_span else 0
        if key not in best or record["mapq"] > best[key]["mapq"]:
            best[key] = record
    if process.wait() != 0:
        raise RuntimeError(f"samtools view failed for {bam}")
    return best


def read_chrom_sizes(path):
    sizes = {}
    with open(path) as handle:
        for line in handle:
            if line.strip():
                chrom, size = line.split()[:2]
                sizes[chrom] = int(size)
    return sizes


def auto_mitochondrial_chromosomes(chrom_sizes):
    aliases = {"chrm", "chrmt", "mt", "m", "mitochondrion_genome", "mitochondrial"}
    return {name for name in chrom_sizes if name.lower() in aliases}


def resolve_chromosome_sets(chrom_sizes, primary_arg, mito_arg):
    mito_text = (mito_arg or "auto").strip()
    if mito_text.lower() == "auto":
        mito = auto_mitochondrial_chromosomes(chrom_sizes)
    elif mito_text.lower() in {"none", "null", "off"}:
        mito = set()
    else:
        mito = {name for name in mito_text.split(",") if name in chrom_sizes}
    primary_text = (primary_arg or "auto").strip()
    if primary_text.lower() in {"", "auto", "all"}:
        primary = [name for name in chrom_sizes if name not in mito]
    else:
        primary = [name for name in primary_text.split(",") if name in chrom_sizes and name not in mito]
    if not primary:
        raise ValueError("No primary chromosomes match the reference chromosome sizes")
    if primary_text.lower() not in {"", "auto", "all"}:
        missing = set(primary_text.split(",")) - set(chrom_sizes)
        if missing:
            raise ValueError(f"Configured primary chromosomes absent from reference: {sorted(missing)}")
    return primary, mito


def merge_intervals(intervals):
    merged = []
    for start, end in sorted(intervals):
        if not merged or start > merged[-1][1]:
            merged.append([start, end])
        else:
            merged[-1][1] = max(merged[-1][1], end)
    return [(start, end) for start, end in merged]


def load_gene_annotation(path, chrom_sizes, promoter_window):
    if not path:
        return {}, {}
    genes, transcripts = collections.defaultdict(list), collections.defaultdict(list)
    with open_text(path) as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip().split("\t")
            if len(fields) < 9 or fields[0] not in chrom_sizes:
                continue
            feature = fields[2].lower()
            if feature not in {"gene", "transcript", "mrna"}:
                continue
            start, end, strand = int(fields[3]) - 1, int(fields[4]), fields[6]
            target = genes if feature == "gene" else transcripts
            target[fields[0]].append((start, end, strand))
    selected = genes if genes else transcripts
    promoter_intervals, gene_intervals = {}, {}
    for chrom, features in selected.items():
        chrom_size = chrom_sizes[chrom]
        promoters, bodies = [], []
        for start, end, strand in features:
            tss = end - 1 if strand == "-" else start
            promoters.append((max(0, tss - promoter_window), min(chrom_size, tss + promoter_window + 1)))
            bodies.append((max(0, start), min(chrom_size, end)))
        promoter_intervals[chrom] = merge_intervals(promoters)
        gene_intervals[chrom] = merge_intervals(bodies)
    return promoter_intervals, gene_intervals


def interval_overlaps(intervals, start, end):
    if not intervals:
        return False
    index = bisect.bisect_left(intervals, (end,))
    return index > 0 and intervals[index - 1][1] > start


def annotation_region_counts(records, promoters, genes, primary_names, mapq, promoter_window):
    """Return observed promoter/gene-body/distal fractions for report visualization only."""
    primary = set(primary_names)
    eligible = [record for record in records.values() if record["mapq"] >= mapq and record["chrom"] in primary]
    promoter_count = gene_count = 0
    for record in eligible:
        chrom = record["chrom"]
        if interval_overlaps(promoters.get(chrom), record["start"], record["end"]):
            promoter_count += 1
        elif interval_overlaps(genes.get(chrom), record["start"], record["end"]):
            gene_count += 1
    total = len(eligible)
    distal = total - promoter_count - gene_count
    return {
        "annotation_mapped_reads": total,
        "promoter_proximal_reads": promoter_count,
        "promoter_proximal_rate": promoter_count / total if total and (promoters or genes) else None,
        "gene_body_reads": gene_count,
        "gene_body_rate": gene_count / total if total and (promoters or genes) else None,
        "intergenic_reads": distal,
        "intergenic_rate": distal / total if total and (promoters or genes) else None,
        "promoter_window_bp": promoter_window,
        "annotation_region_definition": "mutually exclusive: promoter proximal takes precedence over gene body",
    }


def write_karyoploter(records, chrom_sizes, chromosome_names, mapq, bins_per_longest,
                      rscript, script, genes_bed, output, title, assembly, gene_window_bp=250000):
    """Render MAPQ-pass primary alignments with Bioconductor karyoploteR."""
    names = [name for name in chromosome_names if name in chrom_sizes]
    if not names:
        raise ValueError("Karyotype plotting requires at least one configured chromosome")
    longest = max(chrom_sizes[name] for name in names)
    window_size = max(1, math.ceil(longest / max(1, bins_per_longest)))
    eligible = []
    name_set = set(names)
    for record in records.values():
        if record["chrom"] not in name_set or record["mapq"] < mapq:
            continue
        midpoint = min(chrom_sizes[record["chrom"]],
                       max(1, (record["start"] + record["end"]) // 2 + 1))
        eligible.append((record["chrom"], midpoint))
    with tempfile.TemporaryDirectory(prefix="ont_karyoploter_") as directory:
        directory = Path(directory)
        chromosome_table = directory / "chromosomes.tsv"
        reads_table = directory / "reads.tsv"
        chromosome_table.write_text("chromosome\tsize\n" + "".join(
            f"{name}\t{chrom_sizes[name]}\n" for name in names))
        reads_table.write_text("chromosome\tposition\n" + "".join(
            f"{chromosome}\t{position}\n" for chromosome, position in eligible))
        result = subprocess.run([
            str(rscript), str(script), str(reads_table), str(chromosome_table),
            str(genes_bed), str(output), str(title), str(window_size), str(mapq),
            str(assembly), str(gene_window_bp),
        ], text=True, capture_output=True)
        if result.returncode != 0:
            detail = result.stderr.strip() or result.stdout.strip() or "unknown R error"
            raise RuntimeError(f"karyoploteR failed for {title}: {detail}")
    if not Path(output).exists():
        raise RuntimeError(f"karyoploteR did not create {output}")
    return {"mapped_reads": len(eligible), "window_size_bp": window_size}

def bin_stats(records, sizes, bin_size):
    return normalized_collision_stats(records, sizes, bin_size)


def wilson(success, total, z=1.96):
    if total == 0:
        return None, None
    p = success / total
    d = 1 + z * z / total
    center = (p + z * z / (2 * total)) / d
    half = z * math.sqrt(p * (1 - p) / total + z * z / (4 * total * total)) / d
    return max(0, center - half), min(1, center + half)


def read_scope_profiles(path, include_total=False):
    """Collect profiles; optional total includes unassigned/ambiguous units."""
    profiles = collections.defaultdict(lambda: {"reads": {}, "insert_lengths": [], "unit_rows": []})
    molecule_units = collections.Counter()
    with open_text(path) as handle:
        header = handle.readline().rstrip().split("\t")
        for line in handle:
            row = dict(zip(header, line.rstrip("\n").split("\t")))
            molecule_units[row["read_id"]] += 1
            scope_ids = ["total"] if include_total else []
            if row["assignment_status"] == "assigned" and row["scope_id"]:
                scope_ids.append(row["scope_id"])
            for scope_id in scope_ids:
                profile = profiles[scope_id]
                profile["reads"][row["read_id"]] = (int(row["raw_length"]), float(row["mean_q"]))
                profile["insert_lengths"].append(int(row["insert_length"]))
                profile["unit_rows"].append(row)
    result = {}
    for scope_id, profile in profiles.items():
        raw_lengths = [value[0] for value in profile["reads"].values()]
        read_q = [value[1] for value in profile["reads"].values()]
        insert_lengths = profile["insert_lengths"]
        result[scope_id] = {
            "raw_lengths": raw_lengths, "insert_lengths": insert_lengths,
            "unit_rows": profile["unit_rows"],
            "concatemer_reads": sum(molecule_units[read_id] > 1 for read_id in profile["reads"]),
            "raw": length_summary(raw_lengths), "inserts": length_summary(insert_lengths),
            "mean_read_q": statistics.mean(read_q) if read_q else None,
        }
    return result


def full_library_profile(unit_rows, adapter_summary):
    """Model full PCR-library lengths without treating missing handles as absent.

    Assigned physical index pairs take precedence over sample index-length
    consensus. Missing index information is reported, never replaced by a fixed
    136 bp assumption. Both weighting modes use exactly these modeled units.
    """
    pairs = {p['scope_id']: p for p in adapter_summary.get('index_pairing', {}).get('pairs', [])}
    modeled = []
    excluded = collections.Counter()
    for row in unit_rows:
        family = row.get('family')
        if family not in FAMILIES:
            excluded['unknown_adapter_family'] += 1
            continue
        pair = pairs.get(row.get('scope_id')) if row.get('assignment_status') == 'assigned' else None
        if pair is not None and all(isinstance(pair.get(k), str) for k in ('i5', 'i7')):
            index_lengths = [len(pair[k]) for k in ('i5', 'i7')]
            basis = 'assigned_index_pair_consensus'
        else:
            index_lengths = [adapter_summary.get(k, {}).get('length') for k in ('i5', 'i7')]
            basis = 'sample_index_length_consensus'
        if any(isinstance(v, bool) or not isinstance(v, int) or v < 0 for v in index_lengths):
            excluded['unknown_index_length'] += 1
            continue
        insert = int(row['insert_length'])
        if insert < 0:
            raise ValueError('Full-library model requires non-negative insert lengths')
        added = len(P5) + len(P7) + sum(map(len, FAMILIES[family])) + sum(index_lengths)
        modeled.append(dict(unit_id=row['unit_id'], read_id=row['read_id'],
            insert_length_bp=insert, model_adapter_bases=added,
            model_library_length_bp=insert + added, model_basis=basis,
            outer_handles_detected=int(str(row.get('p5_detected')) == '1' and str(row.get('p7_detected')) == '1')))
    return dict(source_units=len(unit_rows), modeled_units=len(modeled),
                excluded_units=sum(excluded.values()), exclusion_reasons=dict(excluded),
                adapter_offsets_bp=sorted({r['model_adapter_bases'] for r in modeled}),
                model_basis_counts=dict(collections.Counter(r['model_basis'] for r in modeled)),
                outer_incomplete_modeled_units=sum(not r['outer_handles_detected'] for r in modeled),
                lengths=[r['model_library_length_bp'] for r in modeled], rows=modeled)


def write_library_size_assets(directory, prefix, label, unit_rows, adapter_summary,
                              maximum_bp=2000, bandwidth=.045):
    """Emit paired full-library plots plus their numeric profiles and model data."""
    import csv
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    model = full_library_profile(unit_rows, adapter_summary)
    names = dict(weighted=f'{prefix}.full_library_length_weighted.svg',
                 count=f'{prefix}.full_library_count_based.svg',
                 profiles=f'{prefix}.full_library_profiles.tsv',
                 units=f'{prefix}.full_library_units.tsv',
                 model=f'{prefix}.full_library_model.json')
    length_data = size_trace_data(model['lengths'], maximum_bp, bandwidth, weighting='length')
    count_data = size_trace_data(model['lengths'], maximum_bp, bandwidth, weighting='count')
    for name, weighting, caption in [('weighted', 'length', 'length-weighted'), ('count', 'count', 'count-based')]:
        svg = svg_size_trace(model['lengths'], f'{label}: ONT full-library profile ({caption})',
                             maximum_bp, bandwidth, weighting=weighting)
        (directory / names[name]).write_text(svg)
    with (directory / names['profiles']).open('w') as handle:
        writer = csv.writer(handle, delimiter='\t')
        writer.writerow(['size_bp', 'length_weighted_relative_signal', 'count_based_relative_signal'])
        writer.writerows(zip(length_data['lengths'], length_data['signal'], count_data['signal']))
    fields = ['unit_id', 'read_id', 'insert_length_bp', 'model_adapter_bases',
              'model_library_length_bp', 'model_basis', 'outer_handles_detected']
    with (directory / names['units']).open('w') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter='\t')
        writer.writeheader()
        writer.writerows(model['rows'])
    info = {k: v for k, v in model.items() if k not in ('rows', 'lengths')}
    info.update(length_definition='insert + P5 + family-specific inner adapters + P7 + physical i5/i7 lengths',
                coordinate='log1p(bp/50)', bandwidth=bandwidth, maximum_bp=maximum_bp,
                main_weighting='length', count_based_default_collapsed=True,
                normalization='each curve peak = 100 within the display window',
                measured_fluorescence=False, shown_units=count_data['shown'],
                above_range_units=count_data['above_range'],
                shown_model_bases=length_data['shown_weight'],
                above_range_model_bases=length_data['above_range_weight'])
    json_dump(info, directory / names['model'])
    names['model_summary'] = info
    return names


def read_index_assignments(path):
    assignments = {}
    with open_text(path) as handle:
        header = handle.readline().rstrip().split("\t")
        for line in handle:
            row = dict(zip(header, line.rstrip().split("\t")))
            assignments[row["unit_id"]] = row
    return assignments


def base_unit_id(qname):
    return qname[:-2] if qname.endswith(("/1", "/2")) else qname


def five_prime(record):
    """0-based coordinate of the outer aligned 5-prime base."""
    return record["start"] if record["strand"] == "+" else record["end"] - 1


def contact_call(cis, total, args):
    """Compatibility name: evidence count only, never a biological classification."""
    return contact_evidence_status(total, args.microc_minimum_contacts)


def microc_contacts(records, primary_sizes, args, unit_scopes=None):
    """Order-invariant endpoint screen; near inward pairs are a heuristic, not proven artifacts."""
    groups = collections.defaultdict(dict)
    for qname, rec in records.items():
        if qname.endswith(("/1", "/2")):
            groups[qname[:-2]][qname[-1]] = rec
    both = [(name, pair["1"], pair["2"]) for name, pair in groups.items() if "1" in pair and "2" in pair]
    mapq_pairs = [(name, a, b) for name, a, b in both
                  if all(args.mapq <= r["mapq"] < 255 and r["chrom"] in primary_sizes
                         and r.get("query_coverage", 1) >= args.microc_minimum_query_coverage for r in (a, b))]
    filtered = []
    same_fragment = 0
    for name, a, b in mapq_pairs:
        distance = abs(five_prime(a) - five_prime(b)) if a["chrom"] == b["chrom"] else None
        left, right = sorted((a, b), key=five_prime)
        if distance is not None and left["strand"] == "+" and right["strand"] == "-" and distance < args.microc_short_inward_bp:
            same_fragment += 1
            continue
        filtered.append((name, a, b, distance))
    seen, unique = set(), []
    for item in filtered:
        name, a, b, _ = item
        scope = (unit_scopes or {}).get(name) or "unassigned"
        # Exact strand-aware 5-prime coordinates avoid artificial 5-bp bin boundaries.
        key = (scope, tuple(sorted((r["chrom"], five_prime(r), r["strand"]) for r in (a, b))))
        if key not in seen:
            seen.add(key)
            unique.append(item)
    return both, mapq_pairs, filtered, unique, same_fragment


def contact_metrics(both, passed, filtered, unique, same_fragment, mappable, args):
    distances = [d for _, _, _, d in unique if d is not None]
    return dict(
        both_legs_mapped=len(both), both_legs_mapping_rate=len(both)/mappable if mappable else None,
        mapq_pass_pairs=len(passed), mapq_pass_pair_rate=len(passed)/len(both) if both else None,
        same_fragment_pairs=same_fragment, same_fragment_filter_pass_pairs=len(filtered),
        same_fragment_filter_pass_rate=len(filtered)/len(passed) if passed else None,
        unique_valid_contacts=len(unique), valid_contact_rate=len(unique)/mappable if mappable else None,
        contact_deduplication_retention=len(unique)/len(filtered) if filtered else None,
        **summarize_contacts(distances, len(unique), args.microc_long_range_cis_min_bp, args.microc_minimum_contacts),
        duplication_method=metric_profile()['dedup_method'],
        microc_minimum_contacts=args.microc_minimum_contacts)


def external_frip(records, path, primary_sizes, mapq):
    if not path:
        return {"external_frip_status": "NOT_CONFIGURED"}
    regions = collections.defaultdict(list)
    with open_text(path) as handle:
        for line in handle:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            chrom, start, end = line.split()[:3]
            start, end = int(start), int(end)
            if start < 0 or end <= start:
                raise ValueError(f"Invalid external peak interval: {line.strip()}")
            if chrom in primary_sizes:
                if end > primary_sizes[chrom]:
                    raise ValueError(f"External peak exceeds reference: {line.strip()}")
                regions[chrom].append((start, end))
    if not regions:
        raise ValueError("External peak BED has no intervals on configured primary chromosomes")
    regions = {chrom: merge_intervals(items) for chrom, items in regions.items()}
    usable = [r for r in records.values() if mapq <= r["mapq"] < 255 and r["chrom"] in primary_sizes]
    hits = sum(interval_overlaps(regions.get(r["chrom"]), r["start"], r["end"]) for r in usable)
    low, high = wilson(hits, len(usable))
    return dict(external_frip_status="DESCRIPTIVE" if usable else "INCONCLUSIVE_NO_USABLE_READS",
                external_peak_bed=str(path), 
                external_frip_reads=hits, external_frip_denominator=len(usable),
                external_frip=hits/len(usable) if usable else None,
                external_frip_ci95_low=low, external_frip_ci95_high=high)


def contributing_molecules(records):
    """Count parent ONT molecules, retaining the ordinary O/E evidence count."""
    parents = set()
    for record in records:
        name = record.get("qname", record.get("unit_id", ""))
        if name.endswith(("/1", "/2")):
            name = name[:-2]
        parent = record.get("read_id") or name.rsplit("|u", 1)[0]
        if parent:
            parents.add(parent)
    return len(parents)


def scope_qc(args, pair, profile, records, chrom_sizes, primary_names, mito, primary_sizes,
             sample_raw_reads):
    units = pair["support_units"]
    mappable = pair.get("mappable_units", units)
    assigned_raw_reads = pair.get("raw_reads", units)
    dimer_units = pair.get("adapter_dimer_units", 0)
    dimer_rate = pair.get("adapter_dimer_raw_reads", dimer_units) / assigned_raw_reads if assigned_raw_reads else None
    dimer_call = "inconclusive_low_units" if units < 20 else ("pass" if dimer_rate <= args.dimer_pass_max else ("warning" if dimer_rate <= args.dimer_warning_max else "fail"))
    result = {
        "scope_id": pair["scope_id"], "scope_name": pair.get("name", pair["scope_id"]),
        "i5": pair.get("i5"), "i7": pair.get("i7"), "index_pair_support_fraction": pair.get("support_fraction"),
        "raw_reads": assigned_raw_reads, "assigned_raw_reads": assigned_raw_reads,
        "complete_reads": assigned_raw_reads,
        "complete_read_rate": assigned_raw_reads / sample_raw_reads if sample_raw_reads else None,
        "complete_rate_denominator_raw_reads": sample_raw_reads,
        "library_units": units, "units_per_complete_read": units / assigned_raw_reads if assigned_raw_reads else None,
        "mappable_library_units": mappable,
        "adapter_dimer_units": dimer_units, "adapter_dimer_rate": dimer_rate,
        "adapter_dimer_call": dimer_call,
        "mean_read_q": profile["mean_read_q"],
        "raw_median_bp": profile["raw"]["median"],
        "raw_n50_bp": profile["raw"]["n50"],
        "insert_median_bp": profile["inserts"]["median"],
        "mapped_primary": len(records),
        "primary_mapping_rate": len(records) / (mappable * (2 if args.assay == "microc" else 1)) if mappable else None,
    }
    collision_records = [r for r in records.values() if args.mapq <= r["mapq"] < 255 and r["chrom"] in primary_sizes]
    result["metric_settings"] = metric_settings(args)
    if args.assay == "chromatin":
        result.update(collision_metric_fields(bin_stats(collision_records, primary_sizes, args.collision_bin_size)))
        result["collision_contributing_molecules"] = contributing_molecules(collision_records)
        mitochondrial = [r for r in records.values() if args.mapq <= r["mapq"] < 255 and r["chrom"] in mito]
        result.update(usable_fragments=len(collision_records),
                      usable_fragment_rate=len(collision_records) / mappable if mappable else None,
                      mitochondrial_fraction=len(mitochondrial) / (len(collision_records) + len(mitochondrial)) if collision_records or mitochondrial else None,
                      sequencing_yield_call="pass" if assigned_raw_reads >= args.minimum_raw and units >= args.minimum_units else "inconclusive_low_yield",
                      biological_qc_call="not_classified",
                      chromatin_qc_method="collision_oe_descriptive")
    else:
        result.update(contact_metrics(*microc_contacts(records, primary_sizes, args), mappable, args))
        result["sequencing_yield_call"] = "pass" if assigned_raw_reads >= args.minimum_raw and units >= args.minimum_units else "inconclusive_low_yield"
    result["adapter_dimer_unit_rate"] = dimer_units / units if units else None
    return result



def scope_structure_metrics(profile, assay, sample_units):
    """Count assigned units; molecule statistics deduplicate raw IDs within a scope."""
    rows = profile.get("unit_rows", [])
    units = len(rows)
    raw_ids = {r["read_id"] for r in rows}
    dimers = [r for r in rows if r["is_adapter_dimer"] == "1"]
    dimer_raw_ids = {r["read_id"] for r in dimers}
    outer = sum(r["p5_detected"] == "1" and r["p7_detected"] == "1" for r in rows)
    return dict(
        outer_complete_units=outer, outer_complete_unit_rate=outer/units if units else None,
        canonical_adapter_dimer_units=sum(r["dimer_class"] == "canonical_adapter_dimer" for r in rows),
        probable_adapter_dimer_units=sum(r["dimer_class"] == "probable_adapter_dimer" for r in rows),
        short_unmapped_units=sum(r["is_adapter_dimer"] != "1" and assay == "microc" and int(r["insert_length"]) < 50 for r in rows),
        raw_structure_dimer_rate=len(dimer_raw_ids)/len(raw_ids) if raw_ids else None,
        concatemer_reads=profile.get("concatemer_reads", 0),
        concatemer_rate=profile.get("concatemer_reads", 0)/len(raw_ids) if raw_ids else None,
        scope_unit_fraction=units/sample_units if sample_units else None,
        metric_unavailable_reasons={
            "index_assignment_rate": "Sample-level: assignment is required to define a sublibrary.",
            "index_ambiguous_fraction": "Not attributable: ambiguous units have no confident sublibrary identity.",
            "index_unassigned_fraction": "Not attributable: unassigned units have no sublibrary identity.",
            "unreconstructed_probable_dimer_reads": "Not attributable: these raw molecules have no reconstructed, assigned unit.",
        })


def mycoplasma_call(reads, loci, total, minimum):
    fraction = reads/total if total else 0
    if reads >= 10 and loci >= 3 or fraction >= .005 and reads >= 3 and loci >= 2:
        return "SIGNAL_DETECTED"
    if reads >= 3 and loci >= 2:
        return "SUSPECT"
    if reads:
        return "TRACE"
    return "NOT_DETECTED_AT_THIS_DEPTH" if total > 0 else "INCONCLUSIVE_LOW_SCREEN_DEPTH"


def scoped_mycoplasma(screen, allowed):
    """Partition competitive full-insert screening evidence by physical index scope."""
    if screen.get("status") == "NOT_CONFIGURED":
        return dict(mycoplasma_status="NOT_CONFIGURED", mycoplasma_screened_reads=None,
                    mycoplasma_qualifying_reads=None, mycoplasma_qualifying_loci=None,
                    mycoplasma_qualifying_fraction=None)
    evidence = screen.get("qualifying_read_loci")
    if evidence is None and screen.get("qualifying_reads") != 0:
        return dict(mycoplasma_status="REQUIRES_SCREEN_RERUN")
    evidence = evidence or {}
    hits = allowed & set(evidence)
    loci = {tuple(locus) for name in hits for locus in evidence[name]}
    total = len(allowed)
    return dict(mycoplasma_status=mycoplasma_call(len(hits), len(loci), total, screen.get("minimum_screened_reads", 500)),
                mycoplasma_screened_reads=total, mycoplasma_qualifying_reads=len(hits),
                mycoplasma_qualifying_loci=len(loci),
                mycoplasma_qualifying_fraction=len(hits)/total if total else None)

def cmd_metrics(args):
    with open(args.adapter_summary) as handle:
        adapter = json.load(handle)
    with open(args.mycoplasma) as handle:
        mycoplasma = json.load(handle)
    with open(args.spike_in) as handle:
        spike_in = json.load(handle)
    records = sam_primary_records(args.bam, args.samtools)
    chrom_sizes = read_chrom_sizes(args.chrom_sizes)
    primary_names, mito = resolve_chromosome_sets(chrom_sizes, args.primary, args.mito)
    primary_sizes = {name: chrom_sizes[name] for name in primary_names}
    if args.microc_minimum_contacts < 1 or args.microc_short_inward_bp < 0:
        raise ValueError("Micro-C requires positive minimum contacts and non-negative short inward cutoff")
    if not 0 <= args.microc_minimum_query_coverage <= 1:
        raise ValueError("Micro-C query coverage must be in [0, 1]")
    if not 0 <= args.dimer_pass_max <= args.dimer_warning_max <= 1:
        raise ValueError("Dimer thresholds must satisfy 0 <= pass <= warning <= 1")
    if args.microc_long_range_cis_min_bp < 0:
        raise ValueError("Micro-C long-range cis threshold must be non-negative")
    units = adapter["library_units"]
    mappable_units = adapter.get("mappable_library_units", units)
    adapter_dimer_units = adapter.get("adapter_dimer_units", 0)
    adapter_dimer_rate = adapter.get("adapter_dimer_rate", 0)
    if units < 20:
        adapter_dimer_call = "inconclusive_low_units"
    elif adapter_dimer_rate <= args.dimer_pass_max:
        adapter_dimer_call = "pass"
    elif adapter_dimer_rate <= args.dimer_warning_max:
        adapter_dimer_call = "warning"
    else:
        adapter_dimer_call = "fail"
    raw_count = adapter["raw"]["count"]
    common = {
        "sample": args.sample, "assay": args.assay, "reference_assembly": args.assembly, "raw_reads": raw_count,
        "complete_reads": adapter["complete_reads"], "complete_read_rate": adapter["complete_read_rate"],
        "library_units": units,
        "units_per_complete_read": units / adapter["complete_reads"] if adapter["complete_reads"] else None,
        "mappable_library_units": mappable_units,
        "adapter_dimer_max_bp": adapter.get("adapter_dimer_max_bp"),
        "outer_complete_units": adapter.get("outer_complete_units"),
        "outer_complete_unit_rate": adapter.get("outer_complete_unit_rate"),
        "short_unmapped_units": adapter.get("short_unmapped_units"),
        "unreconstructed_probable_dimer_reads": adapter.get("unreconstructed_probable_dimer_reads"),
        "index_assignment_rate": adapter.get("index_pairing", {}).get("all_unit_assignment_rate"),
        "index_ambiguous_fraction": adapter.get("index_pairing", {}).get("ambiguous_fraction"),
        "index_unassigned_fraction": adapter.get("index_pairing", {}).get("unassigned_fraction"),
        "adapter_dimer_units": adapter_dimer_units, "adapter_dimer_rate": adapter_dimer_rate,
        "adapter_dimer_unit_rate": adapter.get("adapter_dimer_unit_rate"),
        "canonical_adapter_dimer_units": adapter.get("canonical_adapter_dimer_units", 0),
        "probable_adapter_dimer_units": adapter.get("probable_adapter_dimer_units", 0),
        "raw_structure_dimer_rate": adapter.get("raw_structure_dimer_rate"),
        "adapter_dimer_raw_window_bp": adapter.get("adapter_dimer_raw_window_bp"),
        "adapter_dimer_call": adapter_dimer_call, "concatemer_rate": adapter["concatemer_rate"],
        "raw_median_bp": adapter["raw"]["median"], "raw_n50_bp": adapter["raw"]["n50"],
        "insert_median_bp": adapter["inserts"]["median"], "mean_read_q": adapter["mean_read_q"],
        "adapter_family": max(adapter["adapter_families"], key=adapter["adapter_families"].get)
                          if adapter["adapter_families"] else "not_detected",
        "index_architecture": adapter["index_architecture"], "index_pairing": adapter.get("index_pairing", {}),
        "i5": adapter["i5"]["sequence"], "i5_purity": adapter["i5"]["purity"],
        "i7": adapter["i7"]["sequence"], "i7_purity": adapter["i7"]["purity"],
        "aligner": args.aligner, "mycoplasma_screen": {k: v for k, v in mycoplasma.items() if k != "qualifying_read_loci"}, "mycoplasma_status": mycoplasma.get("status"),
        "mycoplasma_screened_reads": mycoplasma.get("screened_reads"),
        "mycoplasma_qualifying_reads": mycoplasma.get("qualifying_reads"),
        "mycoplasma_qualifying_loci": mycoplasma.get("qualifying_loci"),
        "mycoplasma_qualifying_fraction": mycoplasma.get("qualifying_fraction"),
        "spike_in_screen": {key: value for key, value in spike_in.items() if key not in {"qualifying_read_ids", "screened_read_ids", "candidate_read_species", "qualifying_read_loci"}},
        "spike_in_label": spike_in.get("label"),
        "spike_in_status": spike_in.get("status"),
        "spike_in_screened_reads": spike_in.get("screened_reads"),
        "spike_in_reads": spike_in.get("qualifying_reads"),
        "spike_in_fraction": spike_in.get("qualifying_fraction"),
        "primary_chromosomes_used": primary_names,
        "mitochondrial_chromosomes_used": sorted(mito),
        "mapped_primary": len(records),
        "primary_mapping_rate": len(records) / (mappable_units * (2 if args.assay == "microc" else 1)) if mappable_units else None,
    }
    Path(args.fragments).parent.mkdir(parents=True, exist_ok=True)
    Path(args.pairs).parent.mkdir(parents=True, exist_ok=True)

    collision_records = [r for r in records.values() if args.mapq <= r["mapq"] < 255 and r["chrom"] in primary_sizes]
    common["metric_settings"] = metric_settings(args)
    if args.assay == "chromatin":
        collision = collision_metric_fields(bin_stats(collision_records, primary_sizes, args.collision_bin_size))
        collision["collision_contributing_molecules"] = contributing_molecules(collision_records)
        usable = collision_records
        mitochondrial = [r for r in records.values() if args.mapq <= r["mapq"] < 255 and r["chrom"] in mito]
        metrics = dict(common, **collision, usable_fragments=len(usable),
                       usable_fragment_rate=len(usable) / mappable_units if mappable_units else None,
                       mitochondrial_fraction=len(mitochondrial) / (len(usable) + len(mitochondrial)) if usable or mitochondrial else None,
                       sequencing_yield_call="pass" if raw_count >= args.minimum_raw and units >= args.minimum_units else "inconclusive_low_yield",
                       biological_qc_call="not_classified",
                       chromatin_qc_method="collision_oe_descriptive",
                       classification_note="Chromatin enrichment is not classified. Collision O/E is the depth-normalized concentration diagnostic (random baseline 1); the CPM-normalized TSS-to-TES profile is visual only.")
        with gzip.open(args.fragments, "wt") as out:
            for r in sorted(usable, key=lambda x: (x["chrom"], x["start"], x["end"])):
                out.write(f"{r['chrom']}\t{r['start']}\t{r['end']}\t{r['qname']}\t{r['mapq']}\t{r['strand']}\n")
        with gzip.open(args.pairs, "wt") as out:
            out.write("# chromatin sample; no contact pairs\n")
    else:
        assignment_scopes = {u: row.get("scope_id", "") for u, row in read_index_assignments(args.index_assignments).items()}
        both, passed, filtered, unique, same_fragment = microc_contacts(records, primary_sizes, args, assignment_scopes)
        metrics = dict(common, **contact_metrics(both, passed, filtered, unique, same_fragment, mappable_units, args),
                       sequencing_yield_call="pass" if raw_count >= args.minimum_raw and units >= args.minimum_units else "inconclusive_low_yield",
                       classification_note=f"Micro-C endpoint composition is not assigned a biological success/failure label. The {args.microc_minimum_contacts}-contact minimum describes evidence count only. Cis and long-range cis use all unique valid contacts; Wilson 95% intervals assume independent contacts and do not include molecule dependence or recovery bias. Short inward pairs below {args.microc_short_inward_bp} bp are excluded heuristically and may include real local contacts. Full-insert ligation structure is reported separately.")
        with gzip.open(args.fragments, "wt") as out:
            out.write("# Micro-C sample; see pairs file\n")
        with gzip.open(args.pairs, "wt") as out:
            out.write("## pairs format v1.0\n#shape: whole matrix\n")
            out.write(f"#genome_assembly: {args.assembly}\n")
            for chrom, size in chrom_sizes.items():
                out.write(f"#chromsize: {chrom} {size}\n")
            out.write("#columns: readID chrom1 pos1 chrom2 pos2 strand1 strand2 mapq1 mapq2\n")
            for name, a, b, _ in unique:
                out.write(f"{name}\t{a['chrom']}\t{five_prime(a)+1}\t{b['chrom']}\t{five_prime(b)+1}\t{a['strand']}\t{b['strand']}\t{a['mapq']}\t{b['mapq']}\n")
    if metrics.get("mycoplasma_screened_reads") and metrics.get("mycoplasma_qualifying_fraction") == 0:
        metrics["mycoplasma_status"]="NOT_DETECTED_AT_THIS_DEPTH"
    if spike_in.get("mode") == "auto":
        metrics.update(auto_spikein_summary(spike_in))
    ligation_path=getattr(args, "ligation", "")
    ligation=json.loads(Path(ligation_path).read_text()) if ligation_path else None
    if ligation is not None:
        metrics.update({key:value for key,value in ligation.items() if key.startswith("ligation_")})
        metrics["ligation_method"]=ligation["method"]
        metrics["ligation_evidence_status"]=ligation["status"]
    metrics["metric_schema_version"] = 3
    pairing = adapter.get("index_pairing", {})
    metrics["sub_libraries"] = []
    if pairing.get("mode") == "multiplexed_library":
        profiles = read_scope_profiles(args.read_table)
        assignment = read_index_assignments(args.index_assignments)
        for pair in pairing.get("pairs", []):
            allowed = {unit for unit, row in assignment.items() if row.get("scope_id") == pair["scope_id"] and row.get("assignment_status") == "assigned" and row.get("is_adapter_dimer") != "1"}
            scoped_records = {qname: rec for qname, rec in records.items() if base_unit_id(qname) in allowed}
            scoped = scope_qc(args, pair, profiles.get(pair["scope_id"], {"mean_read_q": None, "raw": length_summary([]), "inserts": length_summary([])}), scoped_records, chrom_sizes,
                              primary_names, mito, primary_sizes, raw_count)
            scoped.update(scope_structure_metrics(profiles.get(pair["scope_id"], {}), args.assay, units))
            scoped.update(scoped_mycoplasma(mycoplasma, allowed))
            spike_hits = set(spike_in.get("qualifying_read_ids", []))
            scoped_spike_reads = len(allowed & spike_hits)
            if spike_in.get("status") == "NOT_CONFIGURED":
                scoped_spike_status = "NOT_CONFIGURED"
            elif scoped_spike_reads:
                scoped_spike_status = "DETECTED"
            elif len(allowed) < spike_in.get("minimum_screened_reads", 500):
                scoped_spike_status = "INCONCLUSIVE_LOW_SCREEN_DEPTH"
            else:
                scoped_spike_status = "NOT_DETECTED"
            scoped.update(
                spike_in_label=spike_in.get("label"),
                spike_in_status=scoped_spike_status,
                spike_in_screened_reads=len(allowed) if scoped_spike_status != "NOT_CONFIGURED" else None,
                spike_in_reads=scoped_spike_reads if scoped_spike_status != "NOT_CONFIGURED" else None,
                spike_in_fraction=scoped_spike_reads / len(allowed) if allowed and scoped_spike_status != "NOT_CONFIGURED" else None,
            )
            if spike_in.get("mode") == "auto":
                scoped.update(auto_spikein_summary(spike_in, allowed))
            if ligation is not None:
                from microc_ligation import summarize
                scope_ligation=ligation.get("scopes", {}).get(pair["scope_id"], summarize([]))
                scoped.update({key:value for key,value in scope_ligation.items() if key.startswith("ligation_")})
            metrics["sub_libraries"].append(scoped)
    json_dump(metrics, args.metrics)
    print(json.dumps(metrics, indent=2))


def svg_profile_placeholder(title, message="No usable mapped reads"):
    return (
        "<svg xmlns=\"http://www.w3.org/2000/svg\" viewBox=\"0 0 820 380\" style=\"font-family:Arial,Helvetica,sans-serif\">"
        "<rect width=\"100%\" height=\"100%\" fill=\"white\"/>"
        f"<text x=\"410\" y=\"32\" text-anchor=\"middle\" font-size=\"17\">{html.escape(title)}</text>"
        f"<text x=\"410\" y=\"190\" text-anchor=\"middle\" fill=\"#687986\">{html.escape(message)}</text></svg>"
    )



def write_sparse_cpm_bedgraph(samtools, bam, mapq, chrom_sizes, primary_names, path):
    """Write sparse CPM coverage from mapped alignment spans without scanning the whole genome."""
    events=collections.defaultdict(list)
    mapped=0
    process=subprocess.Popen([samtools, "view", "-q", str(mapq), "-F", "3844", str(bam)],
                             text=True, stdout=subprocess.PIPE)
    for line in process.stdout:
        fields=line.rstrip().split("\t")
        if len(fields) < 6:
            continue
        chrom=fields[2]
        if chrom not in primary_names or chrom not in chrom_sizes:
            continue
        if int(fields[4]) == 255:
            continue
        mapped += 1
        import re
        start=max(0, int(fields[3])-1)
        for number, op in re.findall(r"(\d+)([MIDNSHP=X])", fields[5]):
            length = int(number)
            if op in "M=X":
                end=min(chrom_sizes[chrom], start+length)
                if end > start:
                    events[chrom].extend(((start, 1), (end, -1)))
                start += length
            elif op in "DN":
                start += length
    if process.wait() != 0:
        raise RuntimeError(f"samtools view failed while building sparse coverage for {bam}")
    if mapped == 0 or not events:
        path.write_text("")
        return mapped, False
    scale=1e6 / mapped
    covered=False
    with path.open("w") as handle:
        for chrom in sorted(primary_names):
            if chrom not in events:
                continue
            depth=0
            previous=None
            for position, delta in sorted(events[chrom], key=lambda item: item[0]):
                if previous is not None and position > previous and depth > 0:
                    handle.write(f"{chrom}\t{previous}\t{position}\t{depth*scale:.8g}\n")
                    covered=True
                depth += delta
                previous=position
    return mapped, covered

def cmd_track(args):
    sizes = read_chrom_sizes(args.chrom_sizes)
    primary, _ = resolve_chromosome_sets(sizes, args.primary, "auto")
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="ont_track_") as directory:
        bedgraph = Path(directory) / "coverage.bedGraph"
        _, covered = write_sparse_cpm_bedgraph(args.samtools, args.bam, args.mapq, sizes, primary, bedgraph)
        if not covered:
            bedgraph.write_text("".join(f"{chrom}\t0\t{sizes[chrom]}\t0\n" for chrom in sorted(primary)))
        subprocess.run([args.converter, str(bedgraph), args.chrom_sizes, args.output], check=True)


def build_sublibrary_gene_profile(out, sample, scope_id, label, bam, genes, qnames,
                                  samtools, bedgraph_to_bigwig, compute_matrix, plot_profile,
                                  mapq, chrom_sizes_path, chrom_sizes, primary_names,
                                  upstream, downstream, body, bin_size, threads):
    """Build a CPM-normalized TSS–TES profile from one assigned sublibrary."""
    profile_dir = out / "qc" / "gene_profile"
    plot_dir = out / "plots" / "gene_profile"
    asset_dir = out / "reports" / "assets"
    profile_dir.mkdir(parents=True, exist_ok=True)
    plot_dir.mkdir(parents=True, exist_ok=True)
    asset_dir.mkdir(parents=True, exist_ok=True)
    plot_path = plot_dir / f"{sample}.{scope_id}.TSS_TES.svg"
    data_path = profile_dir / f"{sample}.{scope_id}.profile.tsv"
    matrix_path = profile_dir / f"{sample}.{scope_id}.matrix.gz"
    asset_path = asset_dir / f"{sample}.{scope_id}.TSS_TES.svg"
    if not qnames:
        svg = svg_profile_placeholder(f"{label}: CPM-normalized TSS–TES profile", "No assigned, non-dimer reads")
        plot_path.write_text(svg)
        asset_path.write_text(svg)
        data_path.write_text("No assigned, non-dimer reads\n")
        return asset_path
    with tempfile.TemporaryDirectory(prefix=f"ont_{sample}_{scope_id}_", dir=str(profile_dir)) as tmp_name:
        tmp = Path(tmp_name)
        qname_path = tmp / "qnames.txt"
        qname_path.write_text("\n".join(sorted(qnames)) + "\n")
        scoped_bam = tmp / "scoped.bam"
        with scoped_bam.open("wb") as handle:
            subprocess.run([samtools, "view", "-b", "-N", str(qname_path), str(bam)], stdout=handle, check=True)
        sparse_bedgraph = tmp / "scoped.CPM.bedGraph"
        mapped, covered = write_sparse_cpm_bedgraph(samtools, scoped_bam, mapq, chrom_sizes, primary_names, sparse_bedgraph)
        if mapped == 0 or not covered:
            svg = svg_profile_placeholder(f"{label}: CPM-normalized TSS–TES profile", "No MAPQ-pass primary mapped reads")
            plot_path.write_text(svg)
            asset_path.write_text(svg)
            data_path.write_text("No MAPQ-pass primary mapped reads\n")
            return asset_path
        scoped_bw = tmp / "scoped.CPM.bw"
        subprocess.run([bedgraph_to_bigwig, str(sparse_bedgraph), str(chrom_sizes_path), str(scoped_bw)], check=True)
        subprocess.run([compute_matrix, "scale-regions", "-S", str(scoped_bw), "-R", str(genes),
                        "--beforeRegionStartLength", str(upstream), "--afterRegionStartLength", str(downstream),
                        "--regionBodyLength", str(body), "--binSize", str(bin_size), "--missingDataAsZero",
                        "--sortRegions", "keep", "--samplesLabel", label, "--startLabel", "TSS", "--endLabel", "TES",
                        "-p", str(threads), "-o", str(matrix_path)], check=True)
        subprocess.run([plot_profile, "-m", str(matrix_path), "-o", str(plot_path),
                        "--outFileNameData", str(data_path), "--plotType", "se",
                        "--startLabel", "TSS", "--endLabel", "TES",
                        "--yAxisLabel", "CPM-normalized signal", "--samplesLabel", label], check=True)
        asset_path.write_bytes(plot_path.read_bytes())
    return asset_path

def read_numbers(path):
    with open(path) as handle:
        return [float(x) for x in handle if x.strip()]


def svg_hist(values, title, x_label, log_x=False):
    width, height, left, top, plot_w, plot_h = 760, 380, 70, 65, 650, 245
    if not values:
        return svg_profile_placeholder(title, "No reconstructed observations")
    # Never pile the >99th percentile tail into the final visible bin.
    upper = quantile(values, .99)
    shown = [v for v in values if v <= upper]
    transformed = [math.log10(1 + v) if log_x else v for v in shown]
    lo, hi = min(transformed), max(transformed)
    if hi <= lo:
        lo, hi = max(0, lo - .05), hi + .05
    bins = 45
    counts = [0] * bins
    for value in transformed:
        index = min(bins - 1, int((value - lo) / (hi - lo) * bins))
        counts[index] += 1
    fractions = [n / len(values) for n in counts]
    peak = max(fractions)
    elements = [f"<svg xmlns='http://www.w3.org/2000/svg' viewBox='0 0 {width} {height}' style='font-family:Arial,Helvetica,sans-serif'>",
                "<rect width='100%' height='100%' fill='white'/>",
                f"<text x='{width/2}' y='25' text-anchor='middle' font-size='17'>{html.escape(title)}</text>",
                f"<text x='{width/2}' y='46' text-anchor='middle' font-size='11'>n={len(values):,}; zero bp={sum(v == 0 for v in values):,}; above P99={len(values)-len(shown):,} (omitted, not folded into last bin)</text>"]
    for index, fraction in enumerate(fractions):
        x = left + index / bins * plot_w
        h = fraction / peak * plot_h
        elements.append(f"<rect x='{x:.2f}' y='{top+plot_h-h:.2f}' width='{plot_w/bins-.5:.2f}' height='{h:.2f}' fill='#4C78A8'><title>{counts[index]} observations; {fraction:.2%} of all observations</title></rect>")
    for index in range(5):
        value = lo + index * (hi-lo) / 4
        label = 10 ** value - 1 if log_x else value
        x = left + index * plot_w / 4
        elements.append(f"<text x='{x:.1f}' y='{top+plot_h+23}' text-anchor='middle' font-size='11'>{label:.0f}</text>")
    for index in range(4):
        fraction = index * peak / 3
        y = top + plot_h - index * plot_h / 3
        elements.append(f"<text x='{left-8}' y='{y+4:.1f}' text-anchor='end' font-size='11'>{fraction:.1%}</text>")
    scale = "log10(1 + bp) bins" if log_x else "linear bins"
    elements.append(f"<text x='{left+plot_w/2}' y='{height-12}' text-anchor='middle' font-size='12'>Insert/read length (bp); {scale}; fraction of observations per bin</text></svg>")
    return "".join(elements)


def svg_bars(labels, values, title):
    width, height, left, top, plot_w, plot_h = 760, 360, 70, 42, 650, 250
    peak = max(values) if values else 1
    peak = peak or 1
    bar_w = plot_w / max(1, len(values)) * .65
    bars = []
    for i, (label, value) in enumerate(zip(labels, values)):
        x = left + (i + .5) * plot_w / len(values) - bar_w / 2
        h = value / peak * plot_h
        bars.append(f"<rect x='{x:.1f}' y='{top+plot_h-h:.1f}' width='{bar_w:.1f}' height='{h:.1f}' fill='#4C78A8'/><text x='{x+bar_w/2:.1f}' y='{top+plot_h-h-5:.1f}' text-anchor='middle' font-size='11'>{value:.3g}</text><text x='{x+bar_w/2:.1f}' y='{top+plot_h+20}' text-anchor='middle' font-size='11'>{html.escape(label)}</text>")
    return f"<svg xmlns='http://www.w3.org/2000/svg' viewBox='0 0 {width} {height}' style='font-family:Arial,Helvetica,sans-serif'><rect width='100%' height='100%' fill='white'/><text x='{width/2}' y='24' text-anchor='middle' font-size='18'>{html.escape(title)}</text><line x1='{left}' y1='{top+plot_h}' x2='{left+plot_w}' y2='{top+plot_h}' stroke='#333'/>{''.join(bars)}</svg>"


def call_color(call):
    call = call or ""
    if "pass" in call or "high_quality" in call:
        return "#16835D"
    if "low_quality" in call or "borderline" in call or "inconclusive" in call or "warning" in call:
        return "#C77A18"
    return "#C23B3B"


METRIC_META = {
    "collision_contributing_molecules": ("Contributing ONT molecules", "Distinct original ONT molecules contributing eligible chromatin fragments; units from one molecule are dependent."),
    "biological_qc_call": ("Cis-ratio evidence", "Confidence-aware endpoint cis screen only; cannot establish complete Micro-C library success."),
    "outer_complete_unit_rate": ("Both outer handles detected", "Reconstructed units with both outer P5 and P7 handles / all reconstructed units."),
    "outer_complete_units": ("Units with both outer handles", "Both outer P5/P7 handles detected, beyond the inner anchors required for reconstruction."),
    "index_assignment_rate": ("All-unit index assignment", "Assigned reconstructed units / all reconstructed units, including units missing index information."),
    "index_ambiguous_fraction": ("Ambiguous index fraction", "Units tied between index pairs / all reconstructed units."),
    "index_unassigned_fraction": ("Unassigned index fraction", "Units without an assignment / all reconstructed units."),
    "unreconstructed_probable_dimer_reads": ("Unresolved dimer candidates", "Unreconstructed raw reads matching a zero-insert template; excluded from confirmed burden because true short inserts can match."),
    "short_unmapped_units": ("Too short for two legs", "Non-dimer Micro-C inserts below 50 bp, excluded from mapping denominators."),
    "raw_reads": ("Raw reads", "ONT molecules supplied to the pipeline before adapter interpretation."),
    "sequencing_yield_call": ("ONT sequencing QC", "Technical sequencing-yield call based on the configured minimum raw-read and reconstructed-unit requirements."),
    "mean_read_q": ("Mean read Q", "Mean across raw molecules of -10 log10(mean base error probability), using all bases; differs from arithmetic mean base Q and basecaller-specific trimming."),
    "raw_median_bp": ("Raw median length", "Median full ONT molecule length before splitting."),
    "raw_n50_bp": ("Raw read N50", "Length above which half of all sequenced bases occur."),
    "insert_median_bp": ("Insert median length", "Median length of adapter-bounded NGS library inserts reconstructed from ONT molecules."),
    "complete_read_rate": ("Inner-anchor reconstruction rate", "Total: raw molecules containing both inner adapters; outer P5/P7 handles are not required. Sublibrary: contribution to all sample raw molecules, not a reconstruction success rate."),
    "library_units": ("Library units", "Number of adapter-bounded library units recovered from ONT molecules."),
    "units_per_complete_read": ("Units per complete read", "Number of inner-adapter-bounded units per reconstructed ONT molecule."),
    "mappable_library_units": ("Mappable units", "Non-dimer query units actually emitted for mapping; Micro-C also requires an insert of at least 50 bp (two nonoverlapping legs of at least 25 bp)."),
    "adapter_dimer_max_bp": ("Adapter-dimer cutoff", "Maximum adapter-bounded insert length counted as adapter dimer."),
    "adapter_dimer_units": ("Adapter-dimer units", "Units called by a canonical central gap or pair-specific zero-insert full-molecule template."),
    "canonical_adapter_dimer_units": ("Canonical adapter dimers", "Units with a central insert gap at or below the configured canonical cutoff."),
    "probable_adapter_dimer_units": ("Structure-matched adapter dimers", "Additional units matching their i5–i7-specific zero-insert full-molecule template."),
    "raw_structure_dimer_rate": ("Raw structure-dimer rate", "Fraction of raw ONT molecules supported as zero-insert dimers by reconstructed library structure."),
    "adapter_dimer_rate": ("Adapter-dimer burden", "Total: raw molecules with a reconstructed dimer-like unit / all raw molecules. Sublibrary: assigned raw molecules with a dimer-like unit / assigned raw molecules. Unreconstructed template candidates are excluded."),
    "adapter_dimer_unit_rate": ("Reconstructed-unit dimer rate", "Fraction of reconstructed library units called as zero-insert dimers; raw-molecule burden is used for the traffic-light call."),
    "adapter_dimer_call": ("Adapter-dimer QC", "Screening call based on adapter-dimer burden and available unit count."),
    "concatemer_rate": ("Concatemer rate", "Total: fraction of all raw molecules with multiple reconstructed units. Sublibrary: fraction of its assigned raw molecules with multiple units, including units from other scopes; shared molecules may contribute to multiple scopes."),
    "scope_unit_fraction": ("Sublibrary share of all units", "Assigned units in this sublibrary / all reconstructed sample units, including unassigned and ambiguous units."),
    "adapter_family": ("Adapter family", "Best-supported Illumina-compatible P5/P7 sequence family."),
    "index_architecture": ("Index architecture", "Detected single-index, dual-index, mixed or unresolved architecture."),
    "i5": ("i5 consensus", "Consensus i5 index sequence normalized to P5-to-P7 orientation."),
    "i5_purity": ("i5 purity", "Fraction supporting the dominant i5 consensus among informative units."),
    "i7": ("i7 consensus", "Consensus i7 index sequence normalized to P5-to-P7 orientation."),
    "i7_purity": ("i7 purity", "Fraction supporting the dominant i7 consensus among informative units."),
    "aligner": ("Aligner", "Reference aligner selected for this assay."),
    "mycoplasma_status": ("Mycoplasma signal screen", "Opportunistic competitive-alignment signal screen; NOT_DETECTED_AT_THIS_DEPTH does not exclude contamination."),
    "mycoplasma_screened_reads": ("Mycoplasma-screened inserts", "Adapter-trimmed full inserts competitively aligned against host, eight common cell-culture Mollicutes and bacterial decoys."),
    "mycoplasma_qualifying_reads": ("Mycoplasma-supporting reads", "Unique inserts passing target-contig, MAPQ, aligned-length, coverage and identity filters."),
    "mycoplasma_qualifying_loci": ("Mycoplasma-supporting loci", "Distinct target genome 100-bp windows supported by qualifying inserts."),
    "mycoplasma_qualifying_fraction": ("Mycoplasma signal fraction", "Qualifying target reads divided by all screened adapter-trimmed inserts; a zero shallow result is not proof of absence."),
    "spike_in_status": ("Spike-in candidate screen", "Automatic competitive reference screen. Foreign-sequence evidence does not prove intentional spike-in; none was added to the current libraries."),
    "spike_in_label": ("Candidate species", "Species with qualifying foreign-sequence assignments; same-species spike-in cannot be separated from host without specific barcodes."),
    "ligation_screened_units": ("Ligation-screened units", "All adapter-bounded non-dimer units, including unresolved cases. This is not the number of independent raw ONT molecules."),
    "ligation_resolved_units": ("Units with resolved block structure", "Query coverage >=80%, minimum aligned block 40 bp, MAPQ >=20 and <255, identity >=85%, and no unexplained query gap >30 bp."),
    "ligation_resolved_fraction": ("Resolved structure fraction", "Resolved units / all ligation-screened units."),
    "ligation_single_block_units": ("One continuous genomic block", "No resolved ligation junction. Adjacent or collinear ligated fragments can be indistinguishable from one block; this is not proof of unligated DNA."),
    "ligation_two_units": ("2-fragment units", "Units with two resolved genomic blocks; evidence-supported minimum fragment count."),
    "ligation_three_units": ("3-fragment units", "Units with three resolved genomic blocks; evidence-supported minimum fragment count."),
    "ligation_four_units": ("4-fragment units", "Units with four resolved genomic blocks; evidence-supported minimum fragment count."),
    "ligation_five_plus_units": ("≥5-fragment units", "Units with at least five resolved genomic blocks; short or repetitive fragments can remain undetected."),
    "ligation_unresolved_units": ("Unresolved ligation structure", "Insufficient unique mapping or query coverage, unexplained gaps, or large indels; kept in the denominator."),
    "ligation_multifragment_fraction": ("Multi-fragment evidence / all units", "Units with >=2 resolved genomic blocks / all screened library units, including unresolved cases in the denominator."),
    "ligation_multifragment_fraction_resolved": ("Multi-fragment evidence / resolved units", "Units with >=2 resolved blocks / units with resolved structure. Conditional on successful structure resolution."),
    "spike_in_screened_reads": ("Spike-in-screened inserts", "Adapter-trimmed inserts included in the host plus spike-in competitive alignment."),
    "spike_in_reads": ("Spike-in-assigned reads", "Unique inserts assigned to target spike-in contigs after MAPQ, alignment-length, coverage and identity filters."),
    "spike_in_fraction": ("Spike-in fraction", "Spike-in-assigned inserts divided by all inserts screened in the competitive reference."),
    "primary_mapping_rate": ("Primary-alignment rate", "Any mapped primary alignment / emitted query legs; primary alignment is distinct from configured primary chromosomes. MAPQ filtering is applied separately."),
    "usable_fragments": ("Usable fragments", "MAPQ-passing chromatin fragments on primary chromosomes."),
    "usable_fragment_rate": ("Usable fragment rate", "Usable chromatin fragments divided by mappable library units."),
    "mitochondrial_fraction": ("Mitochondrial fraction", "Mitochondrial / (configured primary nuclear + mitochondrial) MAPQ-pass alignments; excludes other contigs."),
    "collision_bin_size": ("Collision bin size", "Genomic bin size used for the primary collision statistic."),
    "collision_oe": ("Collision O/E (depth normalized)", "Observed same-bin pairs / expected pairs under uniform sampling of configured nuclear bases. Expected value is depth independent; 1 is random background. Low-depth estimates remain noisy; reference, bin size and mappability still matter. Not target-specific enrichment."),
    "collision_fragments": ("Collision observations", "MAPQ-pass nuclear chromatin fragments used for concentration."),
    "collision_observed_pairs": ("Observed same-bin pairs", "Sum over bins of n_j choose 2. Rare/zero collision counts at shallow depth may be inconclusive."),
    "collision_null_expected_pairs": ("Expected same-bin pairs", "N choose 2 times the uniform-base pair probability, corrected for chromosome terminal-bin widths."),
    "both_legs_mapped": ("Both legs mapped", "Micro-C units with a primary alignment for both inferred ligation legs."),
    "both_legs_mapping_rate": ("Both-leg mapping rate", "Fraction of mappable Micro-C units with both legs mapped."),
    "mapq_pass_pairs": ("MAPQ-pass pairs", "Both-leg Micro-C pairs passing MAPQ (excluding 255), query-coverage and primary-chromosome filters."),
    "mapq_pass_pair_rate": ("MAPQ-pass pair rate", "MAPQ-pass pairs divided by pairs with both legs mapped."),
    "same_fragment_pairs": ("Short inward pairs excluded", "Heuristic based on strand-aware aligned 5-prime separation; can exclude genuine local contacts and does not prove a shared fragment."),
    "same_fragment_filter_pass_pairs": ("Post-artifact-filter contacts", "MAPQ-pass pairs retained after the short inward-facing same-fragment filter."),
    "same_fragment_filter_pass_rate": ("Artifact-filter retention", "Post-artifact-filter contacts divided by MAPQ-pass pairs."),
    "unique_valid_contacts": ("Unique endpoint contacts", "Screened endpoint pairs after exact 5-prime coordinate deduplication within physical index scope; not junction-validated ligation events."),
    "valid_contact_rate": ("Valid-contact rate", "Unique valid contacts divided by mappable reconstructed library units."),
    "contact_deduplication_retention": ("Contact deduplication retention", "Unique valid contacts divided by post-artifact-filter contacts."),
    "duplicate_contact_fraction": ("Duplicate contact fraction", "Fraction of filtered Micro-C contacts removed as coordinate duplicates."),
    "cis_contacts": ("Cis contacts", "Unique valid Micro-C contacts whose legs map to the same chromosome."),
    "trans_contacts": ("Trans contacts", "Unique valid Micro-C contacts whose legs map to different chromosomes."),
    "cis_ratio": ("Cis ratio", "Fraction of unique valid Micro-C contacts mapping within the same chromosome."),
    "cis_ratio_ci95_low": ("Cis ratio CI low", "Lower bound of the Wilson 95% confidence interval for cis ratio."),
    "cis_ratio_ci95_high": ("Cis ratio CI high", "Upper bound of the Wilson 95% confidence interval for cis ratio."),
    "long_range_cis_ci95_low": ("Long-range cis CI low", "Lower Wilson 95% bound for long-range cis contacts divided by all unique contacts."),
    "long_range_cis_ci95_high": ("Long-range cis CI high", "Upper Wilson 95% bound for long-range cis contacts divided by all unique contacts."),
    "cis_lt1kb": ("Cis contacts <1 kb", "Unique cis endpoint contacts separated by less than 1 kb."),
    "cis_1_10kb": ("Cis contacts 1–10 kb", "Unique cis endpoint contacts separated by at least 1 kb and less than 10 kb."),
    "cis_10_100kb": ("Cis contacts 10–100 kb", "Unique cis endpoint contacts separated by at least 10 kb and less than 100 kb."),
    "cis_ge100kb": ("Cis contacts ≥100 kb", "Unique cis endpoint contacts separated by at least 100 kb."),
    "long_range_cis_min_bp": ("Long-range cis threshold", "Minimum genomic separation used to define a long-range cis contact."),
    "long_range_cis_contacts": ("Long-range cis contacts", "Unique valid cis contacts at or above the configured genomic-separation threshold."),
    "long_range_cis_fraction": ("Long-range cis fraction", "Long-range cis contacts divided by all unique valid contacts."),
    "long_range_within_cis_fraction": ("Long-range fraction within cis", "Long-range cis contacts divided only by unique valid cis contacts."),
}

PERCENT_METRICS = {
    "complete_read_rate", "adapter_dimer_rate", "raw_structure_dimer_rate", "concatemer_rate", "i5_purity", "i7_purity",
    "primary_mapping_rate",
    "usable_fragment_rate", "mitochondrial_fraction",
"both_legs_mapping_rate", "duplicate_contact_fraction",
    "cis_ratio", "cis_ratio_ci95_low", "cis_ratio_ci95_high", "long_range_cis_ci95_low", "long_range_cis_ci95_high", "mycoplasma_qualifying_fraction",
    "spike_in_fraction", "ligation_resolved_fraction", "ligation_multifragment_fraction", "ligation_multifragment_fraction_resolved", "mapq_pass_pair_rate", "same_fragment_filter_pass_rate",
    "valid_contact_rate", "contact_deduplication_retention",
    "long_range_cis_fraction", "long_range_within_cis_fraction",
    "adapter_dimer_unit_rate", "outer_complete_unit_rate", "index_assignment_rate", "index_ambiguous_fraction", "index_unassigned_fraction", "scope_unit_fraction",
}

DEFAULT_REFERENCE_RANGES = {
    "raw_reads": {"min": 500, "warning_min": 250},
    "mean_read_q": None,
    "complete_read_rate": {"min": .50, "max": 1.0, "warning_min": .30, "type": "percent"},
    "adapter_dimer_rate": {"min": 0.0, "max": .05, "warning_max": .10, "type": "percent"},
    "primary_mapping_rate": {"min": .50, "max": 1.0, "warning_min": .30, "type": "percent"},
    "usable_fragment_rate": {"min": .30, "max": 1.0, "warning_min": .15, "type": "percent"},
    "mitochondrial_fraction": {"min": 0.0, "max": .20, "warning_max": .40, "type": "percent"},
    "both_legs_mapping_rate": {"min": .50, "max": 1.0, "warning_min": .30, "type": "percent"},
    "cis_ratio": None,
}


def fmt(value):
    if value is None:
        return "NA"
    if isinstance(value, float):
        return f"{value:.4g}"
    return str(value)


def is_percent_metric(key, ranges):
    return key in PERCENT_METRICS or (ranges.get(key) or {}).get("type") == "percent"


def display_value(key, value, ranges):
    if key == "mycoplasma_status" and value == "NOT_DETECTED_AT_THIS_DEPTH":
        return "NOT_DETECTED"
    if value is None:
        return "NA"
    if is_percent_metric(key, ranges) and isinstance(value, (int, float)):
        if key in {"mycoplasma_qualifying_fraction", "spike_in_fraction"} and 0 < value < .001:
            return f"{value * 100:.3g}%"
        return f"{value:.1%}"
    if key.endswith("_bp") or key == "collision_bin_size":
        return f"{value:,.0f} bp" if isinstance(value, (int, float)) else str(value)
    if isinstance(value, int):
        return f"{value:,}"
    if isinstance(value, float):
        return f"{value:.4g}"
    return str(value)


def reference_interval_text(key, ranges, metric=None):
    spec = ranges.get(key) or {}
    minimum, maximum = spec.get("min"), spec.get("max")
    if minimum is not None and maximum is not None:
        return f"{display_value(key, minimum, ranges)}–{display_value(key, maximum, ranges)}"
    if minimum is not None:
        return f"≥{display_value(key, minimum, ranges)}"
    if maximum is not None:
        return f"≤{display_value(key, maximum, ranges)}"
    return reference_context(key, metric)[0]


def merged_reference_ranges(cfg):
    ranges = {key: dict(value) for key, value in DEFAULT_REFERENCE_RANGES.items() if value is not None}
    params = cfg.get("params", {})
    for key, value in (cfg.get("qc_reference_ranges") or {}).items():
        if value is None:
            ranges.pop(key, None)
        else:
            ranges[key] = {**ranges.get(key, {}), **value}
    ranges.pop("cis_ratio", None)
    return ranges


def range_status(key, value, ranges):
    spec = ranges.get(key)
    if not spec or not isinstance(value, (int, float)):
        return None
    if spec.get("min") is not None and value < spec["min"]:
        return "low"
    if spec.get("max") is not None and value > spec["max"]:
        return "high"
    return "normal"


def metric_light_class(key, value, ranges):
    if key in {"sequencing_yield_call", "adapter_dimer_call", "biological_qc_call"}:
        status = str(value or "").lower()
        if not status:
            return "gray"
        if "pass" in status or "high_quality" in status:
            return "green"
        if any(word in status for word in ("warning", "borderline", "inconclusive")):
            return "amber"
        return "red"
    if key == "spike_in_status":
        status = str(value or "").upper()
        if status == "NOT_CONFIGURED" or not status:
            return "gray"
        if status.startswith("INCONCLUSIVE"):
            return "amber"
        if status in {"CANDIDATE_SIGNAL", "TRACE"}:
            return "amber"
        if status in {"DETECTED", "NOT_DETECTED"}:
            return "gray"
        return "red"
    if key == "mycoplasma_status":
        status = str(value or "").upper()
        return "green" if status.startswith("NOT_DETECTED") else "red"
    if key in {"mycoplasma_qualifying_reads", "mycoplasma_qualifying_loci", "mycoplasma_qualifying_fraction"}:
        return "green" if isinstance(value, (int, float)) and value == 0 else "red"
    status = range_status(key, value, ranges)
    if status is None:
        return "gray"
    if status == "normal":
        return "green"
    spec = ranges[key]
    if status == "low" and spec.get("warning_min") is not None and value >= spec["warning_min"]:
        return "amber"
    if status == "high" and spec.get("warning_max") is not None and value <= spec["warning_max"]:
        return "amber"
    return "red"


def summary_overview_svg(metrics, ranges):
    """Technical quick overview: ONT sequencing and NGS library shape only."""
    columns = [
        ("Sample", 30, None),
        ("Assay", 130, None),
        ("Raw reads", 220, "raw_reads"),
        ("Mean Q", 365, "mean_read_q"),
        ("Reconstructed", 495, "complete_read_rate"),
        ("Raw N50", 650, "raw_n50_bp"),
        ("Insert median", 820, "insert_median_bp"),
        ("Mapping", 1000, "primary_mapping_rate"),
        ("Adapter dimer", 1165, "adapter_dimer_rate"),
    ]
    sample_extra = max(0, max((len(m["sample"]) for m in metrics), default=0) * 8 - 85)
    columns = [(label, x + (sample_extra if i else 0), key) for i, (label, x, key) in enumerate(columns)]
    width, row_h = 1370 + sample_extra, 52
    height = 135 + row_h * len(metrics)
    light_colors = {"green": "#16835D", "amber": "#E2A12A", "red": "#C23B3B", "gray": "#9BA8B2"}
    svg = [
        f"<svg xmlns='http://www.w3.org/2000/svg' viewBox='0 0 {width} {height}' role='img' aria-label='All samples technical QC summary' style='font-family:Arial,Helvetica,sans-serif'>",
        "<rect width='100%' height='100%' fill='#F5F8FA'/>",
        "<text x='30' y='40' font-size='28' font-weight='700' fill='#132A3E'>ONT sequencing and NGS library QC — all samples</text>",
        "<text x='30' y='65' font-size='13' fill='#607080'>Quick overview contains technical sequencing and library type/shape metrics only.</text>",
    ]
    for label, x, _ in columns:
        svg.append(f"<text x='{x}' y='102' font-size='12' font-weight='700' fill='#536475'>{html.escape(label)}</text>")
    for index, metric in enumerate(metrics):
        y = 132 + index * row_h
        if index % 2 == 0:
            svg.append(f"<rect x='18' y='{y-23}' width='{width-36}' height='{row_h}' rx='8' fill='white'/>")
        svg.extend([
            f"<text x='30' y='{y+5}' font-size='14' font-weight='700' fill='#18324A'>{html.escape(metric['sample'])}</text>",
            f"<text x='{130+sample_extra}' y='{y+5}' font-size='12' fill='#607080'>{html.escape(metric['assay'])}</text>",
        ])
        for _, x, key in columns[2:]:
            value = metric.get(key)
            color = light_colors["amber" if key == "adapter_dimer_rate" and str(metric.get("adapter_dimer_call", "")).startswith("inconclusive") else metric_light_class(key, value, ranges)]
            svg.append(f"<circle cx='{x+7}' cy='{y}' r='7' fill='{color}'/>")
            svg.append(f"<text x='{x+20}' y='{y+5}' font-size='11.5' fill='#34495E'>{html.escape(display_value(key, value, ranges))}</text>")
    svg.append(f"<text x='30' y='{height-9}' font-size='10' fill='#71808D'>Chromatin- and Micro-C-specific biological metrics are available in each sample report, not in this quick overview.</text></svg>")
    return "".join(svg)


REPORT_CSS = """
.extension-figure{margin:22px 0;max-width:1150px}.extension-figure img{width:100%;height:auto;display:block}.extension-figure figcaption{font-weight:600;margin-bottom:8px}
.table-scroll{overflow-x:auto;max-width:100%}.profile-scope{margin:24px 0 38px}
.plot-card{margin:0;min-width:0;padding:14px;background:#fff;border:1px solid #e3e9ee;border-radius:12px}
.analysis-plots .plot-card img{height:350px;object-fit:contain}
.chromosome-detail{border:1px solid #e3e9ee;border-radius:10px;padding:12px 16px;margin:0 0 24px;background:#fafcfd}
.chromosome-detail summary{cursor:pointer;color:#405e73;font-weight:600}.chromosome-detail p{font-size:13px;color:#667783;line-height:1.5}
.chromosome-detail img{width:100%;display:block;background:white;margin-top:12px}
.metric-status{border-bottom:0!important}.metric-status:focus-visible{outline:2px solid #527a91;outline-offset:4px;border-radius:3px}
@media(max-width:760px){.profile-scope .plots,.plots{grid-template-columns:minmax(0,1fr);gap:14px}.plot-card{padding:8px}.profile-scope .analysis-plots .plot-card img{height:auto}}
:root{font-family:Arial,Helvetica,sans-serif;color:#203040;background:#fff}
body{font-family:Arial,Helvetica,sans-serif;font-size:15px;max-width:1380px;margin:28px auto;padding:0 18px;color:#222}
h1,h2,h3,table,button,input{font-family:Arial,Helvetica,sans-serif}
a{color:#176B91}.plots{display:grid;grid-template-columns:repeat(2,minmax(0,1fr));gap:22px;margin-bottom:22px}.plots img,.summary-plot{width:100%;display:block}
table{border-collapse:collapse;width:100%;font-variant-numeric:tabular-nums}th,td{border-bottom:1px solid #dce2e7;padding:9px 10px;text-align:right;vertical-align:middle}th:first-child,td:first-child{text-align:left}th{background:#f2f4f6;color:#405464}tr:hover{background:#f8fafb}
.traffic-grid{display:grid;grid-template-columns:repeat(auto-fit,minmax(280px,1fr));gap:14px;margin:18px 0}.scope-card{border:1px solid #dce3e8;border-radius:12px;padding:10px 14px;background:#fafcfd;box-shadow:0 2px 10px #1730420d}.scope-card h3{margin:3px 0 9px}.traffic{display:grid;grid-template-columns:18px 1fr auto;gap:8px;align-items:center;padding:6px 0;border-top:1px solid #edf1f3}.traffic em{font-size:12px;font-style:normal;color:#536475}.light{display:inline-block;width:13px;height:13px;border-radius:50%;box-shadow:inset 0 0 0 1px #0002,0 0 0 3px #fff}.light.green{background:#16835D}.light.amber{background:#E2A12A}.light.red{background:#C23B3B}.light.gray{background:#9BA8B2}.note{padding:12px 14px;background:#f2f6f8;border-left:4px solid #527a91;line-height:1.55}.pass{color:#177245;font-weight:700}.fail{color:#b22222;font-weight:700}.warn{color:#b06b00;font-weight:700}code{background:#eee;padding:2px 4px}
.identity-panel{display:grid;grid-template-columns:minmax(260px,.8fr) minmax(420px,1.4fr);gap:14px;margin:16px 0 18px}.identity-card{border:1px solid #dce3e8;border-radius:12px;overflow:hidden;background:white;box-shadow:0 3px 14px #17304210}.identity-card h3{margin:0;padding:11px 14px;background:linear-gradient(100deg,#173f57,#2c7180);color:white;font-size:15px}.identity-grid{display:grid;grid-template-columns:1fr 1fr}.identity-item{padding:11px 14px;border-bottom:1px solid #edf1f3}.identity-item:nth-child(odd){border-right:1px solid #edf1f3}.identity-item small{display:block;color:#6a7a86;margin-bottom:4px}.identity-item strong{color:#173042}.mode-badge{display:inline-block;padding:3px 9px;border-radius:999px;background:#e3f1ec;color:#116145;font-size:12px}.pair-table th,.pair-table td{text-align:left;padding:8px 10px;font-size:13px}.pair-table th:last-child,.pair-table td:last-child{text-align:right}.cutoff-note{padding:8px 12px;color:#667783;font-size:11px;background:#f8fafb}
.metric-tip{cursor:help;border-bottom:1px dotted #607080}.metric-tip b{font-size:11px;color:#668094}
#metric-tooltip{position:fixed;z-index:10000;box-sizing:border-box;width:310px;max-width:calc(100vw - 16px);max-height:calc(100vh - 16px);overflow:auto;padding:10px 12px;border-radius:6px;background:#173042;color:white;font-size:12px;font-weight:400;line-height:1.45;box-shadow:0 4px 15px #0003;text-align:left;pointer-events:none}
.unavailable-reason{display:block;font-size:11px;font-weight:400;color:#667783;max-width:245px;line-height:1.4}
.reference-interval{white-space:nowrap;color:#586A77;font-size:13px}.metric-cell{min-width:155px}.metric-reading{display:flex;align-items:center;gap:8px;justify-content:flex-end;font-weight:650}
@media(max-width:850px){body{margin:12px auto}.identity-panel{grid-template-columns:1fr}.metric-cell{min-width:110px}th,td{padding:7px 5px;font-size:12px}}
"""

REPORT_CSS += Path(__file__).with_name("report_modern.css").read_text()


REPORT_TOOLTIP_SCRIPT = """<script>
(() => {
  const tip = document.createElement('div');
  tip.id = 'metric-tooltip'; tip.setAttribute('role', 'tooltip'); tip.hidden = true;
  document.body.appendChild(tip);
  let active = null;
  function hide() {
    if (active) active.removeAttribute('aria-describedby');
    active = null; tip.hidden = true;
  }
  function show(target) {
    hide(); active = target;
    tip.textContent = target.dataset.tip;
    tip.hidden = false; tip.style.left = '8px'; tip.style.top = '8px';
    target.setAttribute('aria-describedby', tip.id);
    const r = target.getBoundingClientRect(), t = tip.getBoundingClientRect();
    const w = document.documentElement.clientWidth, h = window.innerHeight;
    const left = Math.max(8, Math.min(r.left, w - t.width - 8));
    const below = r.bottom + 8;
    const top = below + t.height <= h - 8 ? below : Math.max(8, r.top - t.height - 8);
    tip.style.left = left + 'px'; tip.style.top = Math.min(top, Math.max(8, h - t.height - 8)) + 'px';
  }
  document.addEventListener('mouseover', e => { const t = e.target.closest('.metric-tip'); if (t && t !== active) show(t); });
  document.addEventListener('mouseout', e => { if (active && active.contains(e.target) && !active.contains(e.relatedTarget)) hide(); });
  document.addEventListener('focusin', e => { if (e.target.matches('.metric-tip')) show(e.target); });
  document.addEventListener('focusout', e => { if (e.target.matches('.metric-tip')) hide(); });
  document.addEventListener('click', e => { const t = e.target.closest('.metric-tip'); if (t) show(t); else hide(); });
  document.addEventListener('keydown', e => { if (e.key === 'Escape') hide(); });
  window.addEventListener('scroll', hide, true); window.addEventListener('resize', hide);
})();
</script>"""


def svg_region_pie(metric, title):
    values = [metric.get("promoter_proximal_rate"), metric.get("gene_body_rate"), metric.get("intergenic_rate")]
    labels = ["Promoter proximal", "Gene body", "Distal"]
    colors = ["#466E8A", "#71A8A0", "#D9BE8C"]
    width, height, cx, cy, radius = 640, 350, 180, 175, 105
    if any(v is None for v in values) or sum(values) <= 0:
        return f"<svg xmlns='http://www.w3.org/2000/svg' viewBox='0 0 {width} {height}' style='font-family:Arial,Helvetica,sans-serif'><rect width='100%' height='100%' fill='white'/><text x='{width/2}' y='28' text-anchor='middle' font-size='17'>{html.escape(title)}</text><text x='{width/2}' y='170' text-anchor='middle' fill='#687986'>No usable annotation data</text></svg>"
    total=sum(values); angle=-math.pi/2; paths=[]
    for label,value,color in zip(labels,values,colors):
        fraction=value/total; end=angle+2*math.pi*fraction
        x1,y1=cx+radius*math.cos(angle),cy+radius*math.sin(angle); x2,y2=cx+radius*math.cos(end),cy+radius*math.sin(end)
        large=1 if fraction>.5 else 0
        if fraction>=.999999:
            d=f"M {cx-radius} {cy} A {radius} {radius} 0 1 0 {cx+radius} {cy} A {radius} {radius} 0 1 0 {cx-radius} {cy}"
        else:
            d=f"M {cx} {cy} L {x1:.2f} {y1:.2f} A {radius} {radius} 0 {large} 1 {x2:.2f} {y2:.2f} Z"
        paths.append(f"<path d='{d}' fill='{color}' stroke='white' stroke-width='2'><title>{html.escape(label)}: {value:.1%}</title></path>")
        angle=end
    legend=''.join(f"<rect x='330' y='{112+i*48}' width='15' height='15' rx='2' fill='{c}'/><text x='354' y='{125+i*48}' font-size='13'>{html.escape(l)} {v:.1%}</text>" for i,(l,v,c) in enumerate(zip(labels,values,colors)))
    return f"<svg xmlns='http://www.w3.org/2000/svg' viewBox='0 0 {width} {height}' style='font-family:Arial,Helvetica,sans-serif'><rect width='100%' height='100%' fill='white'/><text x='{width/2}' y='28' text-anchor='middle' font-size='17'>{html.escape(title)}</text>{''.join(paths)}{legend}<text x='320' y='325' text-anchor='middle' font-size='11'>n={metric.get('annotation_mapped_reads',0):,} alignments; promoter priority</text></svg>"


def svg_contact_funnel(metric, title):
    stages = [
        ("Mappable units", metric.get("mappable_library_units")),
        ("Both legs mapped", metric.get("both_legs_mapped")),
        ("MAPQ pass", metric.get("mapq_pass_pairs")),
        ("Artifact filter pass", metric.get("same_fragment_filter_pass_pairs")),
        ("Unique endpoint contacts", metric.get("unique_valid_contacts")),
    ]
    width, height, left, top, bar_width, row_height = 900, 400, 235, 62, 600, 58
    values = [int(value or 0) for _, value in stages]
    peak = max(values[:1] + [1])
    elements = [
        f"<svg xmlns='http://www.w3.org/2000/svg' viewBox='0 0 {width} {height}' style='font-family:Arial,Helvetica,sans-serif'>",
        "<rect width='100%' height='100%' fill='white'/>",
        f"<text x='{width/2}' y='27' text-anchor='middle' font-size='17'>{html.escape(title)}</text>",
    ]
    previous = None
    colors = ["#356A91", "#3E7C8A", "#4A8D79", "#5A9C69", "#16835D"]
    for index, ((label, _), value, color) in enumerate(zip(stages, values, colors)):
        y = top + index * row_height
        length = bar_width * value / peak
        retention = None if previous in (None, 0) else value / previous
        suffix = "" if retention is None else f" · {retention:.1%} of previous"
        elements.extend([
            f"<text x='{left-14}' y='{y+18}' text-anchor='end' font-size='13' fill='#34495E'>{html.escape(label)}</text>",
            f"<rect x='{left}' y='{y}' width='{bar_width}' height='27' rx='6' fill='#E7ECEF'/>",
            f"<rect x='{left}' y='{y}' width='{length:.1f}' height='27' rx='6' fill='{color}'/>",
            f"<text x='{left+8}' y='{y+43}' font-size='12' fill='#34495E'>{value:,}{html.escape(suffix)}</text>",
        ])
        previous = value
    threshold = metric.get("long_range_cis_min_bp")
    long_fraction = metric.get("long_range_cis_fraction")
    note = "Long-range cis: NA"
    if threshold is not None and long_fraction is not None:
        note = f"Long-range cis (≥{threshold:,} bp / all unique valid contacts): {long_fraction:.1%}"
    elements.append(f"<text x='{left}' y='{height-18}' font-size='12' fill='#536475'>{html.escape(note)}</text></svg>")
    return "".join(elements)


def svg_contact_evidence(metric, title):
    width, height, left, span = 900, 340, 205, 630
    parts = [f"<svg xmlns='http://www.w3.org/2000/svg' viewBox='0 0 {width} {height}' style='font-family:Arial,Helvetica,sans-serif'>",
             "<rect width='100%' height='100%' fill='white'/>",
             f"<text x='450' y='26' text-anchor='middle' font-size='17'>{html.escape(title)}</text>",
             f"<text x='450' y='47' text-anchor='middle' font-size='11'>n={metric.get('unique_valid_contacts',0):,} unique endpoint contacts; Wilson 95% intervals</text>"]
    for i, (label, key, low, high) in enumerate([
        ("Cis / all contacts", "cis_ratio", "cis_ratio_ci95_low", "cis_ratio_ci95_high"),
        ("Long cis / all contacts", "long_range_cis_fraction", "long_range_cis_ci95_low", "long_range_cis_ci95_high")]):
        y = 85 + i*75
        parts.append(f"<text x='{left-12}' y='{y+4}' text-anchor='end' font-size='12'>{label}</text>")
        parts.append(f"<line x1='{left}' x2='{left+span}' y1='{y}' y2='{y}' stroke='#dce3e8' stroke-width='9'/>")
        value = metric.get(key)
        if value is not None and metric.get(low) is not None:
            x, xlo, xhi = [left+span*v for v in [value, metric[low], metric[high]]]
            parts.append(f"<line x1='{xlo:.2f}' x2='{xhi:.2f}' y1='{y}' y2='{y}' stroke='#245d7e' stroke-width='3'/><circle cx='{x:.2f}' cy='{y}' r='5' fill='#245d7e'/>")
            parts.append(f"<text x='{left}' y='{y+24}' font-size='12'>{value:.1%} [{metric[low]:.1%}, {metric[high]:.1%}]</text>")
        else:
            parts.append(f"<text x='{left}' y='{y+24}' font-size='12'>NA: no eligible contacts</text>")
    for i in range(6):
        parts.append(f"<text x='{left+i*span/5}' y='215' text-anchor='middle' font-size='11'>{i*20}%</text>")
    parts.append("<text x='450' y='289' text-anchor='middle' font-size='11'>Cis ratio alone does not validate ligation; short inward exclusion is heuristic.</text>")
    parts.append("<text x='450' y='311' text-anchor='middle' font-size='11'>Intervals assume independent contacts; ONT recovery biases are not represented.</text></svg>")
    return "".join(parts)


def call_class(call):
    call=call or ""
    if not call or "not_classified" in call or "NOT_CONFIGURED" in call: return "gray"
    if "pass" in call or "high_quality" in call: return "green"
    if any(x in call for x in ("warning","low_quality","inconclusive","insufficient","borderline")): return "amber"
    return "red"


def traffic_html(scope, ranges=None):
    ranges=ranges or DEFAULT_REFERENCE_RANGES
    label=scope.get("scope_name",scope.get("scope_id","Total"))
    seq_call=scope.get("sequencing_yield_call")
    shape_keys=["adapter_dimer_rate","primary_mapping_rate"]
    if scope.get("scope_id", "total") == "total": shape_keys.append("complete_read_rate")
    evaluated=[metric_light_class(k,scope.get(k),ranges) for k in shape_keys if scope.get(k) is not None and ranges.get(k)]
    if str(seq_call or "").startswith(("inconclusive", "insufficient")) or scope.get("raw_reads",0) < scope.get("metric_settings",metric_profile())["minimum_raw_reads"] or scope.get("library_units",0) < scope.get("metric_settings",metric_profile())["minimum_complete_units"]:
        shape_call="insufficient_evidence"
    elif "red" in evaluated:shape_call="fail"
    elif "amber" in evaluated:shape_call="warning"
    elif evaluated:shape_call="pass"
    else:shape_call=scope.get("adapter_dimer_call") or "insufficient_evidence"
    items=[("ONT Sequencing QC",seq_call,"Yield screen; inspect read quality and molecule recovery below."),("NGS Library Type / Shape QC",shape_call,"Combined adapter-dimer burden, primary mapping and total-library adapter recovery checks; not a biological enrichment or ligation success certificate.")]
    return f"<section class='scope-card'><h3>{html.escape(label)}</h3>"+''.join(f"<div class='traffic' title='{html.escape(tip,quote=True)}'><span class='light {call_class(v)}'></span><b>{html.escape(k)}</b><em>{html.escape(str(v or 'NA'))}</em></div>" for k,v,tip in items)+"</section>"


def metric_reading_html(rendered, light, state):
    return (f"<div class='metric-cell'><div class='metric-reading metric-tip metric-status' tabindex='0' data-tip='{html.escape(state, quote=True)}'>"
            f"<span class='light {light}' aria-hidden='true'></span><span>{html.escape(rendered)}</span></div></div>")


def scoped_value_html(key, value, ranges, descriptive=False):
    rendered = display_value(key, value, ranges)
    light = "gray" if descriptive else metric_light_class(key, value, ranges)
    if value == "NOT_CONFIGURED":
        state = "Not configured"
    elif descriptive:
        state = "Descriptive only"
    elif key == "mycoplasma_status":
        status = str(value or "").upper()
        if status.startswith("INCONCLUSIVE"):
            state = "Inconclusive: low screen depth"
        elif status.startswith("NOT_DETECTED"):
            state = "No supporting signal at this depth"
        else:
            state = "Mycoplasma-supporting signal detected" if light == "red" else "Not configured"
    elif key in {"mycoplasma_qualifying_reads", "mycoplasma_qualifying_loci", "mycoplasma_qualifying_fraction"}:
        state = "Mycoplasma-supporting signal detected" if light == "red" else "No supporting signal at this depth"
    elif key == "spike_in_status":
        state = "Candidate foreign-sequence evidence; intentional spike-in requires experimental metadata" if str(value) in {"CANDIDATE_SIGNAL", "TRACE"} else "No candidate signal at this depth" if value == "NOT_DETECTED" else "Reference screening result; experimental spike-in intent cannot be inferred"
    else:
        state = {"green":"Within range", "amber":"Borderline", "red":"Outside range", "gray":"Descriptive only"}[light]
    return metric_reading_html(rendered, light, state)



def metric_missing_reason(scope, key):
    reason = scope.get("metric_unavailable_reasons", {}).get(key)
    if reason:
        return reason
    if key.startswith("external_frip") and scope.get("external_frip_status") == "NOT_CONFIGURED":
        return "Not configured"
    if key.startswith("spike_in_") and scope.get("spike_in_status") == "NOT_CONFIGURED":
        return "Not configured"
    if key.startswith("mycoplasma_") and scope.get("mycoplasma_status") == "NOT_CONFIGURED":
        return "Not configured: Mycoplasma competitive screening is disabled."
    if key.startswith("mycoplasma_") and scope.get("mycoplasma_status") == "REQUIRES_SCREEN_RERUN":
        return "Rerun screen: per-insert supporting evidence is missing from the legacy result."
    return "Insufficient observations to estimate this metric."


def scope_metric_cell(scope, key, ranges):
    value = scope.get(key)
    if key in {"cis_ratio", "long_range_cis_fraction", "valid_contact_rate"} and value is not None:
        denominator=scope.get("unique_valid_contacts",0) if key!='valid_contact_rate' else scope.get('mappable_library_units',0)
        numerator=scope.get("cis_contacts",0) if key=='cis_ratio' else scope.get('long_range_cis_contacts',0) if key=='long_range_cis_fraction' else scope.get('unique_valid_contacts',0)
        text=f"{value:.1%}"
        if key!='valid_contact_rate':
            prefix='cis_ratio' if key=='cis_ratio' else 'long_range_cis'
            lo,hi=scope.get(prefix+'_ci95_low'),scope.get(prefix+'_ci95_high')
            if lo is not None and hi is not None:text+=f" [{lo:.1%}, {hi:.1%}]"
        text+=f" · {numerator:,}/{denominator:,}"
        evidence=scope.get('contact_evidence_status','')
        tip='Wilson 95% interval assumes independent contacts; count sufficiency does not establish biological success.'
        if evidence.startswith('insufficient'):tip='Insufficient evidence: too few valid contacts. '+tip
        return f"<td>{metric_reading_html(text,'amber' if evidence.startswith('insufficient') else 'gray',tip)}</td>"
    if key == "collision_oe" and value is not None:
        tip=f"{scope.get('collision_bin_size',10000)/1000:g} kb bins; {scope.get('collision_observed_pairs',0):,} observed pairs; {scope.get('collision_null_expected_pairs',0):.4g} expected under uniform-base sampling. Concentration is not target-specific enrichment."
        return f"<td>{metric_reading_html(display_value(key,value,ranges),'gray',tip)}</td>"
    if key.startswith("mycoplasma_") and value is not None and scope.get("mycoplasma_screened_reads") == 0:
        return f"<td>{metric_reading_html(display_value(key,value,ranges), 'red', 'No inserts screened; mycoplasma cannot be assessed.')}</td>"
    if key == "mycoplasma_status" and scope.get("mycoplasma_screened_reads", 0) and scope.get("mycoplasma_qualifying_fraction") == 0:
        value = "NOT_DETECTED_AT_THIS_DEPTH"
    if value is None:
        reason = metric_missing_reason(scope, key)
        label = "NOT_CONFIGURED" if reason.startswith("Not configured") else "Not attributable" if reason.startswith("Not attributable") else "Sample-level" if reason.startswith("Sample-level") else "Not estimable"
        missing_light='red' if key.startswith('mycoplasma_') else 'gray'
        return f"<td>{metric_reading_html(label, missing_light, reason)}</td>"
    descriptive = ((key == 'complete_read_rate' and scope.get('scope_id') != 'total')
                   or (key == 'cis_ratio' and str(scope.get('biological_qc_call', '')).startswith('inconclusive'))
                   or (key == 'adapter_dimer_rate' and str(scope.get('adapter_dimer_call', '')).startswith('inconclusive')))
    return f"<td>{scoped_value_html(key, value, ranges, descriptive)}</td>"

def scoped_table_html(scopes, keys, ranges):
    heads=''.join(f"<th>{html.escape(x.get('scope_name',x.get('scope_id','Total')))}</th>" for x in scopes)
    rows=[]
    for key in keys:
        label,desc=METRIC_META.get(key,(key,f"Pipeline QC metric: {key}."))
        interval=html.escape(reference_interval_text(key,ranges,scopes[0] if scopes else None))
        _,basis,detail=reference_context(key,scopes[0] if scopes else None)
        if ranges.get(key):
            basis="Local screening setting"
            detail="Numeric range from the effective pipeline configuration, not a universal assay standard."
            for boundary,boundary_label in [("warning_min","Lower warning boundary"),("warning_max","Upper warning boundary")]:
                if ranges[key].get(boundary) is not None:detail+=f" {boundary_label}: {display_value(key,ranges[key][boundary],ranges)}."
            if key=="complete_read_rate":detail+=" The threshold applies to total-library reconstruction. A sublibrary value is its contribution to all sample raw molecules, not a separate reconstruction success rate."
            if key in {"adapter_dimer_rate","primary_mapping_rate","usable_fragment_rate","both_legs_mapping_rate","mitochondrial_fraction","cis_ratio"}:
                detail+=" These ONT measurements and denominators require calibration against matched libraries and controls."
        interval=f"<span class='metric-tip reference-context' tabindex='0' data-tip='{html.escape(basis+': '+detail,quote=True)}'>{interval}</span>"
        cells=''.join(scope_metric_cell(x, key, ranges) for x in scopes)
        rows.append(f"<tr><td><span class='metric-tip' tabindex='0' data-tip='{html.escape(desc,quote=True)}'>{html.escape(label)} <b>ⓘ</b></span></td><td class='reference-interval'>{interval}</td>{cells}</tr>")
    return f"<div class='table-scroll'><table><tr><th>Metric</th><th>Reference interval</th>{heads}</tr>{''.join(rows)}</table></div>"


def pair_label(pair):
    return f"{pair.get('name',pair.get('scope_id','library'))}: i5 {pair.get('i5') or '∅'} / i7 {pair.get('i7') or '∅'} ({(pair.get('support_fraction') or 0):.1%})"


def identity_panel_html(metric, pairing, pairs):
    mode=str(pairing.get("mode","unresolved"))
    mode_label=mode.replace("_"," ").title()
    top=pairing.get("top_pair_fraction")
    top_text="NA" if top is None else f"{top:.1%}"
    pair_rows=[]
    for pair in pairs:
        pair_rows.append("<tr>" + "".join([
            f"<td><strong>{html.escape(str(pair.get('name',pair.get('scope_id','library'))))}</strong></td>",
            f"<td><code>{html.escape(str(pair.get('i5') or '∅'))}</code></td>",
            f"<td><code>{html.escape(str(pair.get('i7') or '∅'))}</code></td>",
            f"<td>{(pair.get('support_fraction') or 0):.1%}</td>"]) + "</tr>")
    if not pair_rows:
        pair_rows=["<tr><td colspan='4'>No confident i5–i7 pair resolved</td></tr>"]
    return ("<div class='identity-panel'><section class='identity-card'><h3>Library identity</h3><div class='identity-grid'>"
            f"<div class='identity-item'><small>Reference assembly</small><strong>{html.escape(str(metric.get('reference_assembly','unspecified')))}</strong></div>"
            f"<div class='identity-item'><small>Assay</small><strong>{html.escape(str(metric.get('assay','unknown')))}</strong></div>"
            f"<div class='identity-item'><small>Library mode</small><span class='mode-badge'>{html.escape(mode_label)}</span></div>"
            f"<div class='identity-item'><small>Top-pair support</small><strong>{html.escape(top_text)}</strong></div></div>"
            f"<div class='cutoff-note'>A dominant i5–i7 pair is single_library at ≥{pairing.get('single_library_fraction', .9):.0%} informative support, or when it is the only supported pair. Pair support is conditional on index-informative units; inspect overall assignment and ambiguous/unassigned fractions below.</div></section>"
            "<section class='identity-card'><h3>Resolved physical index pair(s)</h3><table class='pair-table'><tr><th>Library</th><th>i5 (P5 side)</th><th>i7 (P7 side)</th><th>Support</th></tr>"
            + "".join(pair_rows) + "</table></section></div>")



def profile_scope_html(label, raw_name, insert_name, region_name, karyotype_name, gene_name=None,
                       *, full_library_assets):
    def card(name, alt, extra_class=''):
        name = html.escape(name, quote=True)
        alt = html.escape(alt, quote=True)
        return f"<figure class='plot-card {extra_class}'><a href='assets/{name}' target='_blank' rel='noopener'><img src='assets/{name}' alt='{alt}'></a><figcaption class='figure-link'><a href='assets/{name}' target='_blank' rel='noopener'>Open figure ↗</a></figcaption></figure>"
    assets = full_library_assets
    info = assets['model_summary']
    offsets = info['adapter_offsets_bp']
    if len(offsets) == 1:
        definition = f"Full library = insert + {offsets[0]} bp adapter/index sequence."
    elif offsets:
        definition = "Full library = insert + adapter/index sequence; added lengths: " + ', '.join(map(str, offsets)) + " bp."
    else:
        definition = "Full-library length could not be modeled from the available adapter/index information."
    note = definition + f" {info['modeled_units']:,} modeled units."
    if info['excluded_units']:
        note += f" {info['excluded_units']:,} units lack the adapter/index information needed for this model."
    primary = card(assets['weighted'], 'ONT full-library profile (length-weighted)', 'length-weighted-profile')
    counts = card(assets['count'], 'ONT full-library profile (count-based)', 'full-library-count-card')
    counts += card(raw_name, 'Raw molecule size (count-based)') + card(insert_name, 'Insert size (count-based)')
    analysis = card(region_name, 'Genomic regions')
    if gene_name:
        analysis += card(gene_name, 'CPM-normalized TSS to TES profile')
    downloads = ' · '.join(f"<a href='assets/{html.escape(assets[key], quote=True)}'>{title}</a>" for key, title in
                           [('profiles', 'Download both profiles'), ('units', 'Modeled unit lengths'), ('model', 'Model details')])
    return (f"<section class='profile-scope'><h3>{html.escape(label)}</h3>"
            f"<p class='muted size-model-note'>{html.escape(note)} The main curve weights each unit by its full-library length.</p>"
            f"<div class='plots size-plots primary-size-plots'>{primary}</div>"
            f"<p class='muted size-profile-downloads'>{downloads}</p>"
            "<details class='count-based-details'><summary>Count-based profiles · full library, raw reads and inserts</summary>"
            "<p>Each recovered unit contributes one count to the full-library and insert profiles; each raw read contributes one count to the raw profile. Each curve has its own peak set to 100.</p>"
            f"<div class='plots size-plots count-size-plots'>{counts}</div></details>"
            f"<div class='plots analysis-plots'>{analysis}</div>"
            "<details class='chromosome-detail'><summary>Chromosome coverage · auxiliary diagnostic</summary>"
            "<p>Inspect broad coverage imbalance and isolated concentration. Gene density and shallow coverage alone cannot establish target-specific enrichment or copy-number changes.</p>"
            f"<img src='assets/{html.escape(karyotype_name, quote=True)}' alt='Mapped reads by chromosome karyotype'></details></section>")


def extension_report_html(out, sample, scopes, metric):
    sections=[]
    screen_path=out/"qc"/"spike_in"/f"{sample}.json"
    if screen_path.exists():
        screen=json.loads(screen_path.read_text())
        if screen.get("mode") == "auto":
            rows=[]
            for key,candidate in screen.get("candidates",{}).items():
                label=html.escape(candidate["label"]);status=candidate["status"]
                light="amber" if status in {"CANDIDATE_SIGNAL","TRACE"} else "gray"
                fraction=display_value("spike_in_fraction",candidate.get("fraction"),{})
                tooltip="Same species as the host; a spike-in cannot be distinguished" if status=="HOST_NOT_DISTINGUISHABLE" else "No spike-in was declared. This is sequence evidence only; a positive signal may have another source."
                rows.append(f"<tr><td>{label}</td><td>{html.escape(candidate['role'].replace('_',' '))}</td><td>{fmt(candidate.get('reads'))}</td><td>{fraction}</td><td>{fmt(candidate.get('loci'))}</td><td>{metric_reading_html(status,light,tooltip)}</td></tr>")
            declaration="Spike-in was declared in the sample configuration." if screen.get("declared_spike_in") else "No spike-in was added to these libraries."
            tip=declaration+" Candidate signals are exploratory foreign-sequence evidence; they do not establish intentional spike-in. Host, candidate species and microbial decoys compete in one reference. Same-species spike-in and unlisted synthetic standards are not identifiable by this panel."
            scoped_screen=scoped_table_html(scopes,["spike_in_status","spike_in_screened_reads","spike_in_reads","spike_in_fraction"],{}) if len(scopes)>1 else ""
            sections.append("<section class='report-section' id='spike-screen'><span class='eyebrow'>REFERENCE SCREENING</span><h2>Automatic spike-in candidate screen <span class='metric-tip info-icon' tabindex='0' data-tip='"+html.escape(tip,quote=True)+"'>?</span></h2><p class='muted'>"+fmt(screen.get("screened_reads"))+" screened inserts · "+("Spike-in declared" if screen.get("declared_spike_in") else "No spike-in declared")+"</p><div class='table-scroll'><table><tr><th>Reference candidate</th><th>Role</th><th>Inserts</th><th>Fraction</th><th>Supported loci</th><th>Result</th></tr>"+''.join(rows)+"</table></div>"+scoped_screen+f"<p><a href='../qc/spike_in/{sample}.assignments.tsv.gz'>Per-unit assignments</a> · <a href='../qc/spike_in/{sample}.json'>Screen details</a></p></section>")

    if metric["assay"]=="microc" and metric.get("ligation_screened_units") is not None:
        plots=[]
        for scope in scopes:
            scope_id=scope["scope_id"];source=out/"plots"/"ligation"/f"{sample}.{scope_id}.svg"
            if not source.exists():continue
            asset=out/"reports"/"assets"/f"{sample}.{scope_id}.ligation.svg";asset.write_bytes(source.read_bytes())
            plots.append(f"<figure class='extension-figure'><figcaption>{html.escape(scope.get('scope_name',scope_id))}</figcaption><img src='assets/{asset.name}' alt='Micro-C resolved genomic block counts and insert lengths'></figure>")
        source=out/"plots"/"ligation"/f"{sample}.examples.svg"
        if source.exists():
            asset=out/"reports"/"assets"/f"{sample}.ligation_examples.svg";asset.write_bytes(source.read_bytes());plots.append(f"<figure class='extension-figure'><img src='assets/{asset.name}' alt='Representative query block structures'></figure>")
        sections.append("<section class='report-section' id='ligation'><span class='eyebrow'>MOLECULAR STRUCTURE</span><h2>Micro-C ligation structure</h2><div class='note'>Standard MNase blunt-end ligation has no fixed junction sequence. Whole-insert split alignments estimate the minimum number of resolved genomic fragments within each adapter-bounded library unit. ONT adapter-linked units are counted separately. One block does not establish unligated DNA; adjacent, short or repetitive fragments can be missed. Unresolved units remain in the all-unit denominator. These are exploratory structure measurements, not a calibrated ligation-efficiency assay.</div>"+''.join(plots)+"<details class='metric-group'><summary>Ligation counts and denominators</summary>"+scoped_table_html(scopes,[key for key in METRIC_META if key.startswith("ligation_")],{})+"</details>"+f"<p><a href='../qc/ligation/{sample}.units.tsv.gz'>Every unit and inferred count</a> · <a href='../qc/ligation/{sample}.segments.tsv.gz'>Genomic blocks and query coordinates</a></p></section>")
    mnase_path=out/"qc/ligation"/f"{sample}.mnase.json"
    if metric["assay"]=="microc" and mnase_path.exists():
        mnase=json.loads(mnase_path.read_text());panels=[];rows=[]
        for scope in scopes:
            sid=scope["scope_id"];stats=mnase if sid=="total" else mnase["scopes"].get(sid,{})
            cells=[html.escape(scope.get("scope_name",sid)),fmt(stats.get("mnase_resolved_fragments")),fmt(stats.get("mnase_contributing_molecules")),fmt(stats.get("mnase_median_bp"))]
            for key in ["short_fraction","mono_range_fraction","intermediate_fraction","long_fraction"]:
                value=stats.get("mnase_"+key);cells.append(f"{value:.1%}" if value is not None else "NA")
            rows.append("<tr>"+''.join("<td>"+v+"</td>" for v in cells)+"</tr>")
            images=[]
            for suffix in ['mnase','single_block']:
                source=out/"plots/ligation"/f"{sample}.{sid}.{suffix}.svg"
                asset=out/"reports/assets"/source.name;asset.write_bytes(source.read_bytes())
                caption=("<strong>A · Constituent fragments</strong><p>Each observation is one resolved genomic block inside a unit with ≥2 blocks. For example, a 300-bp unit resolved as two 150-bp blocks contributes two observations near 150 bp.</p>" if suffix=='mnase' else "<strong>B · Whole inserts with one mapped block</strong><p>Each observation is an entire insert resolved as one continuous genomic block. It can be unligated DNA or a product whose junction is not distinguishable, such as adjacent genomic fragments.</p>")
                images.append(f"<figure class='plot-card'><figcaption class='fragment-caption'>{caption}</figcaption><a href='assets/{asset.name}' target='_blank' rel='noopener'><img src='assets/{asset.name}' alt='{suffix.replace('_',' ')} size trace'></a></figure>")
            panels.append(f"<h3>{html.escape(scope.get('scope_name',sid))}</h3><div class='plots size-plots'>"+''.join(images)+"</div>")
        sections.append("<section class='report-section' id='mnase'><span class='eyebrow'>FRAGMENTATION</span><h2>Inferred pre-ligation size &amp; MNase digestion</h2><p class='muted'>Reference-span lengths of resolved fragments in units with ≥2 genomic blocks. Panel A counts individual inferred fragments; panel B counts whole single-block inserts. They describe different observations.</p><div class='table-scroll'><table><tr><th>Library</th><th>Resolved fragments</th><th>Raw molecules</th><th>Median bp</th><th>&lt;100 bp</th><th>100–200 bp</th><th>200–300 bp</th><th>≥300 bp</th></tr>"+''.join(rows)+"</table></div>"+''.join(panels)+"<p class='muted trace-comparison-note'>Each curve has its own peak set to 100; its vertical grayscale lane shows the same count-weighted signal, with larger sizes at the top. Peak height and lane darkness cannot compare abundance between panels. Use the observation counts and ligation table for quantities. These are ONT-derived traces, not measured TapeStation fluorescence.</p><details><summary>How to interpret digestion evidence</summary><p>"+html.escape(mnase['limitations'])+" Size ranges are half-open: 100–200 means 100 ≤ length &lt; 200 bp. The 100–200 bp fraction is mononucleosome-size evidence, not proof of nucleosome identity or a calibrated digestion pass threshold. Blocks below the alignment minimum (40 bp) are not observable. A direct pre-ligation MNase titration trace is still needed to assess the original digestion mixture.</p></details>"+f"<p><a href='../qc/ligation/{sample}.mnase_fragments.tsv.gz'>Fragment lengths and parent molecules</a> · <a href='../qc/ligation/{sample}.mnase.json'>MNase size metrics</a></p></section>")
    return ''.join(sections)


def cmd_report(args):
    import csv
    import yaml
    with open(args.config) as handle: cfg=yaml.safe_load(handle)
    samples=list(cfg["samples"]); ranges=merged_reference_ranges(cfg); out=Path(args.out)
    (out/"reports"/"assets").mkdir(parents=True,exist_ok=True); (out/"plots"/"summary").mkdir(parents=True,exist_ok=True); (out/"summary").mkdir(parents=True,exist_ok=True)
    metrics=[]
    report_tools=cfg.get("tools", {})
    report_params=configured_params(cfg)
    report_samtools=report_tools.get("samtools", "samtools")
    report_rscript=report_tools.get("rscript", "Rscript")
    karyotype_script=Path(__file__).with_name("karyotype_plot.R")
    promoter_window=int(report_params.get("annotation_promoter_window_bp", 2000))
    report_mapq=int(report_params.get("mapq", 20))
    karyotype_bins=int(report_params.get("karyotype_bins_per_longest_chromosome", 250))
    annotation_cache={}
    for sample in samples:
        with open(out/"qc"/"metrics"/f"{sample}.json") as handle: metric=json.load(handle)
        if metric.get("assay") == "chromatin":
            metric["biological_qc_call"] = "not_classified"
            metric["classification_note"] = "Chromatin enrichment is not classified. Collision O/E is depth normalized (random baseline 1). Sampling uncertainty remains at low depth. TSS–TES profiles are CPM normalized and visual only."
        metrics.append(metric)
        raw=read_numbers(out/"qc"/"lengths"/f"{sample}.raw.tsv"); inserts=read_numbers(out/"qc"/"lengths"/f"{sample}.insert.tsv")
        scope_profiles=read_scope_profiles(out/"qc"/"adapters"/f"{sample}.reads.tsv.gz", include_total=True)
        with open(out/"qc"/"adapters"/f"{sample}.summary.json") as handle:
            adapter_summary=json.load(handle)
        metric["scope_id"]="total"; metric["scope_name"]="Total library"
        metric["scope_unit_fraction"] = 1.0 if metric["library_units"] else None
        scopes=[metric]+metric.get("sub_libraries",[])
        if metric.get("assay") == "chromatin":
            for scope in scopes:
                scope["biological_qc_call"] = "not_classified"
        reference_key=cfg["samples"][sample].get("reference", cfg.get("default_reference", "dm6"))
        report_reference=cfg["references"][reference_key]
        report_gene_bed=out/"reference"/f"{reference_key}.genes.bed"
        report_chrom_sizes=read_chrom_sizes(report_reference["chrom_sizes"])
        primary_arg=report_reference.get("primary_chromosomes", "auto")
        primary_arg=primary_arg if isinstance(primary_arg, str) else ",".join(primary_arg)
        mito_arg=report_reference.get("mitochondrial_chromosomes", "auto")
        mito_arg=mito_arg if isinstance(mito_arg, str) else ",".join(mito_arg)
        report_primary, _=resolve_chromosome_sets(report_chrom_sizes, primary_arg, mito_arg)
        report_annotation=report_reference.get("annotation", "")
        annotation_key=(report_annotation, tuple(report_chrom_sizes.items()), tuple(report_primary), promoter_window)
        if annotation_key not in annotation_cache:
            annotation_cache[annotation_key]=load_gene_annotation(report_annotation, report_chrom_sizes, promoter_window)
        report_promoters, report_genes=annotation_cache[annotation_key]
        report_assignments={}
        scope_allowed={}
        if metric.get("sub_libraries"):
            report_assignments=read_index_assignments(out/"qc"/"adapters"/f"{sample}.index_assignments.tsv.gz")
            for sub_scope in metric["sub_libraries"]:
                scope_allowed[sub_scope["scope_id"]]={unit_id for unit_id, row in report_assignments.items() if row.get("scope_id") == sub_scope["scope_id"] and row.get("assignment_status") == "assigned" and row.get("is_adapter_dimer") != "1"}
        report_records=sam_primary_records(str(out/"bam"/f"{sample}.bam"), report_samtools)
        scope_records={"total": report_records}
        for sub_scope in metric.get("sub_libraries", []):
            allowed=scope_allowed.get(sub_scope["scope_id"], set())
            scope_records[sub_scope["scope_id"]]={qname: rec for qname, rec in report_records.items() if base_unit_id(qname) in allowed}
        region_metrics=annotation_region_counts(report_records, report_promoters, report_genes, report_primary, report_mapq, promoter_window)
        sub_region_metrics={}
        for sub_scope in metric.get("sub_libraries", []):
            scoped_records=scope_records[sub_scope["scope_id"]]
            sub_region_metrics[sub_scope["scope_id"]]=annotation_region_counts(scoped_records, report_promoters, report_genes, report_primary, report_mapq, promoter_window)
        gene_profile_assets={}
        if metric["assay"] == "chromatin":
            gene_profile_source=out/"plots"/"gene_profile"/f"{sample}.TSS_TES.svg"
            gene_profile_asset=out/"reports"/"assets"/f"{sample}.TSS_TES.svg"
            gene_profile_asset.write_bytes(gene_profile_source.read_bytes())
            gene_profile_assets["total"]=gene_profile_asset
            if metric.get("sub_libraries"):
                genes_path=report_gene_bed
                profile_threads=int(report_params.get("gene_profile_sublibrary_threads", min(16, os.cpu_count() or 2)))
                profile_upstream=int(report_params.get("gene_profile_upstream_bp", 5000))
                profile_downstream=int(report_params.get("gene_profile_downstream_bp", 5000))
                profile_body=int(report_params.get("gene_profile_body_bp", 5000))
                profile_bin=int(report_params.get("gene_profile_bin_size", 100))
                for sub_scope in metric["sub_libraries"]:
                    scope_id=sub_scope["scope_id"]
                    cached_asset=out/"reports"/"assets"/f"{sample}.{scope_id}.TSS_TES.svg"
                    if getattr(args, "reuse_plot_assets", False) and cached_asset.exists():
                        gene_profile_assets[scope_id]=cached_asset
                        continue
                    gene_profile_assets[scope_id]=build_sublibrary_gene_profile(
                        out, sample, scope_id, sub_scope.get("scope_name", scope_id),
                        out/"bam"/f"{sample}.bam", genes_path, scope_allowed.get(scope_id, set()),
                        report_samtools, report_tools.get("bedGraphToBigWig", "bedGraphToBigWig"),
                        report_tools.get("computeMatrix", "computeMatrix"), report_tools.get("plotProfile", "plotProfile"),
                        report_mapq, report_reference["chrom_sizes"], report_chrom_sizes, report_primary,
                        profile_upstream, profile_downstream, profile_body, profile_bin, profile_threads)
        profile_sections=[]
        for scope in scopes:
            scope_id=scope["scope_id"]; label=scope.get("scope_name",scope_id)
            profile={"raw_lengths":raw,"insert_lengths":inserts} if scope_id=="total" else scope_profiles.get(scope_id, {"raw_lengths":[],"insert_lengths":[]})
            raw_asset=out/"reports"/"assets"/f"{sample}.{scope_id}.raw_lengths.svg"
            insert_asset=out/"reports"/"assets"/f"{sample}.{scope_id}.insert_lengths.svg"
            region_asset=out/"reports"/"assets"/f"{sample}.{scope_id}.regions.svg"
            karyotype_asset=out/"reports"/"assets"/f"{sample}.{scope_id}.karyotype.svg"
            trace_max=float(report_params.get("size_trace_max_bp", 2000))
            trace_bandwidth=float(report_params.get("size_trace_bandwidth", .045))
            raw_asset.write_text(svg_size_trace(profile["raw_lengths"],f"{label}: raw ONT molecule size",trace_max,trace_bandwidth))
            insert_asset.write_text(svg_size_trace(profile["insert_lengths"],f"{label}: adapter-bounded insert size",trace_max,trace_bandwidth))
            full_library_assets=write_library_size_assets(out/"reports"/"assets", f"{sample}.{scope_id}", label,
                scope_profiles.get(scope_id, {}).get("unit_rows", []), adapter_summary, trace_max, trace_bandwidth)
            region_metric=region_metrics if scope_id == "total" else sub_region_metrics.get(scope_id, annotation_region_counts({}, report_promoters, report_genes, report_primary, report_mapq, promoter_window))
            region_asset.write_text(svg_region_pie(region_metric,f"{label}: genomic regions"))
            if not getattr(args, "reuse_plot_assets", False) or not karyotype_asset.exists():
                write_karyoploter(
                    scope_records.get(scope_id, {}), report_chrom_sizes, report_primary, report_mapq,
                    karyotype_bins, report_rscript, karyotype_script, report_gene_bed,
                    karyotype_asset, f"{label}: gene and mapped-read density",
                    metric.get("reference_assembly", "custom genome"),
                    gene_window_bp=int(report_params.get("karyotype_gene_window_bp", 250000)))
            profile_sections.append(profile_scope_html(label, raw_asset.name, insert_asset.name,
                                                       region_asset.name, karyotype_asset.name,
                                                       gene_profile_assets[scope_id].name if scope_id in gene_profile_assets else None,
                                                       full_library_assets=full_library_assets))
        profiles_html=''.join(profile_sections)
        funnel_assets=[]
        if metric["assay"]=="microc":
            for scope in scopes:
                asset=out/"reports"/"assets"/f"{sample}.{scope['scope_id']}.contact_funnel.svg"
                asset.write_text(svg_contact_funnel(scope,f"{scope.get('scope_name',scope['scope_id'])}: valid-contact funnel"))
                funnel_assets.append(asset.name)
                evidence=out/"reports"/"assets"/f"{sample}.{scope['scope_id']}.contact_evidence.svg"
                evidence.write_text(svg_contact_evidence(scope, f"{scope.get('scope_name',scope['scope_id'])}: contact evidence"))
                
                distance_asset=out/"reports"/"assets"/f"{sample}.{scope['scope_id']}.contact_distances.svg"
                distance_asset.write_text(svg_bars(["<1 kb", "1-10 kb", "10-100 kb", ">=100 kb", "trans"],
                    [scope.get(k,0) for k in ["cis_lt1kb","cis_1_10kb","cis_10_100kb","cis_ge100kb","trans_contacts"]],
                    f"{scope.get('scope_name',scope['scope_id'])}: contact counts (not P(s))"))
                funnel_assets.append(distance_asset.name)
        pairing=metric.get("index_pairing",{}); pairs=pairing.get("pairs",[])
        if pairing.get("mode")=="single_library" and pairs: pairs=[max(pairs,key=lambda x:x.get("informative_support",0))]
        top=''.join(traffic_html(x, ranges) for x in scopes)
        ont_keys=["raw_reads","sequencing_yield_call","mean_read_q","raw_median_bp","raw_n50_bp","complete_read_rate"]
        library_keys=(["scope_unit_fraction"] if metric.get("sub_libraries") else [])+["outer_complete_units","outer_complete_unit_rate","index_assignment_rate","index_ambiguous_fraction","index_unassigned_fraction","unreconstructed_probable_dimer_reads","adapter_dimer_unit_rate","concatemer_rate","insert_median_bp","library_units","units_per_complete_read","mappable_library_units","canonical_adapter_dimer_units","probable_adapter_dimer_units","adapter_dimer_rate","raw_structure_dimer_rate","adapter_dimer_call","primary_mapping_rate"]
        screen_keys=["mycoplasma_status","mycoplasma_screened_reads","mycoplasma_qualifying_reads","mycoplasma_qualifying_loci","mycoplasma_qualifying_fraction"]
        if metric["assay"]=="chromatin":
            type_specific_keys=["collision_oe","collision_contributing_molecules","usable_fragments","usable_fragment_rate","mitochondrial_fraction"]+screen_keys
            type_specific_title="Chromatin-specific QC"
        else:
            type_specific_keys=["unique_valid_contacts","cis_ratio","long_range_cis_fraction","valid_contact_rate"]+screen_keys
            type_specific_title="Micro-C-specific QC"
        ont_table=scoped_table_html(scopes,ont_keys,ranges)
        library_key_metrics=(["scope_unit_fraction"] if metric.get("sub_libraries") else [])+["adapter_dimer_call","adapter_dimer_rate","insert_median_bp","outer_complete_unit_rate","index_assignment_rate","primary_mapping_rate","mappable_library_units"]
        library_primary_table=scoped_table_html(scopes,library_key_metrics,ranges)
        library_table=scoped_table_html(scopes,[key for key in library_keys if key not in library_key_metrics],ranges)
        type_specific_table=scoped_table_html(scopes,type_specific_keys,ranges)
        funnel_html=''.join(f"<img src='assets/{x}' alt='Micro-C contact QC plot'>" for x in funnel_assets)
        if metric['assay']=='microc':
            detail_keys=["short_unmapped_units","both_legs_mapped","both_legs_mapping_rate","mapq_pass_pairs","mapq_pass_pair_rate","same_fragment_pairs","same_fragment_filter_pass_pairs","same_fragment_filter_pass_rate","contact_deduplication_retention","cis_contacts","trans_contacts","long_range_cis_contacts","long_range_within_cis_fraction","cis_lt1kb","cis_1_10kb","cis_10_100kb","cis_ge100kb"]
            funnel_html += "<details class='metric-group'><summary>Contact counts and filtering details</summary>"+scoped_table_html(scopes,detail_keys,{})+"</details>"
        settings=metric.get('metric_settings',metric_profile())
        definition_label=f"{settings['metric_profile']} · MAPQ ≥{settings['mapq']}"
        if metric['assay']=='chromatin':definition_label+=f" · {settings['collision_qc_bin_size']/1000:g} kb bins"
        else:definition_label+=f" · {settings['dedup_tolerance_bp']} bp deduplication"
        from custom_regions import report_html as custom_region_report
        custom_html=custom_region_report(out,sample,reference_key,report_chrom_sizes,report_primary,scope_records,scopes,cfg['samples'][sample],Path(__file__).resolve().parents[2],report_mapq) if metric['assay']=='chromatin' else ""
        report=f"""<!doctype html><html lang='en'><head><meta charset='utf-8'><meta name='viewport' content='width=device-width,initial-scale=1'><title>{sample} · ONT library QC</title><style>{REPORT_CSS}</style></head><body>
<header class='report-hero'><span class='eyebrow'>ONT LIBRARY QUALITY CONTROL</span><h1>{sample}</h1><p>Shallow sequencing · whole-molecule library assessment</p><div class='hero-meta'><span>{html.escape(metric['assay'])}</span><span>{html.escape(reference_key)}</span><span>{metric['raw_reads']:,} raw molecules</span><span>{html.escape(definition_label)}</span></div></header>
<nav class='report-nav'><a href='#sequencing'>Sequencing</a><a href='#library-shape'>NGS library QC</a><a href='#assay-qc'>Assay metrics</a>{"<a href='#custom-regions'>Custom regions</a>" if custom_html else "<a href='#ligation'>Ligation</a><a href='#mnase'>MNase</a>"}{"<a href='#spike-screen'>Candidate screen</a>" if metric.get('spike_in_mode')=='auto' or json.loads((out/'qc/spike_in'/f'{sample}.json').read_text()).get('mode')=='auto' else ''}<a href='../summary/index.html'>All samples ↗</a></nav>
<div class='traffic-grid'>{top}</div>
<section class='report-section' id='sequencing'><span class='eyebrow'>01 · SEQUENCING</span><h2>ONT Sequencing QC</h2>{ont_table}</section>
<section class='report-section' id='library-shape'><span class='eyebrow'>02 · NGS LIBRARY QC</span><h2>NGS Library Type / Shape QC</h2>{identity_panel_html(metric,pairing,pairs)}
<div class='adapter-context'><strong>P5/P7-compatible adapter family</strong> · <code>{html.escape(metric['adapter_family'])}</code><p class='muted'>Physical i5 on the P5 side and i7 on the P7 side; sequences are normalized to P5→P7 orientation.</p></div>
<div id='library-key-metrics'>{library_primary_table}</div>
<details class='metric-group' id='library-extra-metrics'><summary>Additional structure, index and reconstruction metrics</summary>{library_table}</details>
{reference_notes_html()}
<details><summary>Size and adapter interpretation</summary><p>The 120–170 bp raw-size window is supportive only. Dimer-like units have a ≤{fmt(metric.get('adapter_dimer_max_bp'))}-bp central gap. Unreconstructed template matches remain candidates; resolved 11–50 bp inserts are retained. Full-library sizes include expected outer adapters and inferred physical indexes, including units with undetected outer handles. The main profile weights each recovered unit by its modeled full-library length; count-based profiles are expandable. Curves are conditional on sequencing and adapter recovery and do not measure TapeStation fluorescence or NGS molarity.</p></details>
<div class='library-profiles' id='profiles'><h3>Library profiles · size &amp; genomic distribution</h3>{profiles_html}</div></section>
<section class='report-section' id='assay-qc'><span class='eyebrow'>03 · ASSAY EVIDENCE</span><h2>{type_specific_title}</h2>{type_specific_table}<details><summary>Assay interpretation</summary><p>{html.escape(metric.get('classification_note',''))}</p></details>{"<h3>Micro-C valid-contact funnel</h3><div class='plots'>"+funnel_html+"</div>" if funnel_html else ""}</section>
{custom_html}{extension_report_html(out, sample, scopes, metric)}
<footer class='report-footer'><details><summary>Methods &amp; interpretation</summary><p>Reference values distinguish local screening settings, biological context and statistical expectations. Gene profiles show mean coverage and SEM across annotated genes; SEM is not biological replicate uncertainty. CPM uses MAPQ-pass configured primary nuclear alignments in each scope. Region fractions and chromosome density use mapped alignments (Micro-C: endpoint legs). PCR-duplicate rate is omitted at this shallow depth; Micro-C contact deduplication remains an internal validity filter.</p></details><p><a href='../summary/index.html'>All samples</a> · <a href='../qc/adapters/{sample}.reads.tsv.gz'>Adapter table</a> · <a href='../qc/adapters/{sample}.index_assignments.tsv.gz'>Index assignments</a></p></footer>{REPORT_TOOLTIP_SCRIPT}</body></html>"""
        (out/"reports"/f"{sample}.html").write_text(report)
    columns=["sample","assay","reference_assembly","mean_read_q","outer_complete_unit_rate","index_assignment_rate","index_ambiguous_fraction","index_unassigned_fraction","cis_ratio_ci95_low","cis_ratio_ci95_high","raw_reads","sequencing_yield_call","complete_read_rate","adapter_dimer_rate","adapter_dimer_call","index_architecture","primary_mapping_rate","spike_in_reads","spike_in_fraction","collision_bin_size","collision_oe","collision_fragments","collision_observed_pairs","collision_null_expected_pairs","usable_fragment_rate","both_legs_mapping_rate","mapq_pass_pairs","same_fragment_filter_pass_pairs","unique_valid_contacts","valid_contact_rate","cis_ratio","long_range_cis_min_bp","long_range_cis_contacts","long_range_cis_fraction","contact_evidence_status","metric_profile","metric_definition_version"]
    for m in metrics:
        m.update({k:m.get("metric_settings",{}).get(k) for k in ("metric_profile","metric_definition_version")})
    with open(out/"summary"/"qc_metrics.tsv","w") as handle:
        handle.write('\t'.join(columns)+'\n')
        for m in metrics: handle.write('\t'.join(fmt(m.get(k)) for k in columns)+'\n')
    # Include structure fields in a dedicated numeric table without rounding source counts.
    ligation_columns=["sample","scope_id","scope_name"]+[key for key in METRIC_META if key.startswith("ligation_")]
    with open(out/"summary"/"ligation_qc_metrics.tsv", "w") as handle:
        handle.write("\t".join(ligation_columns)+"\n")
        for m in metrics:
            if m["assay"] != "microc": continue
            for scope in [dict(m,scope_id="total",scope_name="Total library")]+m.get("sub_libraries",[]):
                row={"sample":m["sample"],**scope};handle.write("\t".join(fmt(row.get(k)) for k in ligation_columns)+"\n")
    mnase_rows=[]
    for path in sorted((out/"qc/ligation").glob("*.mnase.json")):
        value=json.loads(path.read_text())
        for sid,stats in [("total",value)]+list(value["scopes"].items()):
            mnase_rows.append(dict(sample=value["sample"],scope_id=sid,**{k:v for k,v in stats.items() if k.startswith("mnase_")}))
    with open(out/"summary/mnase_qc_metrics.tsv","w") as handle:
        writer=csv.DictWriter(handle,list(mnase_rows[0]) if mnase_rows else ["sample","scope_id"],delimiter="\t");writer.writeheader();writer.writerows(mnase_rows)
    scoped_columns = ["sample", "scope_id", "scope_name", "assay", "reference_assembly"] + list(METRIC_META)
    with open(out/"summary"/"sublibrary_qc_metrics.tsv", "w") as handle:
        handle.write("\t".join(scoped_columns) + "\n")
        for m in metrics:
            for sub in m.get("sub_libraries", []):
                row = {"sample": m["sample"], "assay": m["assay"], "reference_assembly": m["reference_assembly"], **sub}
                handle.write("\t".join(fmt(row.get(k)) for k in scoped_columns) + "\n")
    summary_plot=out/"plots"/"summary"/"all_samples_qc.svg"; summary_plot.write_text(summary_overview_svg(metrics,ranges))
    rows=''.join(f"<tr><td><a href='../reports/{m['sample']}.html'>{html.escape(m['sample'])}</a></td><td>{html.escape(m['assay'])}</td><td>{html.escape(str(m.get('reference_assembly','')))}</td><td>{m['raw_reads']:,}</td><td>{display_value('mean_read_q',m.get('mean_read_q'),ranges)}</td><td>{display_value('complete_read_rate',m.get('complete_read_rate'),ranges)}</td><td>{display_value('raw_n50_bp',m.get('raw_n50_bp'),ranges)}</td><td>{display_value('insert_median_bp',m.get('insert_median_bp'),ranges)}</td><td>{display_value('primary_mapping_rate',m.get('primary_mapping_rate'),ranges)}</td><td>{display_value('adapter_dimer_rate',m.get('adapter_dimer_rate'),ranges)}</td></tr>" for m in metrics)
    index=f"""<!doctype html><html><head><meta charset='utf-8'><title>ONT library QC</title><style>{REPORT_CSS}</style></head><body><h1>ONT full-molecule QC summary</h1><div class='note'>This quick overview is intentionally limited to ONT sequencing QC and NGS library type/shape QC. Chromatin- and Micro-C-specific metrics are organized in the third QC section of each sample report.</div><img class='summary-plot' src='../plots/summary/all_samples_qc.svg' alt='All samples technical QC summary'><table><tr><th>Sample</th><th>Assay</th><th>Assembly</th><th>Raw reads</th><th>Mean Q</th><th>Complete</th><th>Raw N50</th><th>Insert median</th><th>Mapping</th><th>Adapter dimer</th></tr>{rows}</table><p><a href='qc_metrics.tsv'>Download total-library TSV</a> · <a href='sublibrary_qc_metrics.tsv'>Download sublibrary TSV</a></p>{REPORT_TOOLTIP_SCRIPT}</body></html>"""
    (out/"summary"/"index.html").write_text(index)
    json_dump({"created_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
               "config": cfg, "python_version": sys.version,
               "metric_schema_version": 3,
               "metric_definitions": {"native": metric_profile(), "matched_comparison": metric_profile("matched_comparison")},
               "library_size_profiles": {"main_weighting": "length", "count_based_default_collapsed": True,
                   "length_definition": "insert + sequence-derived adapters and physical index lengths",
                   "maximum_bp": float(report_params.get("size_trace_max_bp", 2000)),
                   "bandwidth": float(report_params.get("size_trace_bandwidth", .045)),
                   "coordinate": "log1p(bp/50)", "measured_fluorescence": False}}, out/"summary"/"provenance.json")

def cmd_spikein(args):
    import re
    result = {
        "label": args.label,
        "status": "NOT_CONFIGURED",
        "interpretation": "Optional host+spike-in competitive-alignment screen.",
        "qualifying_reads": None,
        "qualifying_fraction": None,
        "qualifying_read_ids": [],
        "minimum_screened_reads": args.minimum_screened_reads,
    }
    if not args.enabled:
        json_dump(result, args.output)
        print(json.dumps(result))
        return
    if not args.index or not args.target_regex:
        raise ValueError("An enabled spike-in screen requires a combined host+spike-in Bowtie2 index and target contig regex")
    command = [args.bowtie2, "--very-sensitive-local", "-x", args.index, "-U", args.fastq]
    process = subprocess.Popen(command, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    target = re.compile(args.target_regex, re.I)
    hits = set()
    total = 0
    for line in process.stdout:
        if line.startswith("@"):
            continue
        fields = line.rstrip().split("	")
        total += 1
        flag = int(fields[1])
        mapq = int(fields[4])
        cigar = fields[5]
        sequence = fields[9]
        if flag & (4 | 256 | 2048) or not target.search(fields[2]) or mapq < args.mapq:
            continue
        aligned = sum(int(n) for n, op in re.findall(r"(\d+)([M=X])", cigar))
        edit_distance = next((int(field[5:]) for field in fields[11:] if field.startswith("NM:i:")), 999999)
        coverage = aligned / max(1, len(sequence))
        identity = (aligned - edit_distance) / max(1, aligned)
        if aligned >= args.minimum_aligned_bp and coverage >= args.minimum_query_coverage and identity >= args.minimum_identity:
            hits.add(fields[0])
    stderr = process.stderr.read()
    if process.wait() != 0:
        raise RuntimeError(f"Bowtie2 spike-in screen failed: {stderr[-1000:]}")
    count = len(hits)
    fraction = count / max(1, total)
    if count:
        status = "DETECTED"
    elif total < args.minimum_screened_reads:
        status = "INCONCLUSIVE_LOW_SCREEN_DEPTH"
    else:
        status = "NOT_DETECTED"
    result.update(
        status=status,
        screened_reads=total,
        qualifying_reads=count,
        qualifying_fraction=fraction,
        qualifying_read_ids=sorted(hits),
        interpretation="Spike-in fraction is target-assigned inserts divided by all inserts in a host+spike-in competitive alignment.",
    )
    json_dump(result, args.output)
    print(json.dumps(result, indent=2))


def cmd_mycoplasma(args):
    import re
    result = {"status": "NOT_CONFIGURED", "interpretation": "Opportunistic screen only; zero hits never excludes contamination.", "qualifying_reads": None, "qualifying_loci": None, "qualifying_fraction": None,
              "minimum_screened_reads": args.minimum_screened_reads, "qualifying_read_loci": {}}
    if not args.enabled:
        json_dump(result, args.output); print(json.dumps(result)); return
    if not args.index:
        raise ValueError("An enabled Mycoplasma screen requires a manually built combined host+Mollicutes+decoy Bowtie2 index")
    command=[args.bowtie2,"--very-sensitive-local","-x",args.index,"-U",args.fastq]
    process=subprocess.Popen(command,text=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    target=re.compile(args.target_regex,re.I); hits=set(); loci=set(); total=0
    read_loci = collections.defaultdict(set)
    for line in process.stdout:
        if line.startswith("@"): continue
        f=line.rstrip().split("\t"); total += 1
        flag=int(f[1]); mapq=int(f[4]); cigar=f[5]; seq=f[9]
        if flag & (4|256|2048) or not target.search(f[2]) or mapq < args.mapq: continue
        aligned=sum(int(n) for n,op in re.findall(r"(\d+)([M=X])",cigar)); nm=next((int(x[5:]) for x in f[11:] if x.startswith("NM:i:")),999999)
        coverage=aligned/max(1,len(seq)); identity=(aligned-nm)/max(1,aligned)
        if aligned >= args.minimum_aligned_bp and coverage >= args.minimum_query_coverage and identity >= args.minimum_identity:
            hits.add(f[0]); loci.add((f[2],int(f[3])//100))
            read_loci[f[0]].add((f[2],int(f[3])//100))
    stderr=process.stderr.read()
    if process.wait()!=0: raise RuntimeError(f"Bowtie2 Mycoplasma screen failed: {stderr[-1000:]}")
    n=len(hits); nl=len(loci); fraction=n/max(1,total)
    status = mycoplasma_call(n, nl, total, args.minimum_screened_reads)
    result.update(status=status,qualifying_reads=n,qualifying_loci=nl,screened_reads=total,qualifying_fraction=fraction,
                  qualifying_read_loci={name: sorted(values) for name, values in read_loci.items()},
                  interpretation="Opportunistic competitive-alignment screen; a negative result at shallow depth is not a contamination-free certificate.")
    json_dump(result,args.output); print(json.dumps(result,indent=2))


def parser():
    p = argparse.ArgumentParser()
    sub = p.add_subparsers(dest="command", required=True)
    x = sub.add_parser("prepare")
    x.add_argument("--sample", required=True); x.add_argument("--assay", choices=["chromatin", "microc"], required=True)
    x.add_argument("--fastq", nargs="+", required=True); x.add_argument("--inserts", required=True); x.add_argument("--legs", required=True)
    x.add_argument("--leg1", required=True); x.add_argument("--leg2", required=True)
    x.add_argument("--read-table", required=True); x.add_argument("--index-assignments", required=True); x.add_argument("--summary", required=True)
    x.add_argument("--raw-lengths", required=True); x.add_argument("--insert-lengths", required=True)
    x.add_argument("--leg-length", type=int, default=100); x.add_argument("--max-unit", type=int, default=5000)
    x.add_argument("--adapter-dimer-max-bp", type=int, default=10)
    x.add_argument("--ultrashort-insert-max-bp", type=int, default=50)
    x.add_argument("--adapter-dimer-raw-min-bp", type=int, default=120); x.add_argument("--adapter-dimer-raw-max-bp", type=int, default=170)
    x.add_argument("--adapter-dimer-template-identity", type=float, default=.80)
    x.add_argument("--single-index-pair-fraction", type=float, default=.90)
    x.add_argument("--minimum-index-informative-units", type=int, default=20)
    x.add_argument("--minimum-index-cluster-units", type=int, default=20)
    x.add_argument("--minimum-index-cluster-fraction", type=float, default=.01)
    x.add_argument("--index-variant-max-edit", type=int, default=2,
                   help="Maximum total i5+i7 edit distance for ONT index-variant rescue")
    x.add_argument("--index-variant-dominance-ratio", type=float, default=5.0,
                   help="Minimum support ratio of a dominant pair over a rescued variant")
    x.add_argument("--expected-index-pairs", default="")
    x.set_defaults(func=cmd_prepare)
    x = sub.add_parser("metrics")
    x.add_argument("--sample", required=True); x.add_argument("--assay", choices=["chromatin", "microc"], required=True)
    x.add_argument("--bam", required=True); x.add_argument("--samtools", default="samtools"); x.add_argument("--adapter-summary", required=True)
    x.add_argument("--read-table", required=True); x.add_argument("--index-assignments", required=True)
    x.add_argument("--assembly", default="unspecified"); x.add_argument("--mycoplasma", required=True); x.add_argument("--spike-in", required=True)
    x.add_argument("--chrom-sizes", required=True); x.add_argument("--primary", default="auto"); x.add_argument("--mito", default="auto")
    x.add_argument("--aligner", default="unknown")
    x.add_argument("--mapq", type=int, default=metric_profile()["mapq"])
    x.add_argument("--collision-bin-size", type=int, default=metric_profile()["collision_qc_bin_size"])
    x.add_argument("--microc-cis-pass", type=float, default=None, help=argparse.SUPPRESS); x.add_argument("--microc-cis-fail", type=float, default=None, help=argparse.SUPPRESS)
    x.add_argument("--microc-long-range-cis-min-bp", type=int, default=metric_profile()["microc_long_range_cis_min_bp"])
    x.add_argument("--microc-minimum-contacts", type=int, default=metric_profile()["microc_minimum_contacts"])
    x.add_argument("--microc-short-inward-bp", type=int, default=metric_profile()["microc_short_inward_bp"])
    x.add_argument("--microc-minimum-query-coverage", type=float, default=metric_profile()["microc_minimum_query_coverage"])
    x.add_argument("--dimer-pass-max", type=float, default=.05)
    x.add_argument("--dimer-warning-max", type=float, default=.10)
    x.add_argument("--external-peaks", default="")
    x.add_argument("--minimum-raw", type=int, default=metric_profile()["minimum_raw_reads"]); x.add_argument("--minimum-units", type=int, default=metric_profile()["minimum_complete_units"])
    x.add_argument("--ligation", default="", help="Full-insert split-alignment ligation JSON for Micro-C")
    x.add_argument("--metrics", required=True)
    x.add_argument("--fragments", required=True); x.add_argument("--pairs", required=True)
    x.set_defaults(func=cmd_metrics)
    x = sub.add_parser("spikein")
    x.add_argument("--fastq", required=True); x.add_argument("--output", required=True)
    x.add_argument("--enabled", action="store_true"); x.add_argument("--index", default="")
    x.add_argument("--target-regex", default=""); x.add_argument("--label", default="spike-in")
    x.add_argument("--bowtie2", default="bowtie2"); x.add_argument("--mapq", type=int, default=20)
    x.add_argument("--minimum-aligned-bp", type=int, default=80)
    x.add_argument("--minimum-query-coverage", type=float, default=.80)
    x.add_argument("--minimum-identity", type=float, default=.85)
    x.add_argument("--minimum-screened-reads", type=int, default=500)
    x.set_defaults(func=cmd_spikein)
    x = sub.add_parser("mycoplasma")
    x.add_argument("--fastq", required=True); x.add_argument("--output", required=True)
    x.add_argument("--enabled", action="store_true"); x.add_argument("--index", default=""); x.add_argument("--target-regex", default="mycoplasma_target|mycoplasmopsis|mycoplasmoides|metamycoplasma|mesomycoplasma|acholeplasma|spiroplasma|ureaplasma|mollicutes")
    x.add_argument("--bowtie2", default="bowtie2"); x.add_argument("--mapq", type=int, default=20)
    x.add_argument("--minimum-aligned-bp", type=int, default=80); x.add_argument("--minimum-query-coverage", type=float, default=.80); x.add_argument("--minimum-identity", type=float, default=.85)
    x.add_argument("--minimum-screened-reads", type=int, default=500)
    x.set_defaults(func=cmd_mycoplasma)
    x = sub.add_parser("track")
    x.add_argument("--bam", required=True); x.add_argument("--output", required=True)
    x.add_argument("--samtools", default="samtools"); x.add_argument("--converter", default="bedGraphToBigWig")
    x.add_argument("--chrom-sizes", required=True); x.add_argument("--primary", default="auto"); x.add_argument("--mapq", type=int, default=20)
    x.set_defaults(func=cmd_track)
    x = sub.add_parser("report")
    x.add_argument("--reuse-plot-assets", action="store_true", help="Reuse existing karyotype and sublibrary metaplot assets. Remove assets affected by input/config changes first; other plots and tables are refreshed.")
    x.add_argument("--config", required=True); x.add_argument("--out", required=True); x.set_defaults(func=cmd_report)
    return p


if __name__ == "__main__":
    arguments = parser().parse_args()
    arguments.func(arguments)
