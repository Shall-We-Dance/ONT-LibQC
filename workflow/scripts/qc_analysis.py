"""Analysis commands and numerical helpers, independent of report assembly."""
from profile_data import save_profile
import bisect
import collections
import gzip
import html
import json
import math
import statistics
import subprocess
import tempfile
from pathlib import Path
from qc_definitions import profile as metric_profile, metric_settings, summarize_contacts, contact_evidence_status
from report_palette import scope_color
from qc_plot_metrics import normalized_collision_stats, collision_metric_fields, svg_size_trace, size_trace_data
from spikein_auto import scoped_summary as auto_spikein_summary
from library_structure import normalize_structure, umi_lengths, handles, extract_umis, adapter_families, POSITIONS as UMI_POSITIONS
from report_mitochondrial import coverage_html as mitochondrial_coverage_html, integrate_karyotype
from report_provenance import capture_analysis
from library_sequences import microc_leg_span, rc, P5, TRU_R1, FAMILIES, open_text, fastq_records, edit_distance

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

def find_units(seq, max_unit, structure=None):
    all_units = []
    matches_by_family = {}
    for family, (r1, r2) in adapter_families(structure, FAMILIES).items():
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


def single_index_p5_hit(seq, unit, structure):
    """Recognize the overlapping P5/Read-1 junction of TruSeq LT libraries.

    Illumina's single-index universal adapter shares ACAC between P5 and
    Read-1; dual-index adapters instead contain P5 + i5 + Read-1. Missing P5
    alone is never evidence for an absent i5.
    https://support-docs.illumina.com/SHARE/AdapterSequences/Content/SingleIndexes.htm
    """
    if (unit['family'] != 'truseq' or structure['p5_sequence']
            or structure['adapter_i5_sequence'] or structure['i5_length'] not in (None, 0)
            or any(u['position'] in ('before_i5', 'after_i5') for u in structure['umis'])):
        return None
    reverse = unit['orientation'] == 'p7_to_p5'
    oriented = rc(seq) if reverse else seq
    anchor = len(seq) - unit['right']['end'] if reverse else unit['left']['start']
    prefix = P5[:-4]
    lo, hi = max(0, anchor - len(prefix) - 3), anchor + len(TRU_R1) + 3
    matches = approximate_matches(oriented[lo:hi], prefix + TRU_R1, .12)
    for hit in sorted(matches, key=lambda h: h['edits']):
        start, end = lo + hit['start'], lo + hit['end']
        if (abs(anchor - start - len(prefix)) > 2
                or abs(end - anchor - len(TRU_R1)) > 3
                or edit_distance(prefix, oriented[start:anchor], 3) > 3):
            continue
        # Save the physical overlap explicitly for the library-length model.
        a, b = start, anchor + 4
        return dict(start=len(seq)-b if reverse else a,
                    end=len(seq)-a if reverse else b,
                    edits=hit['edits'], identity=hit['identity'])
    return None


def add_indexes(seq, units, structure=None):
    structure = structure or normalize_structure()
    p5, p7_forward = handles(structure)
    p7 = rc(p7_forward)
    lengths = umi_lengths(structure)
    def gap(index):
        if structure['mode'] == 'auto':
            return 18
        return (structure[index + '_length'] if structure[index + '_length'] is not None else 24) + 2 + lengths.get('before_' + index, 0) + lengths.get('after_' + index, 0)
    gap5, gap7 = gap('i5'), gap('i7')
    for unit in units:
        left, right = unit["left"], unit["right"]
        i5_obs = i7_obs = None
        p5_hit = p7_hit = None
        if unit["orientation"] == "p5_to_p7":
            p5_hit = nearest_before(matches_in_window(seq, p5, left["start"] - len(p5) - gap5 - 2, left["start"]), left["start"], max_gap=gap5)
            p7_hit = nearest_after(matches_in_window(seq, p7_forward, right["end"], right["end"] + len(p7) + gap7 + 2), right["end"], max_gap=gap7)
            if p5_hit:
                i5_obs = seq[p5_hit["end"]:left["start"]]
            if p7_hit:
                i7_obs = seq[right["end"]:p7_hit["start"]]
            i5 = i5_obs
            i7 = i7_obs
        else:
            p7_hit = nearest_before(matches_in_window(seq, p7, left["start"] - len(p7) - gap7 - 2, left["start"]), left["start"], max_gap=gap7)
            p5_hit = nearest_after(matches_in_window(seq, rc(p5), right["end"], right["end"] + len(p5) + gap5 + 2), right["end"], max_gap=gap5)
            if p7_hit:
                i7_obs = seq[p7_hit["end"]:left["start"]]
            if p5_hit:
                i5_obs = seq[right["end"]:p5_hit["start"]]
            i5 = rc(i5_obs) if i5_obs is not None else None
            i7 = rc(i7_obs) if i7_obs is not None else None
        overlap = 0 if p5_hit is not None else None
        if p5_hit is None:
            single_hit = single_index_p5_hit(seq, unit, structure)
            if single_hit is not None:
                p5_hit = single_hit
                i5 = i5_obs = ""
                overlap = 4
        unit.update({"p5_adapter_overlap": overlap,
                     "p5_hit": p5_hit, "p7_hit": p7_hit,
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


def index_consensus(values, maximum_length=14):
    values = [v for v in values if v is not None and len(v) <= maximum_length]
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
    if (observed == "") != (expected == ""):
        return False
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
        if (left[key] == "") != (right[key] == ""):
            return maximum_edit + 1
        remaining = maximum_edit - total
        total += edit_distance(left[key], right[key], remaining)
        if total > maximum_edit:
            break
    return total


def cluster_index_pairs(units, single_fraction=.90, minimum_informative=20,
                        minimum_cluster_units=20, minimum_cluster_fraction=.01,
                        expected_text="", variant_max_edit=2,
                        variant_dominance_ratio=5.0, maximum_index_length=14):
    informative = [u for u in units if u["p5_detected"] and u["p7_detected"] and u["i5"] is not None and u["i7"] is not None]
    exact = collections.Counter((u["i5"], u["i7"]) for u in informative
                                if u["i5"] is not None and u["i7"] is not None
                                and len(u["i5"]) <= maximum_index_length and len(u["i7"]) <= maximum_index_length)
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
            observed = {key: unit[key] for key, handle in (("i5", "p5_detected"), ("i7", "p7_detected"))
                        if unit[handle] and unit[key] is not None}
            # Neither unknown ends nor an unindexed end identify a sublibrary.
            if not any(value and seed[key] for key, value in observed.items()):
                if not (seed['i5'] == seed['i7'] == '' and len(observed) == 2):
                    continue
            if all(index_close(value, seed[key], variant_max_edit) for key, value in observed.items()):
                distance = sum(edit_distance(value, seed[key], variant_max_edit)
                               for key, value in observed.items())
                candidates.append((distance, seed))
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
            assigned["informative_support"] += int(unit["p5_detected"] and unit["p7_detected"]
                                                   and unit["i5"] is not None and unit["i7"] is not None)
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
    by_id = {u['unit_id']: u for u in units}
    for seed in seeds:
        overlaps = collections.Counter(by_id[uid].get('p5_adapter_overlap') for uid in seed['unit_ids']
                                       if by_id[uid].get('p5_adapter_overlap') is not None)
        public_pairs.append({
            "index_architecture": index_architecture(seed['i5'], seed['i7']),
            "p5_adapter_overlap": overlaps.most_common(1)[0][0] if overlaps else 0,
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

def index_architecture(i5, i7):
    if i5 is None or i7 is None:
        return 'partial_or_custom'
    if i5 == '' and i7:
        return 'classic_single_index'
    if i7 == '' and i5:
        return 'i5_single_index'
    return 'physical_dual_index' if i5 and i7 else 'unindexed'


def zero_insert_template(family, i5, i7, structure=None, p5_adapter_overlap=0):
    r1, r2 = adapter_families(structure, FAMILIES)[family]
    p5, p7 = handles(structure or normalize_structure())
    if p5_adapter_overlap == 4 and family == 'truseq' and i5 == '' and p5 == P5:
        p5 = p5[:-4]
    return p5 + (i5 or "") + r1 + rc(r2) + (i7 or "") + p7


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
    structure = normalize_structure(json.loads(getattr(args, 'library_structure', '{}')))
    umi_sizes = umi_lengths(structure)
    inner_umi_bases = umi_sizes.get('before_insert', 0) + umi_sizes.get('after_insert', 0)
    structure_rejected = 0
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
        units, _ = find_units(seq, args.max_unit + inner_umi_bases, structure)
        add_indexes(seq, units, structure)
        retained = []
        for unit in units:
            if extract_umis(seq, unit, structure, rc):
                retained.append(unit)
            else:
                structure_rejected += 1
        units = retained
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
                "p5_adapter_overlap": unit["p5_adapter_overlap"],
                "p5_detected": unit["p5_hit"] is not None, "p7_detected": unit["p7_hit"] is not None,
                "is_adapter_dimer": canonical_dimer, "dimer_class": "canonical_adapter_dimer" if canonical_dimer else "not_dimer",
                "template_identity": None, "umis": unit["umis"],
                "umi_sequence": unit["umi_sequence"], "umi_status": unit["umi_status"],
            })
    pairing, assignments = cluster_index_pairs(
        units_all, args.single_index_pair_fraction, args.minimum_index_informative_units,
        args.minimum_index_cluster_units, args.minimum_index_cluster_fraction, args.expected_index_pairs,
        args.index_variant_max_edit, args.index_variant_dominance_ratio,
        maximum_index_length=24 if structure["mode"] == "custom" else 14)
    unit_scope = {row["unit_id"]: row["scope_id"] for row in assignments if row["scope_id"]}
    pair_by_scope = {pair["scope_id"]: pair for pair in pairing["pairs"]}
    raw_dimer_reads = set()
    for unit in units_all:
        scope_id = unit_scope.get(unit["unit_id"])
        pair = pair_by_scope.get(scope_id)
        if pair and not umi_sizes and args.adapter_dimer_raw_min_bp <= unit["raw_length"] <= args.adapter_dimer_raw_max_bp:
            template = zero_insert_template(unit["family"], pair["i5"], pair["i7"], structure, pair.get("p5_adapter_overlap", 0))
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
    family = max(families, key=families.get) if families else next(iter(adapter_families(structure, FAMILIES)))
    for pair in pairing["pairs"] if not umi_sizes else []:
        templates.append((pair["scope_id"], zero_insert_template(family, pair["i5"], pair["i7"], structure, pair.get("p5_adapter_overlap", 0))))
    unreconstructed_probable = None if umi_sizes else 0
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
                    "dimer_class", "dimer_template_identity", "p5_detected", "p7_detected", "p5_adapter_overlap",
                    "i5_observed", "i5_normalized", "i7_observed", "i7_normalized",
                    "umi_status", "umi_sequence"] + ["umi_" + position for position in UMI_POSITIONS]
    with gzip.open(args.inserts, "wt") as insert_out, gzip.open(args.legs, "wt") as leg_out, \
            gzip.open(args.leg1, "wt") as leg1_out, gzip.open(args.leg2, "wt") as leg2_out, \
            gzip.open(args.read_table, "wt") as table:
        table.write("\t".join(table_header) + "\n")
        for unit in units_all:
            assignment = assignment_by_unit[unit["unit_id"]]
            if not unit["is_adapter_dimer"]:
                insert_out.write(f"@{unit['unit_id']}\n{unit['insert_seq']}\n+\n{unit['insert_qual']}\n")
            leg_len = microc_leg_span(unit["insert_length"], args.leg_length) if args.assay == "microc" else 0
            if leg_len and not unit["is_adapter_dimer"]:
                left_seq, left_q = unit["insert_seq"][:leg_len], unit["insert_qual"][:leg_len]
                right_seq, right_q = rc(unit["insert_seq"][-leg_len:]), unit["insert_qual"][-leg_len:][::-1]
                leg_out.write(f"@{unit['unit_id']}/1\n{left_seq}\n+\n{left_q}\n@{unit['unit_id']}/2\n{right_seq}\n+\n{right_q}\n")
                leg1_out.write(f"@{unit['unit_id']}\n{left_seq}\n+\n{left_q}\n")
                leg2_out.write(f"@{unit['unit_id']}\n{right_seq}\n+\n{right_q}\n")
            row = [unit["read_id"], unit["raw_length"], f"{unit['mean_q']:.3f}", unit["unit_id"],
                   assignment["scope_id"], assignment["assignment_status"],
                   unit["family"], unit["orientation"], unit["insert_start"], unit["insert_end"],
                   unit["insert_length"], unit["adapter_edits"], int(unit["is_adapter_dimer"]),
                   unit["dimer_class"], unit["template_identity"], int(unit["p5_detected"]), int(unit["p7_detected"]), unit["p5_adapter_overlap"],
                   unit["i5_observed"], unit["i5"], unit["i7_observed"], unit["i7"],
                   unit["umi_status"], unit["umi_sequence"]] + [unit["umis"].get(position) for position in UMI_POSITIONS]
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
    consensus_limit = 24 if structure["mode"] == "custom" else 14
    i5, i7 = index_consensus(i5_values, consensus_limit), index_consensus(i7_values, consensus_limit)
    architectures = {pair['index_architecture'] for pair in pairing['pairs']}
    architecture = (next(iter(architectures)) if len(architectures) == 1 else
                    'mixed_index_architecture' if architectures else 'partial_or_custom')
    raw_window_count = sum(args.adapter_dimer_raw_min_bp <= x <= args.adapter_dimer_raw_max_bp for x in raw_lengths)
    summary = {
        "sample": args.sample, "assay": args.assay, "raw": length_summary(raw_lengths),
        "library_structure": structure,
        "umi_summary": {"configured_bases": sum(umi_sizes.values()),
                        "complete_units": sum(u['umi_status'] == 'complete' for u in units_all),
                        "partial_units": sum(u['umi_status'] == 'partial' for u in units_all),
                        "rejected_short_units": structure_rejected,
                        "orientation": "P5-to-P7", "deduplication": "unchanged"},
        "unreconstructed_dimer_screen": "not_assessed_for_umi_layout" if umi_sizes else "assessed",
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


def sam_primary_records(bam, samtools, include_ineligible=False):
    process = subprocess.Popen([samtools, "view", bam], text=True, stdout=subprocess.PIPE)
    best = {}
    for line in process.stdout:
        f = line.rstrip().split("\t")
        flag = int(f[1])
        if flag & (4 | 256 | 512 | 1024 | 2048):
            continue
        if not include_ineligible and "ZQ:i:0" in f[11:]:
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



def write_karyoploter(records, chrom_sizes, chromosome_names, mapq, bins_per_longest,
                      rscript, script, genes_bed, output, title, assembly, gene_window_bp=250000, mitochondrial=()):
    """Render nuclear chromosomes and append a matching, independently scaled chrM row."""
    names = list(dict.fromkeys(name for name in chromosome_names if name in chrom_sizes and name not in mitochondrial))
    if not names:
        raise ValueError("Karyotype plotting requires at least one configured chromosome")
    longest = max(chrom_sizes[name] for name in names)
    window_size = max(1, math.ceil(longest / max(1, bins_per_longest)))
    eligible = []
    name_set = set(names)
    for record in records.values():
        if record["chrom"] not in name_set or not mapq <= record["mapq"] < 255:
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
    counts = {name: [0] * math.ceil(chrom_sizes[name] / window_size) for name in names}
    for chromosome, position in eligible:
        counts[chromosome][(position - 1) // window_size] += 1
    metadata = {"coordinate_scale": "shared", "density_scale": "shared", "window_size_bp": window_size,
                "chromosome_lengths": {name: chrom_sizes[name] for name in names},
                "midpoint_counts": counts, "mapq": mapq}
    svg = Path(output).read_text()
    svg = svg.replace('<svg', '<svg data-chromosome-scales="shared"', 1)
    svg = svg.replace('</svg>', '<metadata id="chromosome-coverage">' + html.escape(json.dumps(metadata)) + '</metadata></svg>')
    if mitochondrial:
        panel = mitochondrial_coverage_html(records, chrom_sizes, mitochondrial, mapq, primary=names)
        svg = integrate_karyotype(svg, panel)
    Path(output).write_text(svg)
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
    structure = normalize_structure(adapter_summary.get('library_structure'))
    p5, p7 = handles(structure)
    umi_bases = sum(umi_lengths(structure).values())
    templates = adapter_families(structure, FAMILIES)
    pairs = {p['scope_id']: p for p in adapter_summary.get('index_pairing', {}).get('pairs', [])}
    modeled = []
    excluded = collections.Counter()
    for row in unit_rows:
        family = row.get('family')
        if family not in templates:
            excluded['unknown_adapter_family'] += 1
            continue
        pair = pairs.get(row.get('scope_id')) if row.get('assignment_status') == 'assigned' else None
        if pair is not None and all(isinstance(pair.get(k), str) for k in ('i5', 'i7')):
            index_lengths = [len(pair[k]) for k in ('i5', 'i7')]
            basis = 'assigned_index_pair_consensus'
        else:
            index_lengths = [adapter_summary.get(k, {}).get('length') for k in ('i5', 'i7')]
            basis = 'sample_index_length_consensus'
        if any(structure[k + '_length'] is not None for k in ('i5', 'i7')):
            index_lengths = [structure[k + '_length'] if structure[k + '_length'] is not None else value
                             for k, value in zip(('i5', 'i7'), index_lengths)]
            basis += '_with_configured_lengths'
        if any(isinstance(v, bool) or not isinstance(v, int) or v < 0 for v in index_lengths):
            excluded['unknown_index_length'] += 1
            continue
        insert = int(row['insert_length'])
        if insert < 0:
            raise ValueError('Full-library model requires non-negative insert lengths')
        overlap = row.get('p5_adapter_overlap')
        if overlap in (None, ''):
            overlap = (pair or {}).get('p5_adapter_overlap', 0)
        overlap = int(overlap or 0)
        added = len(p5) + len(p7) + sum(map(len, templates[family])) + sum(index_lengths) + umi_bases - overlap
        modeled.append(dict(unit_id=row['unit_id'], read_id=row['read_id'],
            insert_length_bp=insert, model_adapter_bases=added, model_umi_bases=umi_bases,
            model_library_length_bp=insert + added, model_basis=basis,
            outer_handles_detected=int(str(row.get('p5_detected')) == '1' and str(row.get('p7_detected')) == '1')))
    return dict(source_units=len(unit_rows), modeled_units=len(modeled),
                excluded_units=sum(excluded.values()), exclusion_reasons=dict(excluded),
                adapter_offsets_bp=sorted({r['model_adapter_bases'] for r in modeled}),
                model_basis_counts=dict(collections.Counter(r['model_basis'] for r in modeled)),
                outer_incomplete_modeled_units=sum(not r['outer_handles_detected'] for r in modeled),
                lengths=[r['model_library_length_bp'] for r in modeled], rows=modeled)


def write_library_size_assets(directory, prefix, label, unit_rows, adapter_summary,
                              maximum_bp=2000, bandwidth=.045, color=None):
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
                             maximum_bp, bandwidth, weighting=weighting, color=color)
        (directory / names[name]).write_text(svg)
    with (directory / names['profiles']).open('w') as handle:
        writer = csv.writer(handle, delimiter='\t')
        writer.writerow(['size_bp', 'length_weighted_relative_signal', 'count_based_relative_signal'])
        writer.writerows(zip(length_data['lengths'], length_data['signal'], count_data['signal']))
    fields = ['unit_id', 'read_id', 'insert_length_bp', 'model_adapter_bases', 'model_umi_bases',
              'model_library_length_bp', 'model_basis', 'outer_handles_detected']
    with (directory / names['units']).open('w') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter='\t')
        writer.writeheader()
        writer.writerows(model['rows'])
    info = {k: v for k, v in model.items() if k not in ('rows', 'lengths')}
    info.update(length_definition='insert + P5 + family-specific inner adapters + P7 + physical i5/i7 lengths + configured UMI bases',
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


def rna_mapping_metrics(records, denominator, mapq, mito):
    """Technical alignment counts only; biological RNA evidence is computed separately."""
    usable = [r for r in records.values() if mapq <= r["mapq"] < 255]
    mitochondrial = sum(r["chrom"] in mito for r in usable)
    return dict(usable_fragments=len(usable),
                usable_fragment_rate=len(usable) / denominator if denominator else None,
                mitochondrial_fraction=mitochondrial / len(usable) if usable else None,
                biological_qc_call="beta_descriptive",
                classification_note="RNA/CLIP beta: source composition, strand orientation, splice evidence and transcript profiles are descriptive. Library-selection methods are inferred only when supported and cannot be uniquely identified from shallow reads. No PCR deduplication or IP-specific enrichment is inferred.")


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
    elif args.assay in {"rnaseq", "clip"}:
        result.update(rna_mapping_metrics(records, mappable, args.mapq, mito))
        result["sequencing_yield_call"] = "pass" if assigned_raw_reads >= args.minimum_raw and units >= args.minimum_units else "inconclusive_low_yield"
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
    all_mapped = sam_primary_records(args.bam, args.samtools, include_ineligible=True) if args.assay == "chromatin" else records
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
        "library_structure": adapter.get("library_structure", {"mode": "auto"}),
        "umi_summary": adapter.get("umi_summary", {}),
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
        "mapped_primary": len(all_mapped),
        "primary_mapping_rate": len(all_mapped) / (mappable_units * (2 if args.assay == "microc" else 1)) if mappable_units else None,
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
    elif args.assay in {"rnaseq", "clip"}:
        metrics = dict(common, **rna_mapping_metrics(records, mappable_units, args.mapq, mito),
                       sequencing_yield_call="pass" if raw_count >= args.minimum_raw and units >= args.minimum_units else "inconclusive_low_yield")
        with gzip.open(args.fragments, "wt") as out:
            out.write("# RNA/CLIP: aligned-block evidence is retained in the RNA QC tables\n")
        with gzip.open(args.pairs, "wt") as out:
            out.write("# RNA/CLIP sample; no contact pairs\n")
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
    metrics["analysis_provenance"] = capture_analysis()
    alignment_summary = Path(args.metrics).parent.parent / "chromatin_alignment" / f"{args.sample}.json"
    if args.assay == "chromatin" and alignment_summary.is_file():
        alignment_method = json.loads(alignment_summary.read_text())
        metrics["alignment_method"] = alignment_method
        metrics["analysis_provenance"]["analysis_revision"] = alignment_method["analysis_revision"]
        metrics["metric_settings"].update(chromatin_alignment=alignment_method["method"],
                                          minimum_aligned_bp=30, minimum_query_coverage=.8)
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
            if args.assay == "chromatin":
                mapped_count = sum(base_unit_id(qname) in allowed for qname in all_mapped)
                scoped["mapped_primary"] = mapped_count
                denominator = scoped["mappable_library_units"]
                scoped["primary_mapping_rate"] = mapped_count / denominator if denominator else None
                if "alignment_method" in metrics:
                    scoped["metric_settings"].update(chromatin_alignment=metrics["alignment_method"]["method"],
                                                      minimum_aligned_bp=30, minimum_query_coverage=.8)
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
        if len(fields) < 6 or "ZQ:i:0" in fields[11:]:
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
    data_path.with_suffix(".json").unlink(missing_ok=True)
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
                        "--outFileNameData", str(data_path), "--plotType", "se", "--colors", scope_color(scope_id),
                        "--startLabel", "TSS", "--endLabel", "TES",
                        "--yAxisLabel", "CPM-normalized signal", "--samplesLabel", label], check=True)
        save_profile(matrix_path, data_path, data_path.with_suffix(".json"))
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


