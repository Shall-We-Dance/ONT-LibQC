"""Shared adapter sequences and sequence primitives; no pipeline or report imports."""
import gzip

def microc_leg_span(insert_length, maximum):
    """Length of each non-overlapping endpoint used in the Micro-C contact workflow."""
    return min(maximum, insert_length // 2) if insert_length >= 50 else 0


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


P7_FORWARD = rc(P7)
