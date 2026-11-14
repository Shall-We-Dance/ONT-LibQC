#!/usr/bin/env python3
"""Chromatin local alignment with conservative terminal foldback rescue.

ZQ:i:0 marks alignments failing the length/coverage eligibility filter. They
remain in the BAM for mapping-rate accounting but never contribute to QC tracks.
ZF:i:1 marks sequence-supported terminal correction verified by end-to-end remap.
"""
from pathlib import Path
import argparse, collections, gzip, json, math, re, subprocess, tempfile
from ont_qc import edit_distance, rc, fastq_records

METHOD = "local-foldback-v1"
REVISION = "r2"

def read_sam(path):
    records, headers = {}, []
    with path.open() as handle:
        for line in handle:
            if line.startswith('@'):
                headers.append(line)
                continue
            fields = line.rstrip().split('\t')
            flag = int(fields[1])
            if flag & (4 | 256 | 512 | 1024 | 2048):
                continue
            ops = [(int(n), op) for n, op in re.findall(r'(\d+)([MIDNSHP=X])', fields[5])]
            left = ops[0][0] if ops[0][1] == 'S' else 0
            right = ops[-1][0] if ops[-1][1] == 'S' else 0
            if flag & 16:
                left, right = right, left
            aligned = sum(n for n, op in ops if op in 'M=X')
            query_span = sum(n for n, op in ops if op in 'MI=X')
            reference_span = sum(n for n, op in ops if op in 'MDN=X')
            span = sum(n for n, op in ops if op in 'MID=X')
            tags = {v[:2]: v[5:] for v in fields[11:]}
            records[fields[0]] = {
                'line': line, 'fields': fields, 'chrom': fields[2],
                'start': int(fields[3])-1, 'end': int(fields[3])-1+reference_span,
                'strand': '-' if flag & 16 else '+', 'mapq': int(fields[4]),
                'left': left, 'right': right, 'aligned': aligned,
                'coverage': query_span / len(fields[9]),
                'identity': 1-int(tags.get('NM', '0')) / span, 'seq': fields[9],
            }
    return headers, records


def entropy(s):
 counts=collections.Counter(s)
 return -sum((v/len(s))*math.log2(v/len(s)) for v in counts.values())
def foldback_evidence(clip,mid,side):
 if len(clip)<14 or entropy(clip)<1.5 or set(clip)-set('ACGT'):return None
 errors=len(clip)//10;matches=[]
 for offset in range(4):
  for size in range(max(14,len(clip)-errors),len(clip)+errors+1):
   if offset+size>len(mid):continue
   # The observed structure mirrors the opposite insert end across the insert body.
   piece=mid[len(mid)-offset-size:len(mid)-offset if offset else None] if side=='left' else mid[offset:offset+size]
   dist=edit_distance(clip,rc(piece),errors)
   if dist<=errors:matches.append((dist,offset,size))
 if not matches:return None
 dist,offset,size=min(matches)
 return {'terminal_bp':len(clip),'edit_distance':dist,'opposite_end_offset_bp':offset,'matched_bp':size,'sequence_identity':1-dist/max(size,len(clip)),'terminal_entropy':entropy(clip)}

def eligible(record):
    return record['aligned'] >= 30 and record['coverage'] >= .8


def collect_candidates(records, inputs, chromosomes, mapq):
    candidates = {}
    for uid, r in records.items():
        if not (mapq <= r['mapq'] < 255 and r['chrom'] in chromosomes and r['aligned'] >= 30):
            continue
        left, right = r['left'], r['right']
        seq, quality = inputs[uid]
        side = 'left' if left >= right else 'right'
        n = max(left, right)
        if min(left, right) > 3 or n < 14 or n / len(seq) > .4:
            continue
        clip = seq[:n] if side == 'left' else seq[-n:]
        middle = seq[left:len(seq)-right if right else None]
        evidence = foldback_evidence(clip, middle, side)
        if evidence:
            candidates[uid] = dict(evidence, side=side, original_bp=len(seq))
    return candidates


def verified_remap(prior, corrected, mapq):
    return (mapq <= corrected['mapq'] < 255 and eligible(corrected)
            and corrected['identity'] >= .9
            and (prior['chrom'], prior['strand']) == (corrected['chrom'], corrected['strand'])
            and max(abs(prior['start']-corrected['start']), abs(prior['end']-corrected['end'])) <= 5)


def run(args):
    inputs = {name: (seq, qual) for name, seq, qual in fastq_records([args.fastq])}
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix='chromatin-local-', dir=Path(args.output).parent) as directory:
        tmp = Path(directory)
        local_sam = Path(args.local_sam) if args.local_sam else tmp/'local.sam'
        if not args.local_sam:
            subprocess.run([args.bowtie2, '--very-sensitive-local', '--no-unal', '-p', str(args.threads),
                            '-x', args.index, '-U', args.fastq, '-S', str(local_sam)], check=True)
        headers, records = read_sam(local_sam)
        candidates = collect_candidates(records, inputs, set(args.chromosomes.split(',')), args.mapq)
        accepted = {}
        if candidates:
            fastq = tmp/'corrected.fastq.gz'
            with gzip.open(fastq, 'wt') as handle:
                for uid, evidence in candidates.items():
                    seq, qual = inputs[uid]
                    n = evidence['terminal_bp']
                    seq, qual = (seq[n:], qual[n:]) if evidence['side'] == 'left' else (seq[:-n], qual[:-n])
                    handle.write(f'@{uid}\n{seq}\n+\n{qual}\n')
            corrected_sam = tmp/'corrected.sam'
            subprocess.run([args.bowtie2, '--end-to-end', '--very-sensitive', '--no-unal', '-p', str(args.threads),
                            '-x', args.index, '-U', str(fastq), '-S', str(corrected_sam)], check=True)
            _, corrected = read_sam(corrected_sam)
            accepted = {uid: r for uid, r in corrected.items() if verified_remap(records[uid], r, args.mapq)}
        newly_eligible = sum(not eligible(records[uid]) for uid in accepted)
        with Path(args.output).open('w') as handle:
            handle.writelines(headers)
            handle.write('@CO\tONT-LibQC chromatin local-foldback-v1; MAPQ unchanged; ZQ=0 excludes insufficient aligned length/coverage\n')
            for uid, original in records.items():
                r = accepted.get(uid, original)
                fields = [v for i, v in enumerate(r['fields']) if i < 11 or v[:3] not in {'ZQ:', 'ZF:'}]
                fields.append('ZQ:i:' + str(int(eligible(r))))
                if uid in accepted:
                    fields.append('ZF:i:1')
                handle.write('\t'.join(fields)+'\n')
        summary = dict(method=METHOD, analysis_revision=REVISION, input_units=len(inputs),
                       local_mapped=len(records), mapq=args.mapq, minimum_aligned_bp=30,
                       minimum_query_coverage=.8, foldback_candidates=len(candidates),
                       foldback_verified=len(accepted), foldback_newly_eligible=newly_eligible,
                       foldback_minimum_bp=14, foldback_minimum_identity=.9,
                       foldback_minimum_entropy=1.5, foldback_maximum_fraction=.4,
                       foldback_maximum_other_end_clip_bp=3, foldback_maximum_opposite_end_offset_bp=3,
                       remap_minimum_identity=.9, remap_maximum_locus_shift_bp=5,
                       foldback_model='terminal reverse complement of opposite aligned insert end')
        Path(args.summary).parent.mkdir(parents=True, exist_ok=True)
        Path(args.summary).write_text(json.dumps(summary, indent=2)+'\n')
        return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ['fastq', 'index', 'output', 'summary', 'chromosomes']:
        parser.add_argument('--'+name, required=True)
    parser.add_argument('--bowtie2', default='bowtie2')
    parser.add_argument('--threads', type=int, default=4)
    parser.add_argument('--mapq', type=int, default=20)
    parser.add_argument('--local-sam', help='Reuse an alignment from the same input/reference for a controlled reanalysis')
    args = parser.parse_args()
    if args.mapq < 0 or args.mapq >= 255 or args.threads < 1:
        parser.error('Invalid MAPQ or thread count')
    run(args)


if __name__ == '__main__':
    main()
