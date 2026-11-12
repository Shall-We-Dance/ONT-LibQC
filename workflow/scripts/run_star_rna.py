#!/usr/bin/env python3
"""Bounded STARlong alignment of reconstructed ONT RNA library inserts.

The sparse genome is prepared separately. STAR does not sort BAM or insert
annotations at run time; samtools sorts with one compression thread and a fixed
memory limit. Reads and qualities are preserved, including variable lengths.
"""
import argparse
import gzip
import json
from pathlib import Path
import resource
import subprocess
import tempfile
import time


def prepare_fastq(source, target):
    count = maximum = 0
    opener = gzip.open if str(source).endswith('.gz') else open
    with opener(source, 'rt') as src, open(target, 'w') as dst:
        while True:
            name = src.readline()
            if not name:
                break
            seq, plus, qual = src.readline().rstrip(), src.readline(), src.readline().rstrip()
            if not name.startswith('@') or not plus.startswith('+') or len(seq) != len(qual) or not seq:
                raise ValueError('Malformed or empty FASTQ record supplied to STARlong')
            if len(seq) > 500000:
                raise ValueError('STARlong supports inserts up to 500,000 bases')
            dst.write(name.split()[0] + '\n' + seq + '\n+\n' + qual + '\n')
            count += 1
            maximum = max(maximum, len(seq))
    return count, maximum


def star_command(star, index, fastq, prefix, threads, maximum):
    return [star, '--genomeDir', str(index), '--readFilesIn', str(fastq),
            '--runThreadN', str(threads), '--genomeLoad', 'NoSharedMemory',
            '--twopassMode', 'None', '--outFileNamePrefix', str(prefix),
            '--outSAMtype', 'BAM', 'Unsorted', '--outSAMunmapped', 'Within',
            '--outSAMattributes', 'NH', 'HI', 'AS', 'nM', 'NM', 'MD',
            '--outSAMmapqUnique', '60', '--readNameSeparator', ' ',
            '--alignEndsType', 'Local', '--alignIntronMax', '1000000',
            '--alignSJoverhangMin', '8', '--alignSJDBoverhangMin', '3',
            '--outFilterMultimapNmax', '20', '--outFilterMismatchNmax', '2000',
            '--outFilterMismatchNoverLmax', '0.2', '--outFilterMismatchNoverReadLmax', '0.2',
            '--outFilterMatchNmin', '20', '--outFilterMatchNminOverLread', '0.5',
            '--outFilterScoreMinOverLread', '0', '--scoreDelOpen', '-1',
            '--scoreDelBase', '-1', '--scoreInsOpen', '-1', '--scoreInsBase', '-1',
            '--scoreGapNoncan', '-20', '--seedSearchStartLmax', '50',
            '--seedPerReadNmax', '10000', '--seedPerWindowNmax', '1000',
            '--alignTranscriptsPerReadNmax', '10000', '--alignTranscriptsPerWindowNmax', '1000',
            '--limitIObufferSize', '10000000', '20000000',
            '--limitOutSAMoneReadBytes', str(max(100000, 2 * (maximum + 100) * 20))]


def align(args):
    index = Path(args.index).resolve()
    for name in ('Genome', 'SA', 'SAindex', 'genomeParameters.txt', 'chrNameLength.txt'):
        if not (index / name).is_file():
            raise ValueError(f'Incomplete STAR index: missing {name}')
    output, summary = Path(args.output), Path(args.summary)
    output.parent.mkdir(parents=True, exist_ok=True)
    summary.parent.mkdir(parents=True, exist_ok=True)
    started = time.monotonic()
    version = subprocess.check_output([args.star, '--version'], text=True).strip()
    with tempfile.TemporaryDirectory(prefix='star-rna-', dir=output.parent) as temp:
        tmp = Path(temp)
        query = tmp / 'inserts.fastq'
        count, maximum = prepare_fastq(args.fastq, query)
        command = star_command(args.star, index, query, tmp / 'star.', args.threads, maximum)
        if count:
            subprocess.run(command, check=True)
            unsorted = tmp / 'star.Aligned.out.bam'
        else:
            unsorted = tmp / 'empty.sam'
            with unsorted.open('w') as handle:
                handle.write('@HD\tVN:1.6\tSO:unsorted\n')
                for row in (index / 'chrNameLength.txt').read_text().splitlines():
                    name, size = row.split()[:2]
                    handle.write(f'@SQ\tSN:{name}\tLN:{size}\n')
        subprocess.run([args.samtools, 'sort', '-@', '1', '-m', '256M',
                        '-T', str(tmp / 'sort'), '-o', str(output), str(unsorted)], check=True)
        metrics = {}
        log = tmp / 'star.Log.final.out'
        if log.exists():
            for line in log.read_text().splitlines():
                if '|' in line:
                    key, value = line.split('|', 1)
                    metrics[key.strip()] = value.strip()
        marker = index / 'ontlibqc-star-index.json'
        evidence = {'schema_version': 1, 'status': 'BETA', 'aligner': 'STARlong', 'version': version,
                    'input_units': count, 'maximum_insert_bp': maximum,
                    'parameters': command[1:], 'star_metrics': metrics,
                    'elapsed_seconds': round(time.monotonic() - started, 3),
                    'peak_child_rss_kib': resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss,
                    'sort_memory': '256M', 'sort_threads': 1,
                    'reference_index': json.loads(marker.read_text()) if marker.exists() else {'status': 'operator_supplied'},
                    'interpretation': 'Splice-aware beta tuned for reconstructed, variable-length ONT inserts. Annotation-supported junctions are not independent validation; no full-length transcript or base-resolution crosslink claim.'}
        summary.write_text(json.dumps(evidence, indent=2) + '\n')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--fastq', required=True)
    parser.add_argument('--index', required=True)
    parser.add_argument('--output', required=True)
    parser.add_argument('--summary', required=True)
    parser.add_argument('--star', default='STARlong')
    parser.add_argument('--samtools', default='samtools')
    parser.add_argument('--threads', type=int, default=2)
    args = parser.parse_args()
    if not 1 <= args.threads <= 64:
        parser.error('threads must be in 1..64')
    align(args)


if __name__ == '__main__':
    main()
