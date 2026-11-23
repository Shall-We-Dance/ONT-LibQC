"""Executed only inside an offline, single-job container. No credentials are mounted."""
import gzip, json, os, re, shutil, subprocess, sys, time
from pathlib import Path
import yaml

class InputError(ValueError): pass

def normalize_fastq(source,dest,compressed,max_bytes,max_reads):
    total=0;count=0;seen=set(); max_line=2*1024**2
    opener=gzip.open if compressed else open
    try:
        with opener(source,'rb') as src, gzip.open(dest,'wb',compresslevel=1) as out:
            while True:
                head=src.readline(4097)
                if not head: break
                lines=[head]+[src.readline(max_line+1) for _ in range(3)]
                total+=sum(map(len,lines));count+=1
                if total>max_bytes: raise InputError(f'The expanded FASTQ exceeds {max_bytes/1_000_000:g} MB. Run larger files locally with the pipeline on GitHub.')
                if count>max_reads: raise InputError(f'The FASTQ exceeds {max_reads:,} reads. Run larger files locally with the pipeline on GitHub.')
                if any(not line.endswith(b'\n') or len(line)>max_line for line in lines) or len(head)>4096:
                    raise InputError(f'Invalid or excessively long FASTQ record at read {count}. Use four-line FASTQ.')
                h,seq,plus,qual=[line.rstrip(b'\r\n') for line in lines]
                # Dorado FASTQ metadata can follow the read name after a tab.
                rid=re.split(rb'[ \t]',h[1:],maxsplit=1)[0]
                if not h.startswith(b'@') or not re.fullmatch(rb'[A-Za-z0-9_.:=+-]{1,255}',rid):
                    raise InputError(f'Invalid FASTQ identifier at read {count}. Read names must be 1-255 letters, numbers, or . _ : = + -, followed by optional space- or tab-separated metadata.')
                if rid in seen: raise InputError(f'Duplicate read identifier at read {count}.')
                seen.add(rid)
                if not plus.startswith(b'+') or not re.fullmatch(rb'[ACGTRYSWKMBDHVNacgtryswkmbdhvn]+',seq) or len(seq)!=len(qual) or any(q<33 or q>126 for q in qual):
                    raise InputError(f'Invalid sequence or quality values at read {count}.')
                out.write(b'@'+rid+b'\n'+seq.upper()+b'\n+\n'+qual+b'\n')
    except (OSError,EOFError) as e: raise InputError('The FASTQ is unreadable or the gzip stream is incomplete.') from e
    if not count: raise InputError('The FASTQ contains no reads.')
    return count

def normalize_bed(source,dest,sizes):
    chromosomes={f[0]:int(f[1]) for line in Path(sizes).read_text().splitlines() if (f:=line.split())}
    rows=[]
    try:
        for n,line in enumerate(Path(source).read_text(encoding='ascii').splitlines(),1):
            if not line.strip() or line.startswith(('#','track','browser')):continue
            f=line.split()
            if len(f)<3: raise InputError(f'BED line {n} has an unknown chromosome or missing columns.')
            if f[0] not in chromosomes:
                alternate=f[0][3:] if f[0].startswith('chr') else 'chr'+f[0]
                if alternate not in chromosomes: raise InputError(f'BED line {n} has an unknown chromosome or missing columns.')
                f[0]=alternate
            start,end=int(f[1]),int(f[2])
            if not 0<=start<end<=chromosomes[f[0]]: raise InputError(f'BED line {n} is outside the selected reference coordinates.')
            strand=f[5] if len(f)>5 and f[5] in {'+','-'} else '+'
            rows.append(f'{f[0]}\t{start}\t{end}\tregion_{len(rows)+1}\t0\t{strand}\n')
            if len(rows)>100000:raise InputError('BED exceeds 100,000 regions.')
    except (UnicodeError,ValueError) as e:
        if isinstance(e,InputError):raise
        raise InputError('BED must contain plain text chromosome, start and end coordinates.') from e
    if not rows:raise InputError('BED contains no usable regions.')
    Path(dest).write_text(''.join(rows))

def phase_targets(phase, cfg, options):
    sample = options['sample']
    out = cfg['output_dir']
    if phase == 'prepare':
        return [f'{out}/qc/adapters/{sample}.summary.json']
    if phase == 'align':
        targets = [f'{out}/bam/{sample}.bam.bai', f'{out}/qc/mycoplasma/{sample}.json',
                   f'{out}/qc/spike_in/{sample}.json']
        if cfg.get('species_composition', {}).get('indexes'):
            targets.append(f'{out}/qc/species_composition/{sample}.json')
        if options['assay'] == 'microc':
            targets.extend([f'{out}/qc/ligation/{sample}.json', f'{out}/qc/pairtools/{sample}.json'])
        if options['assay'] in ('rnaseq', 'clip'):
            targets.append(f'{out}/qc/rna/{sample}')
        return targets
    if phase in ('report', 'full'):
        return []
    raise ValueError('Unknown portal workflow phase')


def run_phase(root, request, cfg, phase):
    targets = phase_targets(phase, cfg, request['options'])
    if phase in ('prepare', 'full'):
        count=normalize_fastq('/input/input.upload',root/'input.fastq.gz',request['compressed'],request['max_expanded_bytes'],request['max_reads'])
        if request['custom_regions']:
            normalize_bed('/input/regions.upload',root/'regions.bed',cfg['references'][cfg['default_reference']]['chrom_sizes'])
        shutil.copy2('/app/Snakefile',root/'Snakefile')
        shutil.copytree('/app/workflow',root/'workflow',ignore=shutil.ignore_patterns('__pycache__'))
    else:
        previous = json.loads((root/'result.json').read_text())
        if not previous['ok'] or previous['phase'] != {'align':'prepare', 'report':'align'}[phase]:
            raise ValueError('Missing previous phase result')
        count = previous['raw_reads']
    # Partial invocations must keep all temp outputs for the final workflow.
    args = ['snakemake', *targets, '--cores', str(request['cores']), '--resources', 'portal_slot=1',
            '--default-resources', 'portal_slot=1', '--rerun-incomplete', '--printshellcmds']
    if targets:
        args.append('--notemp')
    subprocess.run(args, cwd=root, check=True)
    result = {'ok': True, 'phase': phase, 'raw_reads': count}
    if phase in ('report', 'full'):
        sys.path.insert(0,str(root/'workflow/scripts'))
        from report_bundle import validate
        report=root/'report'
        (report/'.snakemake_timestamp').unlink(missing_ok=True)
        # Public reports contain only container-local paths and the submitted analysis settings.
        (report/'data/submission.json').write_text(json.dumps({'options':request['options'],'raw_reads':count},indent=2)+'\n')
        from report_entry import direct_report
        entry=direct_report(report,request['options']['sample'])
        validation=validate(report)
        manifest=json.loads((report/'manifest.json').read_text())
        manifest.update(entrypoint=entry,validation=validation)
        (report/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
        result['validation'] = validation
    return result


def main():
    root=Path('/job'); os.chdir(root)
    result={'ok':False,'error':'Analysis did not finish.'}
    try:
        request=json.loads((root/'request.json').read_text())
        cfg=yaml.safe_load((root/'config.yaml').read_text())
        result=run_phase(root, request, cfg, os.environ.get('ONT_LIBQC_PHASE', 'full'))
    except InputError as exc: result['error']=str(exc)
    except Exception as exc:
        print(type(exc).__name__+': '+str(exc),file=sys.stderr,flush=True)
        result['error']='The analysis could not complete. Check the assay and reference settings or contact the service operator.'
    (root/'result.json').write_text(json.dumps(result)+'\n')
    return 0 if result['ok'] else 1

if __name__=='__main__':sys.exit(main())
