"""Operator-owned reference catalog and read-only container mounts."""
import copy, json, re, struct
from pathlib import Path
from .mappers import resolve_mapper, microc_bwa_mapper, microc_default

BT2_PARTS=('1','2','3','4','rev.1','rev.2')
def bt2_ready(prefix):
    return bool(prefix and any(all(Path(str(prefix)+'.'+part+'.'+suffix).is_file() for part in BT2_PARTS) for suffix in ('bt2','bt2l')))
def bwa_ready(prefix):
    return bool(prefix and all(Path(str(prefix)+'.'+part).is_file() for part in ('0123','amb','ann','bwt.2bit.64','pac')))
def bwa_mem_ready(prefix):
    if not prefix or not all(Path(str(prefix)+'.'+part).is_file() for part in ('amb','ann','bwt','pac','sa')):
        return False
    try:
        length=int(Path(str(prefix)+'.ann').read_text().split()[0])
        with Path(str(prefix)+'.bwt').open('rb') as source: header=source.read(40)
        with Path(str(prefix)+'.sa').open('rb') as source: sa=source.read(56)
        return (length>0 and len(header)==40 and len(sa)==56 and header==sa[:40]
                and struct.unpack('<5Q',header)[4]==2*length
                and struct.unpack('<7Q',sa)[6]==2*length
                and struct.unpack('<7Q',sa)[5]>0
                and Path(str(prefix)+'.sa').stat().st_size==56+(2*length//struct.unpack('<7Q',sa)[5])*8
                and Path(str(prefix)+'.pac').stat().st_size==length//4+2)
    except (OSError,ValueError,IndexError,struct.error):return False
def minimap_ready(value):
    if not value:return False
    try:
        with Path(value).open('rb') as source:header=source.read(24)
        if len(header)!=24 or header[:4]!=b'MMI\x02':return False
        window,kmer,buckets,sequences,flags=struct.unpack('<5I',header[4:])
        return window==5 and kmer==9 and 0<buckets<=30 and sequences>0 and flags==0 and Path(value).stat().st_size>24
    except OSError:return False

def public_catalog(catalog):
    result=[]
    for key,r in catalog['references'].items():
        sizes=Path(r['chrom_sizes'])
        chromosomes=[line.split()[0] for line in sizes.read_text().splitlines() if line.strip()] if sizes.is_file() else []
        mappers=['minimap2'] if minimap_ready(r.get('minimap2_index')) else []
        if bt2_ready(r.get('bowtie2_index')): mappers.append('bowtie2')
        if bwa_ready(r.get('bwa_mem2_index')): mappers.append('bwa-mem2')
        if bwa_mem_ready(r.get('bwa_index')): mappers.append('bwa-mem')
        mappers=[m for m in mappers if m not in r.get('disabled_mappers',[])]
        spike=r.get('spike_in',{})
        if spike.get('mode')=='auto':
            manifest=Path(spike['database_manifest'])
            db=json.loads(manifest.read_text()) if manifest.is_file() else {}
            h=db.get('host_indexes',{}).get(key,{})
            spike_ready=bt2_ready(h.get('index')) and Path(h.get('contigs','/missing')).is_file()
        else: spike_ready=bt2_ready(spike.get('combined_bowtie2_index'))
        result.append({'id':key,'label':r.get('label',key),'assembly':r.get('assembly',key),
            'ready':all(Path(r.get(k,'/missing')).is_file() for k in ('fasta','chrom_sizes','annotation')) and 'bowtie2' in mappers and any(m in mappers for m in ('bwa-mem2','bwa-mem')),
            'microc_mapper':microc_default(dict(r,mappers=mappers)),
            'availability_note':r.get('availability_note',''),
            'mappers':mappers,'mycoplasma':bt2_ready(r.get('mycoplasma_index')),'spike_in':spike_ready,
            'chromosomes':chromosomes,'primary_chromosomes':r['primary_chromosomes']})
    return {'references':result}

class ReferenceMounts:
    def __init__(self): self.parents={}
    def path(self,value,prefix=False):
        p=Path(value).absolute()
        # Prefixes need not exist, but their parent must be a real directory.
        parent=p.parent.resolve()
        if not parent.is_dir() or parent==Path('/') or not re.fullmatch(r'[A-Za-z0-9_.-]+',p.name):
            raise ValueError('Unsupported operator reference path')
        # Resolve file symlinks before mounting their directory.
        if p.is_symlink() and not prefix: p=p.resolve();parent=p.parent
        if parent not in self.parents: self.parents[parent]='/refs/r'+str(len(self.parents))
        return self.parents[parent]+'/'+p.name
    def arguments(self):
        args=[]
        for source,dest in self.parents.items():
            if ',' in str(source): raise ValueError('Commas are not supported in reference directory paths')
            args+=['--mount',f'type=bind,src={source},dst={dest},readonly']
        return args

def job_config(options,catalog,jobdir,cores):
    from .schema import PARAMS, RANGES
    data=options.model_dump(); ref=copy.deepcopy(catalog['references'][options.reference]); mounts=ReferenceMounts()
    reference={k:mounts.path(ref[k],prefix=k.endswith('_index')) for k in ('fasta','chrom_sizes','annotation','bowtie2_index','bwa_mem2_index','bwa_index') if ref.get(k)}
    mapper=resolve_mapper(options.mapper,options.assay,ref)
    reference['microc_mapper']=microc_bwa_mapper(mapper,ref)
    reference.update(assembly=ref.get('assembly',options.reference),primary_chromosomes=options.primary_chromosomes or ref['primary_chromosomes'])
    if options.mapper=='minimap2': reference['minimap2_index']=mounts.path(ref['minimap2_index'])
    if options.mycoplasma: reference['mycoplasma_index']=mounts.path(ref['mycoplasma_index'],prefix=True)
    spike={'enabled':options.spike_in,'minimum_screened_reads':options.screen_minimum_reads,
        'minimum_aligned_bp':options.spike_minimum_aligned_bp,'minimum_query_coverage':options.spike_minimum_query_coverage,
        'minimum_identity':options.spike_minimum_identity,'declared_spike_in':options.declared_spike_in}
    if options.spike_in:
        spec=ref['spike_in'];spike.update({k:v for k,v in spec.items() if k not in {'database_manifest','combined_bowtie2_index'}})
        if spec.get('mode')=='auto':
            db=json.loads(Path(spec['database_manifest']).read_text()); h=db['host_indexes'][options.reference]
            safe={'version':db.get('version','operator-panel'),'entries':[{k:e[k] for k in ('id','label','assembly','role') if k in e} for e in db['entries']],
                'host_indexes':{options.reference:{'index':mounts.path(h['index'],prefix=True),'contigs':mounts.path(h['contigs']),'excluded_candidates':h.get('excluded_candidates',[])}}}
            (jobdir/'spike-panel.json').write_text(json.dumps(safe,indent=2)+'\n');spike['database_manifest']='/job/spike-panel.json'
        else: spike['combined_bowtie2_index']=mounts.path(spec['combined_bowtie2_index'],prefix=True)
    sample={'fastq':['/job/input.fastq.gz'],'type':options.assay,'reference':options.reference,'index_pairs':[p.model_dump() for p in options.index_pairs]}
    sample['library_structure']=options.library_structure.model_dump()
    sample['mapper']=mapper
    if options.custom_regions:sample['custom_regions']={'bed':'/job/regions.bed','assembly':reference['assembly'],'label':options.regions_label,'mode':options.regions_mode,**({'bed_name':options.regions_filename} if options.regions_filename else {})}
    ranges=copy.deepcopy(RANGES)
    for k,v in options.qc_reference_ranges.items():
        ranges[k]=None if v is None else dict(v,**({'type':'percent'} if k not in {'raw_reads','mean_read_q'} else {}))
    params={k:v[1] for k,v in PARAMS.items()} | options.params
    params.update(alignment_threads=cores,gene_profile_threads=cores,gene_profile_sublibrary_threads=cores)
    cfg={'metric_profile':'ont_native','output_dir':'.work','report_dir':'report','default_reference':options.reference,
         'alignment':{'bwa_batch_bases':1000000,'samtools_sort_memory':'128M'},
         'references':{options.reference:reference},'samples':{options.sample:sample},'params':params,'qc_reference_ranges':ranges,
         'spike_in':spike,'mycoplasma_screen':{'enabled':options.mycoplasma,'minimum_screened_reads':options.screen_minimum_reads,
          'target_contig_regex':catalog.get('mycoplasma_regex','mycoplasma_target|mycoplasmopsis|mycoplasmoides|metamycoplasma|mesomycoplasma|acholeplasma|spiroplasma|ureaplasma|mollicutes')}}
    return cfg,mounts
