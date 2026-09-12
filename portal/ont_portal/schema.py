"""The public parameter allowlist. Paths, commands and resource limits are operator-only."""
import math, re
from typing import Literal
from pydantic import BaseModel, ConfigDict, Field, field_validator
from email_validator import validate_email, EmailNotValidError

# label, default, minimum, maximum, numeric kind, UI group
PARAMS = {
 'mapq': ('Minimum mapping quality',20,0,60,'int','Alignment'),
 'minimum_raw_reads': ('Minimum raw reads',500,1,1000000,'int','QC evidence'),
 'minimum_complete_units': ('Minimum complete units',200,1,1000000,'int','QC evidence'),
 'microc_leg_length': ('Micro-C leg length (bp)',100,30,1000,'int','Micro-C'),
 'microc_minimum_contacts': ('Minimum contacts',200,1,1000000,'int','Micro-C'),
 'microc_short_inward_bp': ('Short inward contact limit (bp)',200,0,10000,'int','Micro-C'),
 'microc_minimum_query_coverage': ('Minimum leg coverage',0.8,0.1,1,'float','Micro-C'),
 'microc_long_range_cis_min_bp': ('Long-range cis distance (bp)',10000,1,10000000,'int','Micro-C'),
 'ligation_minimum_segment_bp': ('Minimum ligation segment (bp)',40,20,1000,'int','Micro-C'),
 'ligation_minimum_query_coverage': ('Ligation query coverage',0.8,0.1,1,'float','Micro-C'),
 'ligation_maximum_query_gap': ('Maximum ligation gap (bp)',30,0,1000,'int','Micro-C'),
 'adapter_dimer_max_bp': ('Adapter dimer insert limit (bp)',10,0,100,'int','Library structure'),
 'ultrashort_insert_max_bp': ('Ultrashort insert limit (bp)',50,1,500,'int','Library structure'),
 'adapter_dimer_raw_min_bp': ('Adapter dimer raw minimum (bp)',120,1,1000,'int','Library structure'),
 'adapter_dimer_raw_max_bp': ('Adapter dimer raw maximum (bp)',170,1,2000,'int','Library structure'),
 'adapter_dimer_template_identity': ('Adapter dimer template identity',0.8,0.1,1,'float','Library structure'),
 'single_index_pair_fraction': ('Single index pair fraction',0.9,0.5,1,'float','Index discovery'),
 'minimum_index_informative_units': ('Minimum informative index units',20,1,10000,'int','Index discovery'),
 'minimum_index_cluster_units': ('Minimum index cluster units',20,1,10000,'int','Index discovery'),
 'minimum_index_cluster_fraction': ('Minimum index cluster fraction',0.01,0.001,1,'float','Index discovery'),
 'index_variant_max_edit': ('Index variant maximum edit distance',2,0,3,'int','Index discovery'),
 'index_variant_dominance_ratio': ('Index variant dominance ratio',5,1,100,'float','Index discovery'),
 'annotation_promoter_window_bp': ('Promoter window (bp)',2000,100,10000,'int','Chromatin'),
 'gene_profile_upstream_bp': ('Gene profile upstream (bp)',5000,100,20000,'int','Chromatin'),
 'gene_profile_downstream_bp': ('Gene profile downstream (bp)',5000,100,20000,'int','Chromatin'),
 'gene_profile_body_bp': ('Scaled gene body (bp)',5000,100,20000,'int','Chromatin'),
 'gene_profile_bin_size': ('Gene profile bin size (bp)',100,25,1000,'int','Chromatin'),
 'karyotype_bins_per_longest_chromosome': ('Karyotype bins',250,50,1000,'int','Chromatin'),
 'karyotype_gene_window_bp': ('Karyotype gene window (bp)',250000,1000,1000000,'int','Chromatin'),
 'collision_qc_bin_size': ('Collision O/E bin size (bp)',10000,1000,1000000,'int','Chromatin'),
 'size_trace_max_bp': ('Size plot maximum (bp)',2000,200,20000,'int','Plots'),
 'size_trace_bandwidth': ('Size plot bandwidth',0.045,0.01,0.5,'float','Plots'),
}
RANGES = {
 'raw_reads': {'min':500,'warning_min':250},
 'mean_read_q': None,
 'complete_read_rate': {'min':0.5,'max':1.,'warning_min':0.3,'type':'percent'},
 'adapter_dimer_rate': {'min':0.,'max':0.05,'warning_max':0.1,'type':'percent'},
 'primary_mapping_rate': {'min':0.5,'max':1.,'warning_min':0.3,'type':'percent'},
 'usable_fragment_rate': {'min':0.3,'max':1.,'warning_min':0.15,'type':'percent'},
 'mitochondrial_fraction': {'min':0.,'max':0.2,'warning_max':0.4,'type':'percent'},
 'both_legs_mapping_rate': {'min':0.5,'max':1.,'warning_min':0.3,'type':'percent'},
 'cis_ratio':None,'valid_contact_rate':None,'long_range_cis_fraction':None,
}
SAFE_NAME = r'[A-Za-z][A-Za-z0-9_-]{0,47}'

class StrictModel(BaseModel):
    model_config = ConfigDict(extra='forbid')

class IndexPair(StrictModel):
    name: str = Field(pattern='^'+SAFE_NAME+'$', max_length=48)
    i5: str = Field(pattern=r'^[ACGT]{4,24}$')
    i7: str = Field(pattern=r'^[ACGT]{4,24}$')

class Options(StrictModel):
    sample: str = Field(pattern='^'+SAFE_NAME+'$', max_length=48)
    assay: Literal['chromatin','microc']
    reference: str = Field(pattern=r'^[a-zA-Z0-9_-]{1,32}$')
    mapper: Literal['auto','bowtie2','bwa-mem2','minimap2'] = 'auto'
    spike_in: bool = False
    declared_spike_in: bool = False
    mycoplasma: bool = False
    custom_regions: bool = False
    regions_label: str = Field(default='Uploaded regions',pattern=r'^[A-Za-z0-9 _().-]{1,80}$')
    regions_mode: Literal['auto','reference_point','region'] = 'auto'
    regions_assembly: str | None = Field(default=None,pattern=r'^[A-Za-z0-9_-]{1,32}$')
    index_pairs: list[IndexPair] = Field(default_factory=list,max_length=96)
    params: dict[str,float|int] = Field(default_factory=dict)
    qc_reference_ranges: dict[str,dict[str,float]|None] = Field(default_factory=dict)
    primary_chromosomes: list[str] | None = Field(default=None,max_length=64)
    screen_minimum_reads: int = Field(default=500,ge=1,le=1000000)
    spike_minimum_aligned_bp: int = Field(default=80,ge=20,le=1000)
    spike_minimum_query_coverage: float = Field(default=0.8,ge=0.1,le=1)
    spike_minimum_identity: float = Field(default=0.85,ge=0.1,le=1)

    @field_validator('params')
    @classmethod
    def parameter_allowlist(cls, values):
        for k,v in values.items():
            if k not in PARAMS: raise ValueError('Unsupported parameter: '+k)
            _,_,lo,hi,kind,_ = PARAMS[k]
            if isinstance(v,bool) or not math.isfinite(v) or not lo <= v <= hi:
                raise ValueError(f'{k} must be between {lo} and {hi}')
            if kind == 'int' and v != int(v): raise ValueError(k+' must be an integer')
            values[k] = int(v) if kind == 'int' else float(v)
        full = {k:v[1] for k,v in PARAMS.items()} | values
        if full['adapter_dimer_raw_min_bp'] > full['adapter_dimer_raw_max_bp']:
            raise ValueError('Adapter dimer minimum exceeds maximum')
        if full['adapter_dimer_max_bp'] > full['ultrashort_insert_max_bp']:
            raise ValueError('Adapter dimer limit exceeds ultrashort insert limit')
        size = full['gene_profile_bin_size']
        for k in ['gene_profile_upstream_bp','gene_profile_downstream_bp','gene_profile_body_bp']:
            if full[k] % size: raise ValueError(k+' must be a multiple of gene_profile_bin_size')
        return values

    @field_validator('qc_reference_ranges')
    @classmethod
    def range_allowlist(cls, values):
        for k,limits in values.items():
            if k not in RANGES: raise ValueError('Unsupported QC range: '+k)
            if limits is None: continue
            maximum = 10000000 if k == 'raw_reads' else 60 if k == 'mean_read_q' else 1
            for name,number in limits.items():
                if name not in {'min','max','warning_min','warning_max'} or not math.isfinite(number) or not 0<=number<=maximum:
                    raise ValueError('Invalid cutoff for '+k)
            if limits.get('min',0)>limits.get('max',maximum): raise ValueError('Minimum cutoff exceeds maximum')
            if limits.get('warning_min',0)>limits.get('min',maximum): raise ValueError('Warning minimum exceeds pass minimum')
            if limits.get('warning_max',maximum)<limits.get('max',0): raise ValueError('Warning maximum is below pass maximum')
        return values

    def check_catalog(self, catalog, *, require_regions_confirmation=False):
        refs = {r['id']:r for r in catalog['references']}
        ref = refs.get(self.reference)
        if not ref or not ref['ready']: raise ValueError('Reference is unavailable')
        mapper = self.mapper if self.mapper != 'auto' else ('bwa-mem2' if self.assay=='microc' else 'bowtie2')
        if mapper not in ref['mappers'] or (self.assay=='microc' and mapper=='bowtie2'):
            raise ValueError('Mapper is unavailable for this assay/reference')
        if self.mycoplasma and not ref['mycoplasma']: raise ValueError('Mycoplasma reference is unavailable')
        if self.spike_in and not ref['spike_in']: raise ValueError('Spike-in reference is unavailable')
        if self.custom_regions and self.assay!='chromatin': raise ValueError('Custom regions require a chromatin library')
        # Require confirmation for new submissions while retaining queued legacy jobs.
        if self.custom_regions and (require_regions_confirmation or self.regions_assembly is not None):
            assembly=ref.get('assembly',self.reference)
            if self.regions_assembly!=assembly: raise ValueError(f'Confirm that the BED uses {assembly} coordinates')
        if self.primary_chromosomes is not None:
            if not self.primary_chromosomes or not set(self.primary_chromosomes)<=set(ref['chromosomes']):
                raise ValueError('Choose chromosomes from the selected reference')
        names=[p.name for p in self.index_pairs]; pairs=[(p.i5,p.i7) for p in self.index_pairs]
        if len(set(names))!=len(names) or len(set(pairs))!=len(pairs): raise ValueError('Index pair names and pairs must be unique')
        return self

class Submission(StrictModel):
    email: str = Field(max_length=254)
    filename: str = Field(min_length=1,max_length=180)
    size: int = Field(gt=0)
    options: Options
    consent: Literal[True]
    submission_key: str = Field(default='',max_length=128)
    website: str = Field(default='',max_length=100)

    @field_validator('email')
    @classmethod
    def email_address(cls,v):
        try: return validate_email(v,check_deliverability=False,allow_smtputf8=False).normalized
        except EmailNotValidError as e: raise ValueError('Enter a valid email address') from e

    @field_validator('filename')
    @classmethod
    def filename_only(cls,v):
        if not v.lower().endswith(('.fastq','.fq','.fastq.gz','.fq.gz')):
            raise ValueError('Choose a .fastq, .fq, .fastq.gz or .fq.gz file. POD5 and FAST5 must be basecalled first; changing the extension does not convert the file to FASTQ')
        if not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_. -]*\.(?:fastq|fq)(?:\.gz)?',v,re.I|re.ASCII):
            raise ValueError('Rename the file before uploading: start with an English letter or number and use only English letters, numbers, spaces, dots, hyphens (-) or underscores (_). Parentheses such as (1) are not allowed; rename sample(1).fastq to sample_1.fastq')
        return v

class BatchSubmission(StrictModel):
    # The configured limit can be adjusted up to this defensive parsing ceiling.
    files: list[Submission] = Field(min_length=1, max_length=64)

def schema_payload():
    return {'params':[{ 'key':k,'label':v[0],'default':v[1],'min':v[2],'max':v[3],'kind':v[4],'group':v[5]} for k,v in PARAMS.items()],
            'ranges':RANGES,'metric_profile':'ont_native'}
