"""Shared definitions for native ONT QC, molecule downsampling and paired QC."""
import json,math
from pathlib import Path

DEFINITIONS_PATH=Path(__file__).with_name('metric_definitions.json')
DEFINITIONS=json.loads(DEFINITIONS_PATH.read_text())

def profile(name='ont_native'):
    return dict(DEFINITIONS['common'],**DEFINITIONS['profiles'][name],metric_profile=name,metric_definition_version=DEFINITIONS['version'])

def configured_params(config):
    if config.get('metric_profile','ont_native') != 'ont_native':
        raise ValueError('03/05 use ont_native extraction and deduplication. Run 07 for matched_comparison; changing only a profile label cannot change the extraction method.')
    return dict(profile(config.get('metric_profile','ont_native')),**config.get('params',{}))

def wilson(success,total,z=1.96):
    if not total:return None,None
    p=success/total;den=1+z*z/total
    center=(p+z*z/(2*total))/den
    half=z*math.sqrt(p*(1-p)/total+z*z/(4*total*total))/den
    return max(0.,center-half),min(1.,center+half)

def contact_evidence_status(total,minimum=None):
    minimum=DEFINITIONS['common']['microc_minimum_contacts'] if minimum is None else minimum
    return 'insufficient_evidence_no_contacts' if not total else 'insufficient_evidence_low_contacts' if total<minimum else 'sufficient_contact_count'

def summarize_contacts(distances,total,long_range=None,minimum=None):
    long_range=DEFINITIONS['common']['microc_long_range_cis_min_bp'] if long_range is None else long_range
    distances=list(distances);cis=len(distances);long_cis=sum(d>=long_range for d in distances)
    low,high=wilson(cis,total);lr_low,lr_high=wilson(long_cis,total)
    return dict(cis_contacts=cis,trans_contacts=total-cis,cis_ratio=cis/total if total else None,
                cis_ratio_ci95_low=low,cis_ratio_ci95_high=high,long_range_cis_min_bp=long_range,
                long_range_cis_contacts=long_cis,long_range_cis_fraction=long_cis/total if total else None,
                long_range_cis_ci95_low=lr_low,long_range_cis_ci95_high=lr_high,
                long_range_within_cis_fraction=long_cis/cis if cis else None,
                cis_lt1kb=sum(d<1000 for d in distances),cis_1_10kb=sum(1000<=d<10000 for d in distances),
                cis_10_100kb=sum(10000<=d<100000 for d in distances),cis_ge100kb=sum(d>=100000 for d in distances),
                contact_evidence_status=contact_evidence_status(total,minimum),biological_qc_call='not_classified')

def metric_settings(args):
    result=profile(getattr(args,'metric_profile','ont_native'))
    for dest,src in [('mapq','mapq'),('collision_qc_bin_size','collision_bin_size'),('microc_minimum_contacts','microc_minimum_contacts'),('microc_short_inward_bp','microc_short_inward_bp'),('microc_minimum_query_coverage','microc_minimum_query_coverage'),('microc_long_range_cis_min_bp','microc_long_range_cis_min_bp'),('minimum_raw_reads','minimum_raw'),('minimum_complete_units','minimum_units')]:
        if hasattr(args,src):result[dest]=getattr(args,src)
    return result
