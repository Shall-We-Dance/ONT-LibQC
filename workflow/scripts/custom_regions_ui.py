"""Shared presentation for newly generated and already published BED reports."""
import html
import json
from pathlib import Path
import re


def render_section(payload, initial='', requested_mode='auto', bed_name=None):
    reference=html.escape(payload['reference'])
    payload=dict(payload, submitted_bed_name=bed_name)
    submitted=''
    if initial:
        source=(': <strong>'+html.escape(bed_name)+'</strong>.') if bed_name else '. The original filename was not retained for this report.'
        submitted=("<div id='custom-submitted'><h3>Submitted BED</h3>"
                   "<p class='muted'>The plot below uses the BED submitted with this analysis"+source+
                   '</p>'+initial+'</div>')
    options=''.join(f"<option value='{html.escape(s,quote=True)}'>{html.escape(x['name'])}</option>" for s,x in payload['scopes'].items())
    modes=''.join(f"<option value='{value}'{' selected' if value==requested_mode else ''}>{label}</option>" for value,label in [('auto','Auto'),('reference_point','Reference point'),('region','Region (scaled)')])
    heading='Analyze another BED' if initial else 'Analyze a BED file'
    explanation=('The submitted-BED plot stays above for comparison. ' if initial else '')+'Choose a BED, confirm its coordinates, and click Draw metaplot to calculate a new profile from this report’s library data. The file stays in your browser.'
    script=Path(__file__).with_name('custom_regions.js').read_text()
    encoded=json.dumps(payload,separators=(',',':')).replace('<',r'\u003c')
    return f"""<section class='report-section' id='custom-regions'><div class='section-heading'><span class='eyebrow'>TARGET CONTEXT</span><h2>Custom region metaplot</h2></div>{submitted}
<h3 style='margin:24px 0 8px'>{heading}</h3><p class='muted'>{explanation}</p>
<div class='custom-controls'><label class='custom-field'><span class='custom-field-title'>Region BED</span><span class='custom-file-picker'><input type='file' id='custom-bed' accept='.bed,.tsv,.txt'><span class='custom-file-button' aria-hidden='true'>Choose BED</span><span class='custom-file-name' id='custom-bed-name' aria-hidden='true'>No new BED selected</span></span></label><label class='custom-field'><span class='custom-field-title'>Library</span><span class='custom-select-shell'><select id='custom-scope'>{options}</select></span></label><label class='custom-field'><span class='custom-field-title'>Plot mode</span><span class='custom-select-shell'><select id='custom-mode'>{modes}</select></span></label><div class='custom-actions'><label class='custom-confirm'><input type='checkbox' id='custom-assembly'><span> I confirm this BED uses {reference} coordinates</span></label><button id='custom-run' type='button'>Draw metaplot</button></div></div>
<p id='custom-message' role='status' aria-live='polite'></p><div id='custom-local' hidden><h3>New BED result</h3><p id='custom-current-source' class='muted'></p><div id='custom-output'></div><div id='custom-downloads'></div></div>
<details><summary>Normalization and region handling</summary><p>BED uses zero-based, half-open coordinates. BED6 minus strands are reversed; BED3/4 and dot strands use genomic orientation. Auto uses reference point mode when all retained intervals are 1 bp; otherwise it uses region mode, including mixed-width BED files. You can override this with Plot mode. Reference point mode centers each interval at its midpoint (start + end)/2 and shows ±5 kb in 100 bins without scaling a body. Region mode scales bodies to 100 bins with 50 bins across each 5 kb flank. Exact aligned-base coverage (CIGAR M/= /X only; D/N excluded) is divided by each bin's available genomic width, then CPM uses all eligible primary nuclear alignments in the selected scope. Regions are equally weighted. Duplicate intervals are removed; non-primary and absent-contig regions are skipped and counted. Flanks are clipped to chromosome boundaries; bins entirely beyond the chromosome are missing. SEM is across regions. Assembly identity must be confirmed; chromosome names alone cannot verify it. No data are uploaded to a server.</p></details>
<script id='custom-payload' type='application/json'>{encoded}</script><script>{script}</script></section>"""


def upgrade_html(text, bed_name=None):
    """Refresh controls only; preserve original plots, records and calculations."""
    pattern=r"<section\b[^>]*\bid=['\"]custom-regions['\"][^>]*>.*?<script\b[^>]*\bid=['\"]custom-payload['\"][^>]*>(.*?)</script>\s*<script>.*?</script></section>"
    def replace(match):
        section=match.group(0);payload=json.loads(match.group(1))
        if "id='custom-submitted'" in section:
            start=section.index("id='custom-submitted'")
            end=section.index('</div>',start)
            old=section[start:end]
            initial=old[old.index('</p>')+4:]
        else:
            heading=re.search(r'<h2>Custom region metaplot</h2></div>',section)
            controls=re.search(r"<p class=['\"]muted['\"]>Choose a BED file",section)
            initial=section[heading.end():controls.start()].strip() if heading and controls else ''
        mode=re.search(r"<option value=['\"]([^'\"]+)['\"] selected",section)
        return render_section(payload,initial,mode.group(1) if mode else 'auto',bed_name or payload.get('submitted_bed_name'))
    return re.sub(pattern,replace,text,flags=re.S)
