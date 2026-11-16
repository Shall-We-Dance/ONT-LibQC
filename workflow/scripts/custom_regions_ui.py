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
    explanation='Choose a BED file, confirm the genome assembly, then click Draw metaplot. Runs in your browser; files can contain up to 100,000 regions and be up to 10 MiB.'
    script=Path(__file__).with_name('custom_regions.js').read_text()
    encoded=json.dumps(payload,separators=(',',':')).replace('<',r'\u003c')
    return f"""<section class='report-section' id='custom-regions'><div class='section-heading'><span class='eyebrow'>TARGET CONTEXT</span><h2>Custom region metaplot</h2></div>{submitted}
<style>#custom-output,.custom-static{{max-width:850px;margin:12px auto 0}}#custom-output svg{{height:auto}}@media print{{#custom-output{{max-width:100%}}}}</style>
<h3 style='margin:24px 0 8px'>{heading}</h3><p class='muted'>{explanation}</p>
<div class='custom-controls'><label class='custom-field'><span class='custom-field-title'>Region BED</span><span class='custom-file-picker'><input type='file' id='custom-bed'><span class='custom-file-button' aria-hidden='true'>Choose BED</span><span class='custom-file-name' id='custom-bed-name' aria-hidden='true'>No new BED selected</span></span></label><label class='custom-field'><span class='custom-field-title'>Library</span><span class='custom-select-shell'><select id='custom-scope'>{options}</select></span></label><label class='custom-field'><span class='custom-field-title'>Plot mode</span><span class='custom-select-shell'><select id='custom-mode'>{modes}</select></span></label><div class='custom-actions'><label class='custom-confirm'><input type='checkbox' id='custom-assembly'><span> I confirm this BED uses {reference} coordinates</span></label><button id='custom-run' type='button'>Draw metaplot</button></div></div>
<p id='custom-message' role='status' aria-live='polite'></p><div id='custom-local' hidden><h3>New BED result</h3><p id='custom-current-source' class='muted'></p><div id='custom-output'></div><div id='custom-downloads'></div></div>
<details><summary>About this plot</summary>
<p>Use standard BED coordinates from the same genome assembly as this report. Chromosome names work with or without “chr”.</p>
<p>Auto centers the plot on sites when every region is 1 bp wide. Otherwise, it scales regions to the same width. Both modes include 5 kb on either side; you can change the mode above. Regions marked with a minus strand are flipped to follow their direction.</p>
<p>The line shows average coverage in counts per million (CPM); the shaded band shows the standard error across regions. Each region has equal weight. Duplicate regions and regions outside the reference’s primary chromosomes are skipped. Windows are trimmed at chromosome ends.</p></details>
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
