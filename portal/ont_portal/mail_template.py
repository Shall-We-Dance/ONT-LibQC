"""Self-contained branded mail: plain text, inline-styled HTML and an embedded logo."""
from datetime import datetime, timezone
from email.message import EmailMessage
from email.utils import formatdate
from html import escape
import json
from pathlib import Path
from urllib.parse import urlsplit
from .store import token_for


def build_message(config,smtp,row,secret,rows=None):
    origin=config['public_url'].rstrip('/');jid=row['id']
    items=rows if rows is not None else [row]
    multiple=len(items)>1
    ready=all(item['status']=='succeeded' for item in items)
    link=origin+'/jobs/'+jid+'#'+token_for(secret,jid)
    title='Your report is ready.' if ready else 'Your analysis needs attention.'
    intro='Your library analysis is complete. View your report and download the results.' if ready else 'Open your task for details and next steps.'
    action='View your report' if ready else 'View your task'
    if rows is not None:
        jid='batch-'+row['batch_id']
    if multiple:
        link=origin+'/#jobs='+','.join(r['id']+':'+token_for(secret,r['id']) for r in items)
        title='Your reports are ready.' if ready else 'Your submission has finished.'
        intro='All results are collected below. Open them together to review and compare your libraries.'
        action='View all results'
    labels={'succeeded':'Report ready','failed':'Failed','cancelled':'Cancelled','expired':'Expired'}
    colors={'succeeded':('#eaf5ef','#27664d'),'failed':('#fff0ec','#963d2e'),'cancelled':('#f0f3f5','#5d717a'),'expired':('#f0f3f5','#5d717a')}
    counts={status:sum(item['status']==status for item in items) for status in labels}
    overview=' · '.join(f'{count} '+{'succeeded':'ready','failed':'failed','cancelled':'cancelled','expired':'expired'}[status] for status,count in counts.items() if count)
    summary=f'{len(items)} '+('library' if len(items)==1 else 'libraries')+' · '+overview
    details=[];cards=[];single_expiry=''
    privacy='These links are private. Share them only with people you trust.' if multiple else 'This link is private. Share it only with people you trust.'
    for index,item in enumerate(items,1):
        opts=json.loads(item['options'] or '{}')
        label=str(opts.get('sample') or item['filename'] or 'Expired task')
        filename=str(item['filename'] or '')
        assay={'microc':'Micro-C','chromatin':'Chromatin'}.get(opts.get('assay'),'')
        meta=' · '.join(str(v) for v in (assay,opts.get('reference')) if v)
        item_link=origin+'/jobs/'+item['id']+'#'+token_for(secret,item['id'])
        until=datetime.fromtimestamp(item['expires'],timezone.utc).strftime('%b %d, %Y · %H:%M UTC')
        status=labels.get(item['status'],item['status'])
        bg,fg=colors.get(item['status'],('#f0f3f5','#5d717a'))
        item_action='View report & downloads' if item['status']=='succeeded' else 'View task details'
        expiry=('Expired on ' if item['status']=='expired' else 'Available until ')+until
        details.append(f'{label} - {status}\n{filename}\n{meta}\n{expiry}\n{item_link}')
        single_expiry=expiry
        badge=f'<span style="display:inline-block;padding:2px 9px;border-radius:4px;background-color:{bg};color:{fg};font-size:12px;font-weight:600;line-height:18px;">{escape(status)}</span>'
        name_html=f'<p style="margin:0 0 6px;font-size:17px;font-weight:700;line-height:25px;color:#152f39;overflow-wrap:anywhere;word-break:break-word;">{escape(label)}</p>'
        file_html=f'<p style="margin:0 0 14px;font-size:12px;line-height:19px;color:#5d717a;overflow-wrap:anywhere;word-break:break-all;">{escape(filename)}</p>' if filename else ''
        if multiple:
            cards.append(f'''<tr><td style="padding:20px 0;border-top:1px solid #dce5e9;">
<p style="margin:0 0 8px;font-size:11px;line-height:17px;color:#5d717a;letter-spacing:1px;">LIBRARY {index:02d}</p>
{name_html}{file_html}
<p style="margin:0 0 12px;font-size:13px;line-height:24px;color:#5d717a;">{badge} &nbsp; {escape(meta)}</p>
<p style="margin:0 0 8px;font-size:14px;font-weight:600;line-height:22px;"><a href="{escape(item_link,quote=True)}" style="color:#117478;text-decoration:underline;">{escape(item_action)} &nbsp;→</a></p>
<p style="margin:0;font-size:12px;line-height:19px;color:#5d717a;">{escape(expiry)}</p>
</td></tr>''')
        else:
            cards.append(f'''<tr><td style="padding:18px;background-color:#f6f9fa;border:1px solid #dce5e9;border-radius:8px;">
<p style="margin:0 0 10px;font-size:11px;line-height:17px;color:#5d717a;letter-spacing:1px;">YOUR LIBRARY</p>
{name_html}{file_html}
<p style="margin:0 0 16px;">{badge}</p>
<table role="presentation" width="100%" cellspacing="0" cellpadding="0" border="0" style="table-layout:fixed;border-top:1px solid #dce5e9;"><tr>
<td width="50%" valign="top" style="padding:12px 10px 0 0;"><p style="margin:0 0 3px;font-size:11px;line-height:17px;color:#5d717a;">Library assay</p><p style="margin:0;font-size:14px;line-height:21px;color:#152f39;overflow-wrap:anywhere;word-break:break-word;">{escape(assay or 'Not available')}</p></td>
<td width="50%" valign="top" style="padding:12px 0 0 10px;"><p style="margin:0 0 3px;font-size:11px;line-height:17px;color:#5d717a;">Reference genome</p><p style="margin:0;font-size:14px;line-height:21px;color:#152f39;overflow-wrap:anywhere;word-break:break-word;">{escape(str(opts.get('reference') or 'Not available'))}</p></td>
</tr></table>
</td></tr>''')
    text=f'{title}\n\n{summary}\n{intro}\n\n{action}:\n{link}\n\n'+'\n\n'.join(details)+'\n\n'+privacy+'\n\nONT-LibQC\nJongmin Kim Lab\n'
    href=escape(link,quote=True)
    action_html=f'''<table role="presentation" cellspacing="0" cellpadding="0" border="0" style="margin:20px 0;"><tr><td align="center" bgcolor="#117478" style="border-radius:6px;mso-padding-alt:12px 20px;"><a href="{href}" style="display:inline-block;padding:12px 20px;border:1px solid #117478;border-radius:6px;color:#ffffff;background-color:#117478;text-decoration:none;font-size:14px;font-weight:700;line-height:20px;">{escape(action)} &nbsp;→</a></td></tr></table>'''
    libraries_html='<table role="presentation" width="100%" cellspacing="0" cellpadding="0" border="0" style="table-layout:fixed;">'+''.join(cards)+'</table>'
    content_html=action_html+libraries_html if multiple else libraries_html+action_html+f'<p style="margin:0 0 16px;font-size:12px;line-height:19px;color:#5d717a;">{escape(single_expiry)}</p>'
    eyebrow=summary if multiple else ('ANALYSIS COMPLETE' if ready else 'ANALYSIS UPDATE')
    html=f'''<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><meta name="color-scheme" content="light"><title>{escape(title)}</title></head>
<body style="margin:0;padding:0;background-color:#f2f6f8;color:#152f39;font-family:Arial,Helvetica,sans-serif;-webkit-text-size-adjust:100%;">
<div style="display:none;font-size:1px;line-height:1px;max-height:0;max-width:0;opacity:0;overflow:hidden;mso-hide:all;">{escape(summary)}. {escape(title)}</div>
<table role="presentation" width="100%" cellspacing="0" cellpadding="0" border="0" style="background-color:#f2f6f8;"><tr><td align="center" style="padding:24px 12px;">
<!--[if mso]><table role="presentation" width="600" align="center"><tr><td><![endif]-->
<table role="presentation" width="100%" cellspacing="0" cellpadding="0" border="0" style="max-width:600px;table-layout:fixed;background-color:#ffffff;border:1px solid #dce5e9;border-radius:12px;overflow:hidden;">
<tr><td style="background-color:#123a45;padding:22px 24px;">
<table role="presentation" cellspacing="0" cellpadding="0" border="0"><tr><td width="36" style="vertical-align:middle;"><img src="cid:ontlibqc-logo" width="32" height="32" alt="" style="display:block;width:32px;height:32px;border:0;"></td><td style="vertical-align:middle;padding-left:8px;color:#f6fffd;font-size:21px;font-weight:700;line-height:28px;">ONT-LibQC</td></tr></table>
</td></tr>
<tr><td style="padding:28px 24px 8px;">
<p style="margin:0 0 10px;color:#117478;font-size:11px;font-weight:600;line-height:18px;letter-spacing:0.5px;">{escape(eyebrow)}</p>
<h1 style="margin:0 0 12px;color:#152f39;font-size:26px;font-weight:700;line-height:33px;letter-spacing:-0.5px;">{escape(title)}</h1>
<p style="margin:0 0 20px;color:#5d717a;font-size:14px;line-height:23px;">{escape(intro)}</p>
{content_html}
</td></tr>
<tr><td style="padding:18px 24px 22px;border-top:1px solid #dce5e9;background-color:#f8fafb;">
<p style="margin:0 0 8px;color:#5d717a;font-size:12px;line-height:19px;">{escape(privacy)}</p>
<p style="margin:0;color:#5d717a;font-size:12px;line-height:19px;">Button not opening? <a href="{href}" style="color:#117478;text-decoration:underline;">{escape(action)}</a></p>
</td></tr></table>
<!--[if mso]></td></tr></table><![endif]-->
<p style="margin:16px 0 0;color:#5d717a;font-size:12px;line-height:19px;">Jongmin Kim Lab</p>
</td></tr></table></body></html>'''
    message=EmailMessage()
    message['Subject']='Your ONT-LibQC report is ready' if ready else 'Your ONT-LibQC analysis needs attention'
    if multiple:message.replace_header('Subject',f'Your ONT-LibQC results: {len(items)} files - '+('reports ready' if ready else 'submission finished'))
    message['From']=smtp['from_address'];message['To']=row['email']
    message['Date']=formatdate(row['created'],usegmt=True)
    message['Message-ID']='<ont-libqc-'+jid+'-'+('complete' if rows is not None else row['status'])+'@'+urlsplit(origin).hostname+'>'
    if smtp.get('provider')=='resend':message['Resend-Idempotency-Key']='ont-libqc/'+jid+'/'+('complete' if rows is not None else row['status'])
    message.set_content(text);message.add_alternative(html,subtype='html')
    part=message.get_payload()[-1]
    logo=(Path(__file__).resolve().parents[1]/'static/mail-logo.png').read_bytes()
    part.add_related(logo,maintype='image',subtype='png',cid='<ontlibqc-logo>',filename='ontlibqc-logo.png',disposition='inline')
    # Stable boundaries keep repeated SMTP attempts identical for Resend.
    message.set_boundary('ontlibqc_alternative_'+jid)
    part.set_boundary('ontlibqc_related_'+jid)
    return message
