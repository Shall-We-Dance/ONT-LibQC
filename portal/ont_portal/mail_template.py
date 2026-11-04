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
    title={'succeeded':'Report ready','failed':'Analysis failed','cancelled':'Analysis cancelled',
           'expired':'Report unavailable'}.get(row['status'],'Analysis update')
    subject='ONT-LibQC - '+title.lower()
    action='Open report' if ready else 'View analysis'
    if rows is not None:
        jid='batch-'+row['batch_id']
    if multiple:
        link=origin+'/#jobs='+','.join(r['id']+':'+token_for(secret,r['id']) for r in items)
        title='Reports ready' if ready else 'Submission finished'
        subject=f'ONT-LibQC - {len(items)} '+('reports ready' if ready else 'analyses finished')
        action='View all results'
    labels={'succeeded':'Report ready','failed':'Failed','cancelled':'Cancelled','expired':'Expired','deleted':'Deleted'}
    colors={'succeeded':'#256548','failed':'#ac3333','cancelled':'#5b6d82','expired':'#5b6d82','deleted':'#5b6d82'}
    def display_status(item):
        if item['status']=='expired' and (item.get('deletion_reason')=='manual' or item['expires']==0):
            return 'deleted'
        return item['status']
    counts={status:sum(display_status(item)==status for item in items) for status in labels}
    overview=' · '.join(f'{count} '+{'succeeded':'ready','failed':'failed','cancelled':'cancelled','expired':'expired','deleted':'deleted'}[status] for status,count in counts.items() if count)
    summary=f'{len(items)} '+('library' if len(items)==1 else 'libraries')+' · '+overview
    details=[];cards=[];single_expiry=''
    privacy='These links are private. Share them only with people you trust.' if multiple else 'This link is private. Share it only with people you trust.'
    mono="'SFMono-Regular',Consolas,'Liberation Mono',monospace"
    for item in items:
        opts=json.loads(item['options'] or '{}')
        state=display_status(item)
        label=str(opts.get('sample') or item['filename'] or ('Deleted report' if state=='deleted' else 'Expired report'))
        filename=str(item['filename'] or '')
        assay={'microc':'Micro-C','chromatin':'Chromatin'}.get(opts.get('assay'),'')
        meta=' · '.join(str(v) for v in (assay,opts.get('reference')) if v)
        item_link=origin+'/jobs/'+item['id']+'#'+token_for(secret,item['id'])
        status=labels.get(state,state)
        fg=colors.get(state,'#5b6d82')
        if state=='deleted':
            expiry='Manually deleted. This report is no longer available.'
        else:
            until=datetime.fromtimestamp(item['expires'],timezone.utc).strftime('%b %d, %Y · %H:%M UTC')
            expiry=('Expired on ' if state=='expired' else 'Available until ' if state=='succeeded' else 'Task details available until ')+until
        details.append(f'{label} - {status}\n{filename}\n{meta}\n{expiry}\n{item_link}')
        single_expiry=expiry
        name_html=f'<p style="margin:0 0 5px;font-family:{mono};font-size:15px;font-weight:500;line-height:23px;color:#24364b;overflow-wrap:anywhere;word-break:break-word;">{escape(label)}</p>'
        file_html=f'<p style="margin:0 0 10px;font-family:{mono};font-size:12px;line-height:19px;color:#5b6d82;overflow-wrap:anywhere;word-break:break-all;">{escape(filename)}</p>' if filename else ''
        status_html=f'<span style="color:{fg};font-size:12px;font-weight:500;">{escape(status)}</span>'
        meta_html=f' &nbsp;·&nbsp; {escape(meta)}' if meta else ''
        item_action='Open report' if state=='succeeded' else 'View analysis'
        item_details=f'''{name_html}{file_html}
<p style="margin:0;font-size:13px;line-height:22px;color:#5b6d82;">{status_html}{meta_html}</p>'''
        if multiple:
            cards.append(f'''<tr><td style="padding:18px 0;border-top:1px solid #dce5f0;">
{item_details}
<p style="margin:10px 0 6px;font-size:13px;line-height:21px;"><a href="{escape(item_link,quote=True)}" style="color:#2a6cbd;text-decoration:underline;">{escape(item_action)}</a></p>
<p style="margin:0;font-size:12px;line-height:19px;color:#5b6d82;">{escape(expiry)}</p>
</td></tr>''')
        else:
            cards.append(f'''<tr><td style="padding:18px;border:1px solid #dce5f0;border-radius:8px;background-color:#f5f8fc;">{item_details}</td></tr>''')
    text=f'{title}\n\n'+(summary+'\n\n' if multiple else '')+f'{action}:\n{link}\n\n'+'\n\n'.join(details)+'\n\n'+privacy+'\n\nONT-LibQC\nJongmin Kim Lab\n'
    href=escape(link,quote=True)
    action_html=f'''<table role="presentation" cellspacing="0" cellpadding="0" border="0" style="margin:18px 0;"><tr><td align="center" bgcolor="#2a6cbd" style="border-radius:6px;mso-padding-alt:11px 18px;"><a href="{href}" style="display:inline-block;padding:11px 18px;border:1px solid #2a6cbd;border-radius:6px;color:#ffffff;background-color:#2a6cbd;text-decoration:none;font-size:14px;font-weight:600;line-height:20px;">{escape(action)}</a></td></tr></table>'''
    libraries_html='<table role="presentation" width="100%" cellspacing="0" cellpadding="0" border="0" style="table-layout:fixed;">'+''.join(cards)+'</table>'
    content_html=action_html+libraries_html if multiple else libraries_html+action_html+f'<p style="margin:0 0 18px;font-size:12px;line-height:19px;color:#5b6d82;">{escape(single_expiry)}</p>'
    overview_html=f'<p style="margin:6px 0 0;color:#5b6d82;font-size:13px;line-height:21px;">{escape(summary)}</p>' if multiple else ''
    html=f'''<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><meta name="color-scheme" content="light"><title>{escape(title)}</title></head>
<body style="margin:0;padding:0;background-color:#f5f8fc;color:#24364b;font-family:Arial,Helvetica,sans-serif;-webkit-text-size-adjust:100%;">
<div style="display:none;font-size:1px;line-height:1px;max-height:0;max-width:0;opacity:0;overflow:hidden;mso-hide:all;">{escape(summary)}. {escape(title)}</div>
<table role="presentation" width="100%" cellspacing="0" cellpadding="0" border="0" style="background-color:#f5f8fc;"><tr><td align="center" style="padding:20px 12px;">
<!--[if mso]><table role="presentation" width="600" align="center"><tr><td><![endif]-->
<table role="presentation" width="100%" cellspacing="0" cellpadding="0" border="0" style="max-width:600px;table-layout:fixed;background-color:#ffffff;border:1px solid #dce5f0;border-radius:8px;overflow:hidden;">
<tr><td style="padding:20px 22px;border-bottom:1px solid #dce5f0;">
<table role="presentation" cellspacing="0" cellpadding="0" border="0"><tr><td width="32" style="vertical-align:middle;"><img src="cid:ontlibqc-logo" width="32" height="32" alt="" style="display:block;width:32px;height:32px;border:0;"></td><td style="vertical-align:middle;padding-left:10px;color:#24364b;font-size:20px;font-weight:600;line-height:28px;">ONT-LibQC</td></tr></table>
</td></tr>
<tr><td style="padding:22px 22px 8px;">
<h1 style="margin:0 0 18px;color:#24364b;font-size:24px;font-weight:600;line-height:31px;letter-spacing:-0.3px;">{escape(title)}</h1>
{overview_html}{content_html}
</td></tr>
<tr><td style="padding:16px 22px 20px;border-top:1px solid #dce5f0;">
<p style="margin:0 0 8px;color:#5b6d82;font-size:12px;line-height:19px;">{escape(privacy)}</p>
<p style="margin:0;color:#5b6d82;font-size:12px;line-height:19px;">Jongmin Kim Lab</p>
</td></tr></table>
<!--[if mso]></td></tr></table><![endif]-->
</td></tr></table></body></html>'''
    message=EmailMessage()
    message['Subject']=subject
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
