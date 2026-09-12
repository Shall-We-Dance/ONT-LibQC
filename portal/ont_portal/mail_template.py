"""Self-contained branded mail: plain text, inline-styled HTML and an embedded logo."""
from datetime import datetime, timezone
from email.message import EmailMessage
from email.utils import formatdate
from html import escape
import json
from pathlib import Path
from urllib.parse import urlsplit
from .store import token_for


def build_message(config,smtp,row,secret):
    origin=config['public_url'].rstrip('/');jid=row['id'];ready=row['status']=='succeeded'
    link=origin+'/jobs/'+jid+'#'+token_for(secret,jid)
    options=json.loads(row['options'] or '{}')
    sample=str(options.get('sample') or 'Your library')
    reference=str(options.get('reference') or '')
    assay='Micro-C' if options.get('assay')=='microc' else 'Chromatin'
    title='Your report is ready.' if ready else 'Your analysis needs attention.'
    intro='Your library analysis is complete. Open your private task page to explore the report or download your results.' if ready else 'We couldn’t complete this analysis. Open your private task page for details and the next steps.'
    action='View your report' if ready else 'View your task'
    expiry=datetime.fromtimestamp(row['expires'],timezone.utc).strftime('%b %d, %Y, %H:%M UTC')
    retention=('Report available until '+expiry+'.') if ready else ('Task available until '+expiry+'.')
    text=f'{title}\n\n{intro}\n\n{sample}\n{assay} · {reference}\n\n{action}:\n{link}\n\n{retention}\nThis private link grants access to your task. Share it only with people you trust.\n\nONT-LibQC\nJongmin Kim Lab\n'
    href=escape(link,quote=True)
    html=f'''<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><meta name="color-scheme" content="light"><title>{escape(title)}</title></head>
<body style="margin:0;padding:0;background-color:#f2f6f8;color:#152f39;font-family:Arial,Helvetica,sans-serif;-webkit-text-size-adjust:100%;">
<div style="display:none;font-size:1px;color:#f2f6f8;line-height:1px;max-height:0;max-width:0;opacity:0;overflow:hidden;mso-hide:all;">{escape(title)} Open your private ONT-LibQC task page.</div>
<table role="presentation" width="100%" cellspacing="0" cellpadding="0" border="0" style="background-color:#f2f6f8;"><tr><td align="center" style="padding:28px 12px;">
<!--[if mso]><table role="presentation" width="560" align="center"><tr><td><![endif]-->
<table role="presentation" width="100%" cellspacing="0" cellpadding="0" border="0" style="max-width:560px;background-color:#ffffff;border:1px solid #dce5e9;border-radius:14px;overflow:hidden;">
<tr><td style="background-color:#123a45;padding:26px 24px;">
<table role="presentation" cellspacing="0" cellpadding="0" border="0"><tr><td width="48" style="vertical-align:middle;"><img src="cid:ontlibqc-logo" width="40" height="40" alt="ONT-LibQC logo" style="display:block;width:40px;height:40px;border:0;"></td><td style="vertical-align:middle;padding-left:8px;color:#f6fffd;font-size:24px;font-weight:700;letter-spacing:-0.5px;">ONT-LibQC</td></tr></table>
<p style="margin:18px 0 0;color:#8ce4ce;font-size:11px;font-weight:700;letter-spacing:1.7px;line-height:1.5;">LIBRARY QUALITY CONTROL</p>
</td></tr>
<tr><td style="padding:30px 24px 26px;">
<h1 style="margin:0 0 14px;color:#152f39;font-size:28px;font-weight:700;line-height:1.2;letter-spacing:-0.6px;">{escape(title)}</h1>
<p style="margin:0 0 24px;color:#5d717a;font-size:16px;line-height:1.65;">{escape(intro)}</p>
<table role="presentation" width="100%" cellspacing="0" cellpadding="0" border="0" style="background-color:#f2f6f8;border:1px solid #dce5e9;border-radius:8px;"><tr><td style="padding:16px 18px;">
<p style="margin:0 0 5px;color:#152f39;font-size:16px;font-weight:700;line-height:1.5;word-break:break-all;">{escape(sample)}</p>
<p style="margin:0;color:#5d717a;font-size:14px;line-height:1.5;">{escape(assay)} &nbsp;·&nbsp; {escape(reference)}</p>
</td></tr></table>
<table role="presentation" cellspacing="0" cellpadding="0" border="0" style="margin:26px 0 22px;"><tr><td align="center" bgcolor="#117478" style="border-radius:7px;mso-padding-alt:14px 22px;"><a href="{href}" style="display:inline-block;padding:14px 22px;border:1px solid #117478;border-radius:7px;color:#ffffff;background-color:#117478;text-decoration:none;font-size:16px;font-weight:700;line-height:1.3;">{action} &nbsp;→</a></td></tr></table>
<p style="margin:0;color:#5d717a;font-size:13px;line-height:1.6;">{escape(retention)}</p>
</td></tr>
<tr><td style="padding:20px 24px 24px;border-top:1px solid #dce5e9;">
<p style="margin:0 0 10px;color:#5d717a;font-size:13px;line-height:1.6;">This private link grants access to your task. Share it only with people you trust.</p>
<p style="margin:0 0 6px;color:#5d717a;font-size:12px;line-height:1.5;">Button not working? Copy this link into your browser:</p>
<p style="margin:0;font-size:12px;line-height:1.7;word-break:break-all;overflow-wrap:anywhere;"><a href="{href}" style="color:#117478;text-decoration:underline;word-break:break-all;overflow-wrap:anywhere;">{escape(link)}</a></p>
</td></tr></table>
<!--[if mso]></td></tr></table><![endif]-->
<p style="margin:20px 0 0;color:#5d717a;font-size:12px;line-height:1.7;">Jongmin Kim Lab</p>
</td></tr></table></body></html>'''
    message=EmailMessage()
    message['Subject']='Your ONT-LibQC report is ready' if ready else 'Your ONT-LibQC analysis needs attention'
    message['From']=smtp['from_address'];message['To']=row['email']
    message['Date']=formatdate(row['created'],usegmt=True)
    message['Message-ID']='<ont-libqc-'+jid+'-'+row['status']+'@'+urlsplit(origin).hostname+'>'
    if smtp.get('provider')=='resend':message['Resend-Idempotency-Key']='ont-libqc/'+jid+'/'+row['status']
    message.set_content(text);message.add_alternative(html,subtype='html')
    part=message.get_payload()[-1]
    logo=(Path(__file__).resolve().parents[1]/'static/mail-logo.png').read_bytes()
    part.add_related(logo,maintype='image',subtype='png',cid='<ontlibqc-logo>',filename='ontlibqc-logo.png',disposition='inline')
    # Stable boundaries keep repeated SMTP attempts identical for Resend.
    message.set_boundary('ontlibqc_alternative_'+jid)
    part.set_boundary('ontlibqc_related_'+jid)
    return message
