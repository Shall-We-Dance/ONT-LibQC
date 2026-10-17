[Portal guide](../README.md) / Email

# Email notifications

Email is optional. Without an address, users get results through their private task links. Those submissions do not enter mail queues or share an empty-address quota. IP/client and server-wide limits still apply. There is no email-based recovery endpoint.

A supplied address must be valid. Omitted, null, empty and whitespace-only values are stored as an empty string. All files in a batch use the same address, or all leave it blank.

## Enable delivery

Email is sent only when an address is provided **and** SMTP is enabled.

1. Edit `runtime/smtp.yaml` with the SMTP host, port, sender and credentials. Set `enabled: true` and use `starttls` (usually port 587) or `ssl` (usually 465). Plaintext SMTP and disabled certificate validation are not supported.
2. Keep this file mode `0600` and out of Git. Never put real credentials in source files or command-line arguments.
3. Set `mail_enabled: true` in `runtime/web.yaml` and in `runtime/frontend.yaml` for a split deployment, or `runtime/worker.yaml` for a single host. Recreate the web container.
4. Install the notification timer on the frontend using the deployment's Python environment:

```bash
python portal/scripts/mail.py --runtime /srv/ont-qc/runtime
```

The timer runs every 30 seconds with a 96 MB memory ceiling. In a split deployment, the frontend sends mail only after receiving and publishing the report; the compute bridge disables email for imported tasks. Jobs created with email disabled are not notified retroactively.

### Resend

Verify the sending domain in Resend first. A sending-only API key is sufficient.

| SMTP setting | Value |
| --- | --- |
| `host` | `smtp.resend.com` |
| `port` | `465` |
| `tls_mode` | `ssl` |
| `username` | `resend` |
| `password` | Your API key, stored only in the private runtime file |
| `provider` | `resend` |
| Sender | For example, `ONT-LibQC <reports@your-domain.org>` |

Disable open/click tracking for private report links. Keep Resend's SPF and DKIM records in DNS.

Publish one DMARC TXT record at `_dmarc` in the sending domain, starting with `v=DMARC1; p=none;`. For `reports@example.org`, the full name is `_dmarc.example.org`. This policy does not request quarantine or rejection; confirm authentication across legitimate senders before choosing a stricter policy. DMARC belongs in DNS, not `smtp.yaml`. See [Resend's DMARC guide](https://resend.com/docs/dashboard/domains/dmarc).

## One email per submission

The service waits until every file in the batch has reached a terminal state (completed, failed, cancelled or expired) and successful reports are published. It then sends one summary with:

- Each file's status, assay, reference, availability date and private task link.
- **View all results** for a batch, or **View your report / task** for a single file.
- A plain-text alternative and an HTML version with a small embedded logo. No remote images or sequencing files are attached.

Separate submissions stay separate, even if they share an email address. Batch metadata is saved when files are reserved, preventing early individual messages. Legacy tasks without batch metadata retain individual delivery; already sent messages are not sent again.

Sample labels are escaped and long filenames wrap. Subjects and body text use singular/plural wording appropriate to the file count, with a plain hyphen separator in batch subjects. Private links open the task page for report viewing and downloads.

## Retries and troubleshooting

Temporary failures retry after five minutes, for at most five attempts. Permanent rejections stop immediately. A process lock prevents competing senders; Resend idempotency keys suppress duplicates within its 24-hour window. Ambiguous deliveries older than 23 hours stop retrying. Generic SMTP cannot guarantee duplicate suppression after a connection failure.

Batch retries share one attempt counter and idempotency key. MIME boundaries remain stable, and the private MIME payload is frozen in the database during retries, then removed after delivery or expiration.

Provider acceptance does not guarantee inbox delivery. Check `mail_status` in the private task database and the service journal:

```bash
journalctl --user -u ont-libqc-mail.service
```

Logs omit recipients, links and credentials. Report deletion requires a second confirmation, and deleted or expired reports cannot be viewed or downloaded. The task page explains whether removal was manual or due to retention.

## Sender logos

The logo inside a message is separate from a mailbox app's sender avatar. Resend cannot set that avatar through SMTP.

- Apple Mail uses [Apple Branded Mail verification](https://resend.com/blog/apple-branded-mail).
- Gmail uses [BIMI with a supported mark certificate](https://knowledge.workspace.google.com/admin/security/set-up-bimi).
- `/static/brand-logo.png` and `/static/mail-logo.svg` are available for these setups. Publishing them, or adding DMARC alone, does not activate a verified inbox logo.
