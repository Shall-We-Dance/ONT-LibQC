[Portal guide](../README.md) / Deployment

# Deploy the portal

Run commands from the repository root as the deployment owner. Paths and domains below are examples; keep runtime files and secrets outside the source tree.

- [Install](#install)
- [Domain HTTPS and renewal](#domain-https-and-renewal)
- [Operations and migration](#operations-and-migration)
- [Separate frontend and compute hosts](#separate-frontend-and-compute-hosts)

## Install

Requirements: Docker Engine and Compose, Python 3.11 or later with venv support, OpenSSL, and systemd user services. The dispatcher runs as an ordinary user with Docker access. Use a dedicated service account on a shared production server. Do not mount the Docker socket into the public web container.

Run the following from the repository root:

```bash
python3 -m venv /srv/ont-qc/venv
/srv/ont-qc/venv/bin/pip install -r portal/requirements.txt

docker build -f portal/deploy/Dockerfile.runner -t ont-qc-runner:1 .
docker build -f portal/deploy/Dockerfile.web -t ont-qc-web:1 .
```

Choose a runtime directory writable by the service account; `/srv/ont-qc` is an example. Keep it outside the source tree. Building the runner installs the pipeline's aligners, fastp, Pairtools, Snakemake, plotting tools, R and karyoploteR. Rebuild the runner image when upgrading a deployment that predates the Pairtools analysis, then update `runner_image` in the compute worker settings and restart the idle worker.

Prepare or adopt [reference indexes](references.md), then initialize the deployment:

```bash
/srv/ont-qc/venv/bin/python portal/scripts/service.py init \
  --runtime /srv/ont-qc/runtime \
  --catalog /srv/ont-qc/references/common/catalog.admin.json \
  --public-url https://qc.example.org:8443 --port 8443 \
  --cert /path/to/fullchain.pem --key /path/to/privkey.pem

/srv/ont-qc/venv/bin/python portal/scripts/service.py up \
  --runtime /srv/ont-qc/runtime
loginctl enable-linger "$(id -un)"
```

The initializer generates an optional submission code in `runtime/secrets.json`. Share that code with invited users. Add `--open-submissions` to initialization if you want bounded public submissions without a code. Private report access is always enforced independently of the submission code.

### Evaluation without a domain

With no domain or certificate, use an HTTPS IP origin, such as `https://192.0.2.10:8443`, and omit `--cert`/`--key`. This creates a **30-day self-signed certificate for initial evaluation**. Browsers will not trust it automatically. Before general external use, install a certificate trusted by your users; do not ask users to disable certificate verification. Distribute a self-signed certificate only through a separate trusted channel if it must be used within a controlled group. Certificate renewal/replacement is the operator's responsibility. After replacing certificates, run `service.py harden` to refresh the backend CA bundle, then recreate both web and edge containers.

### Network boundary

Only the Nginx edge container publishes a port. It handles TLS, limits requests and connections, and closes incomplete headers after 10 seconds. The web application is on a separate internal Docker network with no published port. Edge-to-web traffic uses verified TLS. The edge overwrites forwarded client headers using its socket peer; the backend trusts these headers only because this network contains the trusted edge and web containers. Do not attach untrusted containers or publish the backend port. Reports, credentials and upload bodies are not logged by the proxy.

### Upgrade an existing deployment

Existing deployments can adopt the same layout without moving data or rebuilding references:

```bash
python portal/scripts/service.py harden --runtime /srv/ont-qc/runtime
docker compose -f /srv/ont-qc/runtime/compose.yaml up -d --force-recreate web edge
```

Build the updated web image before applying a source upgrade and retain the existing signing key. `harden` preserves the configured web image, public port and bind address, writes `edge.conf` and a certificate trust bundle, and configures the small `nginx:stable-alpine` container. Run it again after changing limits or certificates. Public TLS and internal backend TLS use separate certificates; the latter has a fixed verified DNS identity even when the public address is an IP. Apply web updates when no files are actively uploading; running analysis containers are independent.

### Private network and address changes

For access through ZeroTier, point the domain's DNS record to the server's ZeroTier address and add `--bind <server-zerotier-ip>` when initializing. Use the domain in `--public-url`, with the chosen HTTPS port. Clients must belong to the authorized ZeroTier network and be able to reach that address. Binding to the ZeroTier address limits the published port to that interface; certificate trust is still required.

To change an existing deployment's address, update `public_url` and `allowed_hosts` in both `web.yaml` and `worker.yaml`, update the host address in `compose.yaml` under `ports`, and install a certificate covering the new hostname. Keep `secrets.json` unchanged to preserve private task credentials. Apply the change when no jobs are active with `docker compose -f /srv/ont-qc/runtime/compose.yaml up -d --force-recreate web edge` and `systemctl --user restart ont-qc-portal-worker.service`. Existing task links need their origin changed to the new URL while retaining the path and fragment.

HTTPS enforcement uses `require_https` (default `true`). Legacy `secure_cookies` configurations are read as this setting during migration; they do not enable cookies.

## Domain HTTPS and renewal

Install Certbot on the host and point the domain's A records at the frontend. Open inbound TCP 80 and 443. Port 80 serves only ACME validation and HTTPS redirects. With Cloudflare proxying enabled, use **Full (strict)**; keep HTTP ACME challenge paths exempt from any custom redirect/WAF challenge rules. No Cloudflare API key is required for HTTP validation.

Run as the deployment owner (not root), using the deployment's Python environment:

```bash
python portal/scripts/https.py prepare --runtime /srv/ont-qc/runtime \
  --domain example.org --alias www.example.org --cloudflare --legacy-port 8443
# Add --docker-sudo above only when this host requires sudo for Docker.
docker compose -f /srv/ont-qc/runtime/compose.yaml up -d --force-recreate web edge
python portal/scripts/https.py issue --runtime /srv/ont-qc/runtime
python portal/scripts/https.py renew-timer --runtime /srv/ont-qc/runtime
```

After issuance, set `public_url: https://example.org` and the domain in `allowed_hosts` in `web.yaml` and `frontend.yaml` (or `worker.yaml` on a single host); recreate web and edge. Preserve `secrets.json` and `state/` so private task links remain valid. The optional old IP port redirects to the domain, retaining the old IP certificate. Omit `--legacy-port` for a fresh domain deployment and omit `--cloudflare` when using DNS-only/direct hosting.

The Cloudflare option retrieves its official IP ranges, trusts `CF-Connecting-IP` only from those networks, and blocks direct access to the domain's application listener. Real client addresses still drive per-source admission limits. Re-run `prepare` and reload the edge when updating the trusted ranges. HTTP validation stays reachable from certificate authorities. Application and report responses use `Cache-Control: no-store`; do not add Cloudflare Cache Everything rules for the portal.

Certificates, ACME credentials and logs live in the private runtime directory. The public web container receives only its internal backend TLS key, not the public certificate's private key or SMTP credentials. The user timer checks renewal twice daily, reloads nginx after successful renewal, and survives logout/reboot when user lingering is enabled. To verify renewal independently:

```bash
certbot renew --dry-run --config-dir /srv/ont-qc/runtime/acme \
  --work-dir /srv/ont-qc/runtime/acme-work --logs-dir /srv/ont-qc/runtime/acme-logs
systemctl --user list-timers ont-libqc-certificate.timer
```

Keep TCP 80 open for HTTP validation, or migrate renewal to a DNS challenge plugin if port 80 must be closed. After moving to another host, update DNS and runtime paths and reinstall the timers; do not copy old host paths blindly.

## Operations and migration

```bash
python portal/scripts/service.py status --runtime /srv/ont-qc/runtime
journalctl --user -u ont-qc-portal-worker.service --since today
python portal/scripts/service.py down --runtime /srv/ont-qc/runtime
```

`down` stops the services without deleting runtime data. The worker service and Docker restart policy survive logout; lingering must be enabled for automatic user-service startup after reboot. Docker Engine itself must be configured to start at boot.

To migrate to another Linux host:

1. Copy the repository and prepare Docker/Python as above. Build the two images, or transfer them using `docker save` / `docker load` for compatible CPU architectures.
2. Copy or rebuild reference files and indexes. Run `adopt` against paths valid on the destination host. Source code contains no server-specific reference paths.
3. For a fresh service, run `service.py init` with the destination URL, runtime directory and certificates.
4. To retain existing reports, stop the old service first, initialize destination settings, then copy its `state/` directory and **original `secrets.json`** to the destination runtime directory. Preserve the service user's file ownership and restrictive permissions. The signing key is required for existing private links. Do not overwrite the newly generated destination configuration with old machine paths. If the hostname changes, update the origin in links given to users.
5. Start the destination services, verify `/healthz` and a small job, and then route users to the destination.

Back up private reports and the signing key only if your retention policy permits it. Encrypt and restrict any such backups; do not include them in source control. The portable deliverables are the source, image definitions, reference preparation script and deployment configuration templates.

## Separate frontend and compute hosts

A small frontend VM can host the public HTTPS UI, upload storage, task database and completed report downloads. The compute host retains all references and runs the same bounded dispatcher for both local and remote submissions. No reference or analysis image needs to be copied to the frontend.

On the frontend, copy this source tree and the compute host's **public** reference catalog, prepare the Python venv, and load/build the web image. Initialize the service with a frontend-only configuration:

```bash
/srv/ont-qc/venv/bin/python portal/scripts/service.py init \
  --runtime /srv/ont-qc/runtime --catalog /srv/ont-qc/catalog.public.json \
  --public-url https://203.0.113.10:8443 --port 8443 --bind 0.0.0.0 \
  --frontend-only --open-submissions
/srv/ont-qc/venv/bin/python portal/scripts/service.py up --runtime /srv/ont-qc/runtime
```

The example IP is a placeholder. Omit certificate arguments only for the intended self-signed evaluation deployment. Allow TCP 8443 in the cloud firewall only for the frontend instance. Keep SSH restricted to the compute host or administrator networks where possible. A frontend without a connected compute dispatcher refuses new submissions and continues serving existing task pages and reports.

On the compute host, configure an existing trusted SSH alias for the frontend, then install a dedicated restricted transfer key and the bridge:

```bash
python portal/scripts/connect_frontend.py --frontend frontend-ssh-alias \
  --runtime /srv/ont-qc/runtime --remote-root /srv/ont-qc
systemctl --user link /srv/ont-qc/runtime/ont-qc-remote-bridge.service
systemctl --user daemon-reload
systemctl --user enable --now ont-qc-remote-bridge.service
```

Run the connector using the compute worker's venv. The remote root must contain `source/`, `venv/` and `runtime/`. It installs a forced-command SSH key with shell, PTY and forwarding disabled. Private keys stay on the compute host; the public web container receives no SSH credentials. Enable user-service lingering on both hosts so the bridge and frontend retention timer survive logout and reboot. Use one bridge per compute runtime. Existing signing keys and private report links are preserved by ordinary code upgrades.

Transfers use encrypted SSH streams with persisted byte offsets. Interrupted FASTQ or BED downloads remain staged and are resumed before entering the compute queue. Reports are also resumed from the last received byte. The frontend validates archive paths, member types, size limits and report links, then publishes a complete report atomically. Repeating a request after a lost acknowledgement does not rerun analysis or duplicate a published report. Partial files are never exposed through report URLs.

An SSH interruption does not stop an analysis already running on the compute host. Its completed report is retained there until the frontend acknowledges its own durable copy; the usual frontend report retention period starts at publication. Outstanding transfers remain subject to queue and storage limits. If the connection stays down and storage fills, further analyses wait instead of exhausting disk space. Reconnection also reconciles progress, cancellations and expired/deleted tasks. The independent frontend cleanup timer removes expired reports even when the compute host is offline.

After changing limits, update the limits files on both hosts, run `service.py harden` on each frontend, and recreate its web/edge services. Restart the compute worker only after active analyses finish; changing a batch/upload limit alone does not require stopping a running analysis.
