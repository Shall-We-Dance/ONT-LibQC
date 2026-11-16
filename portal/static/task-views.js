/* Task card rendering; lifecycle actions are supplied by the controller. */
window.ONTTaskViews = function ({ node, assayLabel, labels, getCatalog, followReportLink, privateLink, saveText, api, refreshTask, renderTasks }) {
const taskViews = new WeakMap();
function taskState(task) {
  const job = task.job;
  if (job?.status === 'expired') return job.deletion_reason === 'manual' ? 'deleted' : 'expired';
  if (task.deletionRequested || job?.deletion_requested) return 'deleting';
  if (task.waitingForUpload) return 'waiting';
  if (task.uploading) return 'uploading';
  return job?.status || 'loading';
}
function taskProgress(task) {
  const job = task.job || {};
  const state = taskState(task);
  const percent = Math.round(Math.min(100, Math.max(0, (task.uploading ? task.uploadProgress * 100 : job.progress) || 0)));
  const message = state === 'deleting' ? 'Deletion requested. Report access has been removed.'
    : state === 'queued' ? (job.queue_ahead ? job.queue_ahead + (job.queue_ahead === 1 ? ' analysis ahead' : ' analyses ahead') : 'Waiting to start')
    : state === 'running' ? (job.stage || 'Starting analysis')
    : labels[state] || state;
  return { percent, message, state };
}
function taskView(task) {
  const job = task.job || { status: 'loading', sample: task.name || 'Loading analysis' };
  const expired = job.status === 'expired';
  const deleted = expired && job.deletion_reason === 'manual';
  const deletionPending = !expired && Boolean(task.deletionRequested || job.deletion_requested);
  const card = node('article', { class: 'task-card' });
  const state = taskState(task);
  const title = node('div', { class: 'task-title' }, node('h2', { class: !expired && job.sample ? 'task-sample' : '' }, expired ? (deleted ? 'Report deleted' : 'Report expired') : job.sample || 'Analysis'), node('span', { class: 'status ' + state }, labels[state] || state));
  card.append(title, node('p', { class: 'task-meta' }, [job.filename ? node('span', { class: 'task-filename' }, job.filename) : null, assayLabel(job.assay), job.reference].filter(Boolean).flatMap((part, index) => index ? [' · ', part] : [part])));
  const { percent, message } = taskProgress(task);
  const bar = node('span'); bar.style.width = percent + '%';
  if (expired) {
    const retention = getCatalog()?.retention_days;
    const reason = deleted
      ? 'This report was manually deleted. It can no longer be viewed or downloaded.'
      : 'This report was automatically deleted when its retention period ended.' + (retention ? ' Reports are kept for ' + retention + (retention === 1 ? ' day.' : ' days.') : '');
    card.append(node('div', { class: 'report-unavailable', role: 'status' }, node('p', {}, reason), node('p', {}, 'This cannot be undone. Submit the FASTQ again to generate a new report.')));
  } else if (deletionPending) {
    card.append(node('p', { class: 'progress-label deletion-pending' }, message));
  } else if (['created', 'uploading', 'waiting', 'uploaded', 'queued', 'running'].includes(state)) {
    card.append(node('div', { class: 'progress', role: 'progressbar', 'aria-label': state === 'uploading' || state === 'waiting' ? 'Upload progress' : 'Analysis progress', 'aria-valuemin': '0', 'aria-valuemax': '100', 'aria-valuenow': String(percent) }, bar));
    card.append(node('div', { class: 'task-progress-detail' }, node('p', { class: 'progress-label' }, message), node('span', { class: 'task-progress-value', hidden: !['uploading', 'running'].includes(state) }, percent + '%')));
  }
  if (task.failure || task.fetchError || job.error) card.append(node('div', { class: 'error' }, task.failure || task.fetchError || job.error));
  const actions = node('div', { class: 'task-actions' });
  if (job.status === 'succeeded' && !deletionPending && job.report_url && job.download_url) {
    const reportClick = event => followReportLink(event, task, 'report_url');
    const downloadClick = event => followReportLink(event, task, 'download_url');
    actions.append(node('a', { class: 'primary button', href: job.report_url, target: '_blank', rel: 'noopener', 'aria-label': 'Open ' + (job.sample || 'sample') + ' report', title: 'Opens in a new tab', onclick: reportClick, onauxclick: reportClick }, 'Open report'), node('a', { class: 'secondary button', href: job.download_url, 'aria-label': 'Download report ZIP', onclick: downloadClick, onauxclick: downloadClick }, 'Download'));
    const expires = new Date(job.expires * 1000);
    card.append(node('p', { class: 'field-note task-expiry' }, 'Available until ', node('time', { datetime: expires.toISOString() }, expires.toLocaleString())));
  }
  actions.append(node('button', { type: 'button', class: 'secondary task-utility copy-private-link', onclick: async event => { try { await navigator.clipboard.writeText(privateLink(task)); event.target.textContent = 'Link copied'; } catch { saveText('ont-libqc-private-link.txt', privateLink(task)); } } }, 'Copy link'));
  if (['created', 'uploaded', 'queued', 'running'].includes(job.status) && !task.uploading) actions.append(node('button', { type: 'button', class: 'secondary task-utility cancel-analysis', onclick: async () => { try { await api('/api/jobs/' + task.id + '/cancel', 'POST', {}, task.token); await refreshTask(task); } catch (e) { task.failure = e.message; renderTasks(); } } }, 'Cancel analysis'));
  if (['succeeded', 'failed', 'cancelled'].includes(job.status)) {
    const warning = node('p', { class: 'deletion-warning', id: 'delete-warning-' + task.id, role: 'status', hidden: true }, 'Permanently deletes this report and its files. This cannot be undone.');
    const deleteButton = node('button', { type: 'button', class: 'secondary task-utility delete-report', onclick: async () => {
      if (task.deleting || task.deletionRequested || job.deletion_requested) return;
      if (!task.confirmDelete) { task.confirmDelete = true; updateDeleteButton(); return; }
      task.deleting = true; task.failure = null; updateDeleteButton();
      try {
        await api('/api/jobs/' + task.id, 'DELETE', undefined, task.token);
        task.deletionRequested = true; task.confirmDelete = false;
      } catch (e) { task.failure = e.message; }
      finally { task.deleting = false; renderTasks(); }
    } });
    function updateDeleteButton() {
      const requested = Boolean(task.deletionRequested || job.deletion_requested);
      const confirming = Boolean(task.confirmDelete && !requested);
      deleteButton.classList.toggle('danger-confirm', confirming);
      deleteButton.textContent = task.deleting ? 'Requesting deletion…' : requested ? 'Deletion requested' : confirming ? 'Confirm permanent deletion' : 'Delete report';
      deleteButton.disabled = Boolean(task.deleting || requested);
      warning.hidden = !confirming;
      if (confirming) deleteButton.setAttribute('aria-describedby', warning.id);
      else deleteButton.removeAttribute('aria-describedby');
    }
    updateDeleteButton(); actions.append(deleteButton, warning);
  }
  card.append(actions);
  if (job.options && job.status !== 'expired') { const detail = node('details', {}, node('summary', {}, 'Submitted settings'), node('pre', { class: 'field-note' }, JSON.stringify(job.options, null, 2))); detail.open = Boolean(task.detailsOpen); detail.addEventListener('toggle', () => { task.detailsOpen = detail.open; }); card.append(detail); }
  return card;
}
// Keep existing elements connected so polling does not restart motion or lose focus.
function syncChildren(parent, children) {
  const keep = new Set(children);
  for (const child of [...parent.children]) if (!keep.has(child)) child.remove();
  children.forEach((child, index) => {
    if (parent.children[index] !== child) parent.insertBefore(child, parent.children[index] || null);
  });
}
function updateTaskView(task) {
  // Progress and timestamps change often; update the bar without rebuilding controls.
  const { progress, updated, stage, ...job } = task.job || {};
  const key = JSON.stringify([job, task.name, task.uploading, task.waitingForUpload, task.failure, task.fetchError, task.confirmDelete, task.deleting, task.deletionRequested, getCatalog()?.retention_days]);
  let view = taskViews.get(task);
  if (!view || view.key !== key) {
    const card = taskView(task);
    if (view) {
      const previousProgress = view.card.querySelector('.progress');
      const nextProgress = card.querySelector('.progress');
      syncChildren(view.card, [...card.children].map(child => child === nextProgress && previousProgress ? previousProgress : child));
      view.key = key;
    } else {
      view = { card, key }; taskViews.set(task, view);
    }
  }
  const { percent, message } = taskProgress(task);
  const track = view.card.querySelector('.progress');
  const phase = task.uploading ? 'upload' : 'analysis';
  if (track) {
    // Upload and analysis have separate scales; reset immediately between them.
    track.firstElementChild.style.transitionDuration = view.phase === phase ? '' : '0s';
    track.firstElementChild.style.width = percent + '%';
    track.setAttribute('aria-valuenow', String(percent));
    track.setAttribute('aria-label', task.uploading || task.waitingForUpload ? 'Upload progress' : 'Analysis progress');
  }
  const label = view.card.querySelector('.progress-label');
  if (label && label.textContent !== message) label.textContent = message;
  const value = view.card.querySelector('.task-progress-value');
  if (value && value.textContent !== percent + '%') value.textContent = percent + '%';
  view.phase = phase;
  return view.card;
}
return { taskState, syncChildren, updateTaskView };
};
