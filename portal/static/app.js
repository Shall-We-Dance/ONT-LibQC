'use strict';
const $ = id => document.getElementById(id);
let catalog, files = [], tasks = [], submitting = false, refreshing = false, retryRequestsAt = 0;
// These IDs only connect controls within this page; task credentials come from the server.
let nextLocalId = 0;
const taskViews = new WeakMap();
const bytes = n => (n / 1_000_000).toLocaleString(undefined, { maximumFractionDigits: 2 }) + ' MB';
function node(tag, attrs = {}, ...children) {
  const n = document.createElement(tag);
  for (const [k, v] of Object.entries(attrs)) {
    if (k.startsWith('on')) n.addEventListener(k.slice(2), v);
    else if (['checked', 'disabled', 'hidden', 'value', 'required'].includes(k)) n[k] = v;
    else n.setAttribute(k, v);
  }
  for (const c of children.flat()) if (c !== null && c !== undefined) n.append(c instanceof Node ? c : document.createTextNode(String(c)));
  return n;
}
function showError(message, target = 'form-error') { $(target).textContent = message; $(target).hidden = !message; }
// Do not infer file type from a .bed suffix: directories may have that suffix too.
async function checkedBedFile(input, maxBytes = 10 * 1024 ** 2) {
  const file = input.files?.[0];
  if ([...(input.webkitEntries || [])].some(entry => entry.isDirectory) || ['inode/directory', 'application/x-directory'].includes(file?.type)) {
    throw Error('Choose a BED file inside the folder, not the folder itself.');
  }
  if (!file) throw Error('Choose a BED file first.');
  if (file.size > maxBytes) throw Error('Use a BED file of at most ' + maxBytes / 1024 ** 2 + ' MiB.');
  // Dialog selections often have no filesystem entry. A read probe rejects
  // directory-like or unreadable selections before parsing or upload.
  try { await file.slice(0, 1).arrayBuffer(); }
  catch { throw Error('Cannot read this selection. Open the folder and choose a readable BED file inside it.'); }
  if (!file.size) throw Error('The BED file is empty. Choose a file containing regions.');
  return file;
}
function bedDropIsDirectory(event) {
  return [...(event.dataTransfer?.items || [])].some(item => (item.getAsEntry?.() || item.webkitGetAsEntry?.())?.isDirectory);
}
async function api(path, method = 'GET', body, token) {
  const headers = { 'X-Portal-Request': '1' };
  if (token) headers.Authorization = 'Bearer ' + token;
  if (body !== undefined) headers['Content-Type'] = 'application/json';
  const response = await fetch(path, { method, headers, body: body === undefined ? undefined : JSON.stringify(body), credentials: 'same-origin' });
  const data = await response.json().catch(() => ({}));
  if (!response.ok) {
    if ([429, 503].includes(response.status)) retryRequestsAt = Date.now() + (Math.min(120, Math.max(5, Number(response.headers.get('Retry-After')) || 5)) + Math.random() * 3) * 1000;
    throw new Error(typeof data.detail === 'string' ? data.detail : 'The service is busy. Please try again shortly.');
  }
  return data;
}
function control(label, input) { const id = 'field-' + ++nextLocalId; input.id = id; return node('div', { class: 'field' }, node('label', { for: id }, label), input); }
function select(values, value) { const input = node('select', {}, values.map(([id, label]) => node('option', { value: id }, label))); input.value = value; return input; }
function requiredChoice(label, input) {
  input.required = true;
  const field = control(label, input);
  field.querySelector('label').append(node('span', { class: 'required-marker' }, 'Required'));
  field.classList.add('required-choice');
  const hint = node('p', { class: 'visually-hidden', id: input.id + '-hint' }, 'Please choose for this file.');
  input.setAttribute('aria-describedby', hint.id);
  const update = () => { field.classList.toggle('needs-choice', !input.value); hint.textContent = input.value ? 'Selected for this file.' : 'Please choose for this file.'; };
  input.addEventListener('change', update); update(); field.append(hint); return field;
}
function check(label, description = '') {
  const input = node('input', { type: 'checkbox' });
  return { input, wrapper: node('label', { class: 'check' }, input, node('span', {}, label, description ? node('em', {}, description) : null)) };
}
function screenCheck(label, description) {
  const field = check(label, description);
  field.preferred = true;
  field.input.addEventListener('change', () => { field.preferred = field.input.checked; });
  return field;
}
function referenceChange(entry) {
  const r = catalog.references.find(x => x.id === entry.ref.value) || { assembly: 'the selected reference', primary_chromosomes: [] };
  entry.bedAssembly = r.assembly || r.id;
  entry.bedConfirmText.textContent = 'I confirm this BED uses ' + entry.bedAssembly + ' coordinates';
  entry.bedHint.textContent = 'Choose a BED file for ' + entry.bedAssembly + '. Up to 4 MiB. The metaplot will be included in your report.';
  entry.bedConfirm.checked = false;
  entry.myco.input.disabled = !r.mycoplasma;
  entry.spike.input.disabled = !r.spike_in;
  entry.myco.input.checked = Boolean(r.mycoplasma && entry.myco.preferred);
  entry.spike.input.checked = Boolean(r.spike_in && entry.spike.preferred);
  entry.chromosomes.placeholder = r.primary_chromosomes.join(', ');
  entry.chromosomes.value = '';
  updateAssay(entry);
}
function updateAssay(entry) {
  const reference = catalog.references.find(x => x.id === entry.ref.value);
  const r = reference || { mappers: [] };
  const microc = entry.assay.querySelector('option[value="microc"]');
  microc.disabled = Boolean(reference && !r.mappers.some(mapper => ['bwa-mem2', 'bwa-mem', 'minimap2'].includes(mapper)));
  microc.textContent = microc.disabled ? 'Micro-C · unavailable on this server' : 'Micro-C';
  entry.availabilityNote.textContent = reference?.availability_note || '';
  entry.availabilityNote.hidden = entry.assay.value !== 'microc' || !reference?.availability_note;
  if (entry.assay.selectedOptions[0]?.disabled) {
    entry.assay.value = '';
    entry.assay.dispatchEvent(new Event('change'));
    return;
  }
  const assay = entry.assay.value;
  entry.mapper.disabled = !assay || !entry.ref.value;
  entry.mapper.replaceChildren(node('option', { value: 'auto' }, !assay ? 'Choose a library assay first' : assay === 'microc' ? `Default · ${(r.microc_mapper || 'bwa-mem2').toUpperCase()}` : 'Default · Bowtie2'),
    ...r.mappers.filter(x => assay !== 'microc' || x !== 'bowtie2').map(x => node('option', { value: x }, x)));
  entry.bedCheck.input.disabled = assay !== 'chromatin';
  if (assay !== 'chromatin') entry.bedCheck.input.checked = false;
  updateCustomRegions(entry);
}
function updateCustomRegions(entry) {
  const enabled = entry.assay.value === 'chromatin' && entry.bedCheck.input.checked;
  entry.bedBox.hidden = !enabled;
  for (const input of [entry.bed, entry.bedLabel, entry.bedMode, entry.bedConfirm]) input.disabled = !enabled;
  entry.bedConfirm.required = enabled;
  if (!enabled) entry.bedConfirm.checked = false;
}
function customField(title, input) {
  return node('label', { class: 'custom-field' }, node('span', { class: 'custom-field-title' }, title), input);
}
function filenameProblem(name) {
  const suffix = ['.fastq.gz', '.fq.gz', '.fastq', '.fq'].find(ext => name.toLowerCase().endsWith(ext));
  if (!suffix) return '“' + name + '” was not added. Choose a .fastq, .fq, .fastq.gz or .fq.gz file. POD5 and FAST5 files must be basecalled first; changing their extension does not convert them to FASTQ.';
  const problems = [];
  if (name.length > 180) problems.push('the filename exceeds 180 characters');
  if (!/^[A-Za-z0-9]/.test(name)) problems.push('the filename must start with an English letter or number');
  const unsupported = [...new Set(name.match(/[^A-Za-z0-9_. -]/gu) || [])];
  if (unsupported.length) problems.push('unsupported characters: ' + unsupported.map(c => JSON.stringify(c)).join(', '));
  if (!problems.length) return '';
  let stem = name.slice(0, -suffix.length).replace(/[^A-Za-z0-9_.-]+/g, '_').replace(/_+/g, '_').replace(/[._-]+$/, '');
  if (!/^[A-Za-z0-9]/.test(stem)) stem = 'sample_' + stem.replace(/^[^A-Za-z0-9]+/, '');
  stem = stem.slice(0, 180 - suffix.length).replace(/[._-]+$/, '') || 'sample';
  return '“' + name + '” was not added: ' + problems.join('; ') + '. Use English letters, numbers, spaces, dots, hyphens (-) or underscores (_). Rename it to “' + stem + name.slice(-suffix.length) + '” and try again.';
}
function showUploadErrors(messages) {
  const box = $('upload-error');
  box.replaceChildren(...messages.map(message => node('p', {}, message)));
  box.hidden = !messages.length;
}
function addFiles(incoming) {
  showError('');
  const errors = [];
  for (const file of incoming) {
    if (files.length >= catalog.max_files) { errors.push('Some files were not added. Choose up to ' + catalog.max_files + ' files per submission.'); break; }
    const problem = filenameProblem(file.name);
    if (problem) { errors.push(problem); continue; }
    if (!file.size) { errors.push('“' + file.name + '” was not added because it is empty. Choose a FASTQ file containing reads.'); continue; }
    if (file.size > catalog.max_upload_bytes) { errors.push('“' + file.name + '” was not added: maximum ' + bytes(catalog.max_upload_bytes) + ' per file. For larger files, download the pipeline from GitHub and run it locally.'); continue; }
    const entry = { file, id: ++nextLocalId };
    entry.sample = node('input', { class: 'sample-name', value: file.name.replace(/\.(fastq|fq)(\.gz)?$/i, '').replace(/[^A-Za-z0-9_-]/g, '_').replace(/^[^A-Za-z]/, 'S_').slice(0, 48), required: true, pattern: '[A-Za-z][A-Za-z0-9_-]{0,47}', maxlength: '48' });
    entry.assay = select([['', 'Choose library assay…'], ['chromatin', 'Chromatin · CUT&RUN / CUT&Tag / ChIP / ATAC, etc.'], ['microc', 'Micro-C']], '');
    const available = catalog.references.filter(r => r.ready);
    entry.ref = select([['', 'Choose reference genome…'], ...available.map(r => [r.id, r.label])], '');
    entry.mapper = select([['auto', 'Default · Bowtie2']], 'auto');
    entry.myco = screenCheck('Mycoplasma screen', 'An opportunistic screen; negative results do not exclude contamination.');
    entry.spike = screenCheck('Spike-in / foreign-sequence screen', 'Use the host-competitive reference panel available for this assembly.');
    entry.declared = check('A spike-in was intentionally added');
    entry.bedCheck = check('Include custom BED regions', 'Optional. Chromatin reports also support BED uploads in the browser.');
    entry.bed = node('input', { type: 'file' }); // Keep .bed folders navigable in Safari.
    entry.bedLabel = node('input', { value: 'Uploaded regions', maxlength: '80', pattern: '[A-Za-z0-9 _().-]{1,80}' });
    entry.bedMode = select([['auto', 'Auto'], ['reference_point', 'Reference point'], ['region', 'Region (scaled)']], 'auto');
    entry.bedConfirm = node('input', { type: 'checkbox' });
    entry.bedConfirmText = node('span');
    entry.bedHint = node('p', { class: 'field-note' });
    entry.bedName = node('span', { class: 'custom-file-name', 'aria-hidden': 'true' }, 'No BED selected');
    const bedPicker = node('span', { class: 'custom-file-picker' }, entry.bed,
      node('span', { class: 'custom-file-button', 'aria-hidden': 'true' }, 'Choose BED'), entry.bedName);
    entry.bedBox = node('div', { class: 'bed-input', hidden: true },
      node('h3', {}, 'Custom region metaplot'), entry.bedHint,
      node('div', { class: 'custom-controls' }, customField('Region BED', bedPicker),
        customField('Region label', entry.bedLabel),
        customField('Plot mode', node('span', { class: 'custom-select-shell' }, entry.bedMode)),
        node('div', { class: 'custom-actions' }, node('label', { class: 'custom-confirm' }, entry.bedConfirm, entry.bedConfirmText))));
    entry.bedError = node('p', { class: 'error', role: 'alert', hidden: true }); entry.bedBox.append(entry.bedError);
    const bedSelection = async () => {
      const version = entry.bedVersion = (entry.bedVersion || 0) + 1;
      const file = entry.bed.files[0], name = file?.name || 'No BED selected';
      entry.bedName.textContent = name; bedPicker.title = name; entry.bedConfirm.checked = false;
      entry.bedError.hidden = true; entry.bedError.textContent = '';
      entry.bedChecking = Boolean(file); updateCount();
      if (!file) return;
      try { await checkedBedFile(entry.bed, 4 * 1024 ** 2); }
      catch (e) {
        if (version === entry.bedVersion) {
          entry.bed.value = ''; entry.bedName.textContent = 'No BED selected'; bedPicker.title = '';
          entry.bedError.textContent = e.message; entry.bedError.hidden = false;
        }
      } finally { if (version === entry.bedVersion) { entry.bedChecking = false; updateCount(); } }
    };
    entry.bed.addEventListener('change', bedSelection);
    entry.bed.addEventListener('drop', event => {
      if (!bedDropIsDirectory(event)) return;
      event.preventDefault(); entry.bedVersion = (entry.bedVersion || 0) + 1; entry.bed.value = ''; entry.bedChecking = false;
      entry.bedName.textContent = 'No BED selected'; bedPicker.title = ''; entry.bedConfirm.checked = false;
      entry.bedError.textContent = 'Choose a BED file inside the folder, not the folder itself.'; entry.bedError.hidden = false; updateCount();
    });
    entry.bedCheck.input.addEventListener('change', () => { updateCustomRegions(entry); updateCount(); });
    entry.pairs = node('textarea', { rows: '3', placeholder: 'libA ACGTACGT TGCATGCA\nlibB GATCTAGC CTAGATCG' });
    entry.chromosomes = node('input', { placeholder: '' });
    const advanced = node('details', {}, node('summary', {}, 'Advanced settings'));
    const info = node('p', { class: 'field-note' }, 'Defaults match the ONT-native pipeline. Leave fields unchanged to use the standard settings. Rates are fractions: 0.8 means 80%.');
    advanced.append(info, node('div', { class: 'fields' }, control('Alignment method', entry.mapper), control('Primary chromosomes (optional override)', entry.chromosomes)),
      node('p', { class: 'advanced-group' }, 'Expected physical index pairs'), entry.pairs,
      node('p', { class: 'field-note' }, 'One pair per line: name i5 i7. Leave empty to discover index pairs automatically.'));
    entry.params = {};
    for (const group of [...new Set(catalog.params.map(p => p.group))]) {
      const grid = node('div', { class: 'fields' });
      for (const p of catalog.params.filter(x => x.group === group)) {
        const input = node('input', { type: 'number', min: String(p.min), max: String(p.max), step: p.kind === 'int' ? '1' : 'any', placeholder: String(p.default), value: '' });
        entry.params[p.key] = input;
        const field = control(p.label, input); field.insertBefore(node('span', { class: 'param-key' }, p.key), input); grid.append(field);
      }
      advanced.append(node('p', { class: 'advanced-group' }, group), grid);
    }
    entry.screen = {};
    const screenGrid = node('div', { class: 'fields' });
    for (const [key, label, value, min, max] of [
      ['screen_minimum_reads', 'Minimum screened reads', 500, 1, 1000000],
      ['spike_minimum_aligned_bp', 'Spike-in minimum aligned bp', 80, 20, 1000],
      ['spike_minimum_query_coverage', 'Spike-in minimum coverage', 0.8, 0.1, 1],
      ['spike_minimum_identity', 'Spike-in minimum identity', 0.85, 0.1, 1]]) {
      const input = node('input', { type: 'number', value: String(value), min: String(min), max: String(max), step: max === 1 ? 'any' : '1' });
      entry.screen[key] = input; screenGrid.append(control(label, input));
    }
    advanced.append(node('p', { class: 'advanced-group' }, 'Screening evidence'), screenGrid);
    const ranges = node('details', {}, node('summary', {}, 'Report reference cutoffs'), node('p', { class: 'field-note' }, 'Edit local screening cutoffs. Check “Descriptive only” to remove the cutoff for a metric. These are not universal biological pass/fail criteria.'));
    entry.ranges = {};
    for (const [key, defaults] of Object.entries(catalog.ranges)) {
      const disable = check('Descriptive only'); disable.input.checked = defaults === null;
      const row = node('div', { class: 'range-row' }, node('span', {}, key));
      const values = {};
      for (const limit of ['min', 'max', 'warning_min', 'warning_max']) {
        values[limit] = node('input', { type: 'number', step: 'any', min: '0', placeholder: limit, 'aria-label': key + ' ' + limit, value: defaults?.[limit] === undefined ? '' : String(defaults[limit]), disabled: defaults === null });
        row.append(values[limit]);
      }
      disable.input.addEventListener('change', () => { for (const input of Object.values(values)) input.disabled = disable.input.checked; });
      ranges.append(row, disable.wrapper); entry.ranges[key] = { values, disable: disable.input };
    }
    const header = node('div', { class: 'file-header' }, node('div', {}, node('strong', {}, file.name), node('br'), node('small', {}, bytes(file.size))),
      node('button', { type: 'button', class: 'remove', 'aria-label': 'Remove ' + file.name, onclick: () => { files = files.filter(f => f.id !== entry.id); entry.card.remove(); updateCount(); } }, '×'));
    entry.structure = LibraryStructure.create(catalog.library_structure, value => applyStructureToAll(entry, value));
    entry.applyNote = node('span', { class: 'apply-note', role: 'status', 'aria-live': 'polite' });
    entry.availabilityNote = node('p', { class: 'field-note', role: 'status', hidden: true });
    entry.applyButton = node('button', { type: 'button', class: 'secondary apply-all', onclick: () => applyToAll(entry), 'aria-label': 'Apply reference genome and library assay from ' + file.name + ' to all libraries' }, 'Apply to all');
    entry.applyRow = node('div', { class: 'apply-row', hidden: true }, entry.applyNote, entry.applyButton);
    entry.card = node('article', { class: 'file-card' }, header, node('div', { class: 'file-body' },
      node('div', { class: 'fields' }, control('Sample name', entry.sample), requiredChoice('Reference genome', entry.ref), node('div', { class: 'full' }, requiredChoice('Library assay', entry.assay))),
      entry.availabilityNote, entry.applyRow, entry.structure.root, node('div', { class: 'check-list' }, entry.myco.wrapper, entry.spike.wrapper, entry.declared.wrapper, entry.bedCheck.wrapper), entry.bedBox, advanced, ranges));
    entry.ref.addEventListener('change', () => referenceChange(entry)); entry.assay.addEventListener('change', () => updateAssay(entry));
    for (const input of [entry.ref, entry.assay]) input.addEventListener('change', () => {
      for (const file of files) file.applyFeedback = '';
      updateApplyControls();
    });
    referenceChange(entry); files.push(entry); $('file-list').append(entry.card);
  }
  showUploadErrors(errors);
  updateCount();
}
function updateApplyControls() {
  for (const entry of files) {
    entry.applyRow.hidden = files.length < 2;
    entry.structure.applyButton.hidden = files.length < 2;
    entry.structure.applyButton.disabled = submitting;
    entry.applyButton.disabled = submitting || !entry.ref.value || !entry.assay.value;
    entry.applyNote.textContent = !entry.ref.value || !entry.assay.value
      ? 'Choose both fields to apply to all libraries.'
      : entry.applyFeedback || 'Apply this reference and assay to all libraries.';
  }
}
function applyStructureToAll(source, value) {
  if (submitting) return;
  for (const entry of files) if (entry !== source) entry.structure.setValue(value);
  source.structure.feedback('Structure applied to ' + (files.length - 1) + ' other ' + (files.length === 2 ? 'library.' : 'libraries.'));
}
function applyToAll(source) {
  if (submitting || !source.ref.value || !source.assay.value) return;
  for (const entry of files) {
    entry.applyFeedback = '';
    if (entry === source) continue;
    if (entry.ref.value !== source.ref.value) {
      entry.ref.value = source.ref.value;
      entry.ref.dispatchEvent(new Event('change'));
    }
    if (entry.assay.value !== source.assay.value) {
      entry.assay.value = source.assay.value;
      entry.assay.dispatchEvent(new Event('change'));
    }
  }
  const count = files.length - 1;
  source.applyFeedback = 'Applied to ' + count + ' other ' + (count === 1 ? 'library.' : 'libraries.');
  updateApplyControls();
}
function updateCount() {
  updateApplyControls();
  $('file-count').textContent = files.length + (files.length === 1 ? ' file' : ' files');
  $('empty-hint').hidden = files.length > 0;
  $('submit-button').disabled = !files.length || submitting || !catalog?.worker_online || !$('data-policy').checked || files.some(e => e.bedCheck.input.checked && e.bedChecking);
  $('submit-note').textContent = files.length && !$('data-policy').checked ? 'Confirm the data policy to continue' : files.length ? files.length + ' independent ' + (files.length === 1 ? 'analysis' : 'analyses') : 'Select a file to begin';
}
function optionsFor(e) {
  if (!e.ref.value || !e.assay.value) throw new Error('Choose a reference genome and library assay for ' + e.file.name + '.');
  const pairs = e.pairs.value.trim() ? e.pairs.value.trim().split(/\n+/).map(line => {
    const columns = line.trim().split(/\s+/);
    if (columns.length !== 3) throw new Error('Index pairs need exactly three columns: name i5 i7.');
    return { name: columns[0], i5: columns[1].toUpperCase(), i7: columns[2].toUpperCase() };
  }) : [];
  if (e.bedCheck.input.checked && (!e.bed.files[0] || !e.bed.files[0].size || e.bed.files[0].size > 4 * 1024 ** 2)) throw new Error('Select a BED file of up to 4 MiB for ' + e.sample.value);
  if (e.bedCheck.input.checked && !e.bedConfirm.checked) throw new Error('Confirm that the BED uses ' + e.bedAssembly + ' coordinates for ' + e.sample.value);
  const qc_reference_ranges = {};
  for (const [key, r] of Object.entries(e.ranges)) qc_reference_ranges[key] = r.disable.checked ? null : Object.fromEntries(Object.entries(r.values).filter(([, i]) => i.value !== '').map(([k, i]) => [k, Number(i.value)]));
  return { sample: e.sample.value, assay: e.assay.value, reference: e.ref.value, mapper: e.mapper.value,
    spike_in: e.spike.input.checked, mycoplasma: e.myco.input.checked, declared_spike_in: e.declared.input.checked,
    custom_regions: e.bedCheck.input.checked, regions_label: e.bedLabel.value, index_pairs: pairs,
    library_structure: e.structure.value(),
    regions_filename: e.bedCheck.input.checked ? e.bed.files[0].name : null,
    regions_mode: e.bedMode.value, regions_assembly: e.bedCheck.input.checked && e.bedConfirm.checked ? e.bedAssembly : null,
    primary_chromosomes: e.chromosomes.value.trim() ? e.chromosomes.value.trim().split(/[,\s]+/) : null,
    params: Object.fromEntries(Object.entries(e.params).filter(([, i]) => i.value !== '').map(([k, i]) => [k, Number(i.value)])), qc_reference_ranges,
    ...Object.fromEntries(Object.entries(e.screen).map(([k, i]) => [k, Number(i.value)])) };
}
function privateLink(task) { return location.origin + '/jobs/' + task.id + '#' + task.token; }
function saveLocation() {
  if (tasks.length === 1) history.replaceState(null, '', '/jobs/' + tasks[0].id + '#' + tasks[0].token);
  else history.replaceState(null, '', '/#jobs=' + tasks.map(t => t.id + ':' + t.token).join(','));
}
function binaryUploadOnce(task, file, endpoint, onprogress) {
  return new Promise((resolve, reject) => {
    const xhr = new XMLHttpRequest(); xhr.open('PUT', '/api/jobs/' + task.id + '/' + endpoint);
    xhr.setRequestHeader('X-Portal-Request', '1'); xhr.setRequestHeader('Authorization', 'Bearer ' + task.token); xhr.setRequestHeader('Content-Type', 'application/octet-stream'); xhr.timeout = 610000;
    xhr.upload.onprogress = e => { if (e.lengthComputable) onprogress(e.loaded / e.total); };
    xhr.onload = () => { let data; try { data = JSON.parse(xhr.responseText); } catch { data = {}; } if (xhr.status >= 200 && xhr.status < 300) resolve(data); else { const error = new Error(data.detail || 'Upload failed.'); error.status = xhr.status; error.retryAfter = Number(xhr.getResponseHeader('Retry-After')) || 5; reject(error); } };
    xhr.onerror = () => reject(new Error('Connection interrupted during upload.')); xhr.ontimeout = () => reject(new Error('Upload timed out.')); xhr.send(file);
  });
}
async function binaryUpload(task, file, endpoint, onprogress) {
  const deadline = Date.now() + 600000;
  while (true) {
    try {
      task.waitingForUpload = false;
      return await binaryUploadOnce(task, file, endpoint, onprogress);
    } catch (error) {
      if (![429, 503].includes(error.status) || Date.now() >= deadline) { task.waitingForUpload = false; throw error; }
      task.waitingForUpload = true; renderTasks();
      await new Promise(resolve => setTimeout(resolve, Math.min(60, Math.max(1, error.retryAfter)) * 1000));
    }
  }
}
const labels = { created: 'Waiting for upload', uploading: 'Uploading', uploaded: 'Upload complete', queued: 'Queued', running: 'Analyzing', succeeded: 'Report ready', failed: 'Analysis failed', cancelling: 'Cancelling', cancelled: 'Cancelled', expired: 'Expired' };
function taskProgress(task) {
  const job = task.job || { status: 'created' };
  const percent = Math.round(Math.min(100, Math.max(0, (task.uploading ? task.uploadProgress * 100 : job.progress) || 0)));
  const deletionPending = job.status !== 'expired' && Boolean(task.deletionRequested || job.deletion_requested);
  const message = deletionPending ? 'Your report is being permanently deleted and is no longer available.' : task.waitingForUpload ? 'Waiting to upload…' : task.uploading ? 'Uploading · ' + percent + '%' : job.status === 'queued' ? (job.queue_ahead ? job.queue_ahead + (job.queue_ahead === 1 ? ' analysis' : ' analyses') + ' ahead of you' : 'Waiting for the analysis worker') : job.status === 'running' ? 'Processing reads and building your report…' : job.status === 'succeeded' ? 'Your interactive report and downloads are ready.' : labels[job.status];
  return { percent, message };
}
function taskView(task) {
  const job = task.job || { status: 'created', sample: task.name || 'Loading task…' };
  const expired = job.status === 'expired';
  const deleted = expired && job.deletion_reason === 'manual';
  const deletionPending = !expired && Boolean(task.deletionRequested || job.deletion_requested);
  const card = node('article', { class: 'task-card' });
  const title = node('div', { class: 'task-title' }, node('h2', { class: !expired && job.sample ? 'task-sample' : '' }, expired ? (deleted ? 'Report deleted' : 'Report expired') : job.sample || 'Analysis'), node('span', { class: 'status ' + (deleted ? 'deleted' : job.status) }, deletionPending ? 'Deletion requested' : deleted ? 'Deleted' : labels[job.status] || job.status));
  card.append(title, node('p', { class: 'task-meta' }, [job.filename ? node('span', { class: 'task-filename' }, job.filename) : null, job.assay === 'microc' ? 'Micro-C' : job.assay === 'chromatin' ? 'Chromatin' : '', job.reference].filter(Boolean).flatMap((part, index) => index ? [' · ', part] : [part])));
  const { percent, message } = taskProgress(task);
  const bar = node('span'); bar.style.width = percent + '%';
  if (expired) {
    const retention = catalog?.retention_days;
    const reason = deleted
      ? 'This report was manually deleted. It can no longer be viewed or downloaded.'
      : 'This report was automatically deleted when its retention period ended.' + (retention ? ' Reports are kept for ' + retention + (retention === 1 ? ' day.' : ' days.') : '');
    card.append(node('div', { class: 'report-unavailable', role: 'status' }, node('p', {}, reason), node('p', {}, 'This cannot be undone. Submit the FASTQ again to generate a new report.')));
  } else {
    if (!deletionPending) card.append(node('div', { class: 'progress', role: 'progressbar', 'aria-label': 'Analysis progress', 'aria-valuemin': '0', 'aria-valuemax': '100', 'aria-valuenow': String(percent) }, bar));
    card.append(node('p', { class: 'progress-label' }, message));
  }
  if (task.failure || task.fetchError || job.error) card.append(node('div', { class: 'error' }, task.failure || task.fetchError || job.error));
  const actions = node('div', { class: 'task-actions' });
  if (job.status === 'succeeded' && !deletionPending) {
    actions.append(node('a', { class: 'primary button', href: job.report_url, target: '_blank', rel: 'noopener' }, 'Open ' + (job.sample || 'sample') + ' report ↗'), node('a', { class: 'secondary button', href: job.download_url }, 'Download report ZIP'));
    card.append(node('p', { class: 'field-note' }, 'Available until ' + new Date(job.expires * 1000).toLocaleString()));
  }
  actions.append(node('button', { type: 'button', class: 'secondary', onclick: async event => { try { await navigator.clipboard.writeText(privateLink(task)); event.target.textContent = 'Link copied'; } catch { saveText('ont-libqc-private-link.txt', privateLink(task)); } } }, 'Copy private link'));
  if (['created', 'uploaded', 'queued', 'running'].includes(job.status) && !task.uploading) actions.append(node('button', { type: 'button', class: 'secondary', onclick: async () => { try { await api('/api/jobs/' + task.id + '/cancel', 'POST', {}, task.token); await refreshTask(task); } catch (e) { task.failure = e.message; renderTasks(); } } }, 'Cancel analysis'));
  if (['succeeded', 'failed', 'cancelled'].includes(job.status)) {
    const warning = node('p', { class: 'deletion-warning', id: 'delete-warning-' + task.id, role: 'status', hidden: true }, 'Permanently deletes this report and its files. This cannot be undone.');
    const deleteButton = node('button', { type: 'button', class: 'secondary delete-report', onclick: async () => {
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
      deleteButton.textContent = task.deleting ? 'Requesting deletion…' : requested ? 'Deletion requested' : confirming ? 'Confirm permanent deletion' : 'Delete report now';
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
  const { progress, updated, ...job } = task.job || {};
  const key = JSON.stringify([job, task.name, task.uploading, task.waitingForUpload, task.failure, task.fetchError, task.confirmDelete, task.deleting, task.deletionRequested, catalog?.retention_days]);
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
  }
  const label = view.card.querySelector('.progress-label');
  if (label && label.textContent !== message) label.textContent = message;
  view.phase = phase;
  return view.card;
}
function renderTasks() {
  const single = tasks.length === 1;
  const ready = tasks.length > 0 && tasks.every(t => t.job?.status === 'succeeded' && !t.deletionRequested && !t.job.deletion_requested);
  const unavailable = tasks.length > 0 && tasks.every(t => t.job?.status === 'expired');
  const needsAttention = tasks.some(t => t.deletionRequested || t.job?.deletion_requested || ['failed', 'cancelled', 'expired'].includes(t.job?.status));
  $('tasks-title').textContent = unavailable
    ? (single ? 'Your report is no longer available.' : 'These reports are no longer available.')
    : ready
    ? (single ? 'Your report is ready.' : 'Your reports are ready.')
    : needsAttention || !tasks.length
      ? (single ? 'Your analysis status.' : 'Your submission status.')
      : (single ? 'Your report is on its way.' : 'Your reports are on their way.');
  $('tasks-eyebrow').textContent = single ? 'YOUR ANALYSIS' : 'YOUR ANALYSES';
  $('tasks-description').textContent = unavailable
    ? 'The reason for removal is shown below. Submit the FASTQ again to generate ' + (single ? 'a new report.' : 'new reports.')
    : single
    ? 'Track your analysis below. Save your private task link to return later.'
    : 'Track each analysis below. Save your private task links to return later.';
  $('save-links').textContent = single ? 'Save private link' : 'Save private links';
  $('submission-view').hidden = true; $('tasks-view').hidden = false;
  syncChildren($('task-list'), tasks.map(updateTaskView));
}
async function refreshTask(task) {
  if (task.uploading) return;
  try {
    task.job = await api('/api/jobs/' + task.id, 'GET', undefined, task.token); task.fetchError = null;
    if (task.job.status === 'succeeded' && !task.job.deletion_requested && (!task.session || Date.now() - task.session > 1200000)) {
      await api('/api/jobs/' + task.id + '/session', 'POST', {}, task.token); task.session = Date.now();
    }
  } catch (e) { task.fetchError = e.message; }
  if (!refreshing) renderTasks();
}
async function refreshTasks() {
  if (refreshing) return;
  refreshing = true;
  try {
    for (const task of tasks) {
      if (Date.now() < retryRequestsAt) break;
      if (!task.job || !['succeeded', 'expired', 'cancelled', 'failed'].includes(task.job.status) || task.fetchError || (task.job?.status !== 'expired' && (task.deletionRequested || task.job?.deletion_requested || task.job?.expires <= Date.now() / 1000)) || (task.job.status === 'succeeded' && (!task.session || Date.now() - task.session > 1200000))) await refreshTask(task);
    }
  } finally { refreshing = false; renderTasks(); }
}
function saveText(name, text) { const url = URL.createObjectURL(new Blob([text], { type: 'text/plain' })); const a = node('a', { href: url, download: name }); document.body.append(a); a.click(); a.remove(); setTimeout(() => URL.revokeObjectURL(url), 1000); }
$('save-links').addEventListener('click', () => saveText(tasks.length === 1 ? 'ont-libqc-private-link.txt' : 'ont-libqc-private-links.txt', tasks.map(t => (t.job?.sample || t.name || t.id) + '\n' + privateLink(t)).join('\n\n')));
$('data-policy').checked = false;
$('data-policy').addEventListener('change', updateCount);
window.addEventListener('pageshow', event => {
  if (event.persisted) { $('data-policy').checked = false; updateCount(); }
});

$('submission-form').addEventListener('submit', async event => {
  event.preventDefault(); if (submitting || !files.length) return;
  showError(''); let entries;
  try {
    if (!$('data-policy').checked) {
      $('data-policy').focus();
      throw new Error('Confirm that this submission contains no human patient data or other sensitive personal data.');
    }
    if (files.length > catalog.max_files) throw new Error('Choose up to ' + catalog.max_files + ' files.');
    entries = files.map(e => ({ entry: e, options: optionsFor(e), bedFile: e.bed.files[0] }));
  } catch (e) { showError(e.message); return; }
  submitting = true; updateCount();
  try {
    for (const { entry, options, bedFile } of entries) {
      if (options.custom_regions && await checkedBedFile(entry.bed, 4 * 1024 ** 2) !== bedFile) throw new Error('The BED selection changed. Review the file and submit again.');
    }
    const batch = await api('/api/batches', 'POST', { files: entries.map(({ entry, options }) => ({
      email: $('email').value, filename: entry.file.name, size: entry.file.size, options,
      consent: $('consent').checked, data_policy_acknowledged: $('data-policy').checked, submission_key: $('submission-code').value, website: $('website').value
    })) });
    tasks = batch.jobs.map((job, i) => ({ id: job.id, token: job.token, name: entries[i].options.sample,
      uploadProgress: 0, job: { ...job, ...entries[i].options, filename: entries[i].entry.file.name } }));
    saveLocation(); renderTasks();
    for (let i = 0; i < entries.length; i++) {
      const { entry, options, bedFile } = entries[i], task = tasks[i];
      if (['cancelled', 'cancelling', 'expired'].includes(task.job.status)) continue;
      try {
        task.uploading = true; renderTasks();
        await binaryUpload(task, entry.file, 'upload', progress => { task.uploadProgress = progress; renderTasks(); });
        if (options.custom_regions) await binaryUpload(task, bedFile, 'regions', () => {});
        task.job = await api('/api/jobs/' + task.id + '/submit', 'POST', {}, task.token);
        task.uploading = false; renderTasks();
      } catch (e) {
        task.failure = e.message + ' Cancel this task and submit the file again.';
        task.uploading = false; await refreshTask(task);
      }
    }
  } catch (e) { showError(e.message); }
  finally { submitting = false; updateCount(); }
});
async function init() {
  try {
    catalog = await api('/api/catalog');
    $('availability').textContent = catalog.worker_online ? 'Ready for analysis' : 'Worker offline'; $('availability').classList.toggle('online', catalog.worker_online);
    $('upload-limit').textContent = 'FASTQ / FASTQ.gz · ' + bytes(catalog.max_upload_bytes) + ' per file · up to ' + catalog.max_files + ' files';
    $('pipeline-link').href = catalog.pipeline_url;
    $('large-file-note').textContent = 'Larger than ' + bytes(catalog.max_upload_bytes) + '? Run the pipeline locally. Gzip files must also expand to at most ' + bytes(catalog.max_expanded_bytes) + '. Online analyses are limited to ' + catalog.max_reads.toLocaleString() + ' reads per file.';
    $('retention').textContent = catalog.retention_days + ' days';
    $('email-note').textContent = catalog.mail_enabled ? 'We’ll send one summary email with all your results when every file in this submission has finished.' : 'Email delivery is not enabled on this server. Save your private task links to collect the reports here.';
    $('privacy-delivery-note').textContent = catalog.mail_enabled ? 'Your private link will also arrive by email. Keep it safe and share it only with people you trust.' : 'Save your task links before leaving. Without email delivery, lost links cannot be recovered here.';
    $('code-field').hidden = !catalog.submission_key_required; $('submission-code').required = catalog.submission_key_required;
    const single = location.pathname.match(/^\/jobs\/([0-9a-f]{32})$/);
    if (single && /^[0-9a-f]{64}$/.test(location.hash.slice(1))) tasks = [{ id: single[1], token: location.hash.slice(1) }];
    else if (location.hash.startsWith('#jobs=')) tasks = location.hash.slice(6).split(',').slice(0, 64).map(t => t.split(':')).filter(([id, token]) => /^[0-9a-f]{32}$/.test(id) && /^[0-9a-f]{64}$/.test(token)).map(([id, token]) => ({ id, token }));
    if (single && !tasks.length) { $('submission-view').hidden = true; showError(new URLSearchParams(location.search).has('report') ? 'Open the complete private task link from your email or saved links to access this report. The report URL alone does not include its access key.' : 'This task needs its complete private link, including the part after #.', 'global-error'); }
    if (tasks.length) { renderTasks(); await refreshTasks(); }
    updateCount();
  } catch (e) { showError('Cannot connect to the analysis service. ' + e.message, 'global-error'); $('submit-button').disabled = true; }
}
$('fastq-input').addEventListener('change', e => { if (catalog) addFiles(e.target.files); e.target.value = ''; });
$('dropzone').addEventListener('keydown', e => { if (e.key === 'Enter' || e.key === ' ') { e.preventDefault(); $('fastq-input').click(); } });
for (const type of ['dragenter', 'dragover']) $('dropzone').addEventListener(type, e => { e.preventDefault(); $('dropzone').classList.add('drag'); });
for (const type of ['dragleave', 'drop']) $('dropzone').addEventListener(type, e => { e.preventDefault(); $('dropzone').classList.remove('drag'); });
$('dropzone').addEventListener('drop', e => { if (catalog && !submitting) addFiles(e.dataTransfer.files); });
async function pollTasks() {
  if (tasks.length) await refreshTasks();
  setTimeout(pollTasks, Math.max(document.hidden ? 30000 : 10000, retryRequestsAt - Date.now()) + Math.random() * 2000);
}
setTimeout(pollTasks, 10000);
init();
