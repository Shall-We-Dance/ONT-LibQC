(() => {
  const dialog = document.getElementById('sc-name-editor');
  const table = document.querySelector('#library-quality .pair-table');
  if (!dialog || !table) return;

  const inputs = new Map(Array.from(dialog.querySelectorAll('input[data-scope-id]'), input => [input.dataset.scopeId, input]));
  const cells = new Map();
  for (const row of table.querySelectorAll('tr[data-sublibrary]')) {
    const id = row.dataset.sublibrary;
    if (!inputs.has(id)) continue;
    const cell = row.insertCell(1);
    cell.className = 'scope-name';
    cell.dataset.scopeId = id;
    cell.textContent = '—';
    cells.set(id, cell);
  }
  if (!cells.size) return;
  const header = document.createElement('th');
  header.scope = 'col';
  header.textContent = 'Name';
  table.rows[0].insertBefore(header, table.rows[0].cells[1]);
  const form = dialog.querySelector('form');
  const fields = form.querySelector('fieldset');
  const save = form.querySelector('[type=submit]');
  const status = form.querySelector('.sc-name-status');
  const endpoint = document.querySelector('meta[name="sublibrary-names-api"]')?.content;
  const key = 'ont-libqc:scrna:names:' + dialog.dataset.reportKey;
  dialog.querySelector('.sc-name-storage').textContent = endpoint ?
    'Names are shared with this report. Leave a name blank to clear it.' :
    'Names are saved in this browser. Leave a name blank to clear it.';

  // Match scope labels only where the report renders a physical library identifier.
  const labels = [];
  for (const el of document.querySelectorAll('main th, main .mixture-label, #chemistry h4, main details > summary, main [data-scope-label]')) {
    const text = Array.from(el.childNodes).filter(n => n.nodeType === Node.TEXT_NODE).map(n => n.textContent).join('').trim();
    const id = el.dataset.scopeLabel || text.replace(/^Physical library · /, '');
    if (!cells.has(id)) continue;
    const label = document.createElement('span');
    label.className = 'scope-display-name';
    el.insertBefore(label, el.querySelector(':scope > strong, :scope > small'));
    labels.push([id, label]);
  }

  function checkedName(id, value) {
    if (typeof value !== 'string') throw new Error('Invalid name for ' + id + '.');
    const name = value.normalize('NFC').replace(/^ +| +$/g, '');
    if (name && ([...name].length > 48 || !/^[\p{L}\p{N}][\p{L}\p{N}\p{M} ._()+-]*$/u.test(name))) {
      throw new Error(id + ': use up to 48 characters, starting with a letter or number. Spaces and . _ ( ) + - are allowed.');
    }
    return name;
  }
  function display(data) {
    for (const [id, cell] of cells) cell.textContent = data.names[id] || '—';
    for (const [id, label] of labels) label.textContent = data.names[id] ? ' · ' + data.names[id] + ' ' : '';
  }
  function localData() {
    return JSON.parse(localStorage.getItem(key) || '{"revision":0,"names":{}}');
  }
  async function readOrSave(update) {
    let data;
    if (endpoint) {
      const url = new URL(endpoint, location.href);
      if (url.origin !== location.origin || !['http:', 'https:'].includes(url.protocol)) throw new Error('Invalid report address.');
      const options = {credentials: 'same-origin', cache: 'no-store', redirect: 'error'};
      if (update) Object.assign(options, {method: 'PUT', headers: {'Content-Type': 'application/json', 'X-Portal-Request': '1'},
        body: JSON.stringify(update)});
      const response = await fetch(url, options);
      if (!response.ok) {
        if (response.status === 409) throw new Error('Names changed in another tab. Cancel and reopen to see them.');
        throw new Error('Could not ' + (update ? 'save' : 'load') + ' names. Please try again.');
      }
      data = await response.json();
    } else {
      data = localData();
      if (update) {
        if (data.revision !== update.revision) throw new Error('Names changed in another tab. Cancel and reopen to see them.');
        data = {revision: update.revision + 1, names: update.names};
        localStorage.setItem(key, JSON.stringify(data));
      }
    }
    if (!Number.isSafeInteger(data?.revision) || data.revision < 0 || !data.names || typeof data.names !== 'object' ||
        Array.isArray(data.names)) throw new Error('Saved names could not be read.');
    const names = Object.create(null);
    for (const id of cells.keys()) {
      if (Object.hasOwn(data.names, id)) names[id] = checkedName(id, data.names[id]);
    }
    return {revision: data.revision, names};
  }

  let revision;
  let request = 0;
  let saving = false;
  async function load(editor) {
    const current = ++request;
    try {
      const data = await readOrSave();
      if (current !== request) return;
      display(data);
      if (editor && dialog.open) {
        revision = data.revision;
        for (const [id, input] of inputs) input.value = data.names[id] || '';
        fields.disabled = false;
        save.disabled = false;
        inputs.values().next().value.focus();
      }
    } catch (error) {
      if (editor && current === request) status.textContent = error.message;
    }
  }
  let opener;
  function openEditor(event) {
    opener = event.currentTarget;
    form.reset();
    fields.disabled = true;
    save.disabled = true;
    status.textContent = '';
    dialog.showModal();
    load(true);
  }
  for (const heading of [table.closest('.identity-card').querySelector('h3'),
                         document.querySelector('#sublibrary-composition > h2')]) {
    if (!heading) continue;
    const edit = document.createElement('button');
    edit.type = 'button';
    edit.className = 'sc-name-edit';
    edit.textContent = 'Edit names';
    edit.setAttribute('aria-haspopup', 'dialog');
    edit.setAttribute('aria-controls', dialog.id);
    edit.addEventListener('click', openEditor);
    heading.append(edit);
  }
  form.querySelector('[type=button]').addEventListener('click', () => dialog.close());
  dialog.addEventListener('cancel', event => { if (saving) event.preventDefault(); });
  dialog.addEventListener('close', () => opener?.focus());
  form.addEventListener('submit', async event => {
    event.preventDefault();
    if (save.disabled) return;
    save.disabled = true;
    saving = true;
    fields.disabled = true;
    const cancel = form.querySelector('[type=button]');
    cancel.disabled = true;
    try {
      const names = Object.create(null);
      for (const [id, input] of inputs) {
        const name = checkedName(id, input.value);
        if (name) names[id] = name;
      }
      const data = await readOrSave({revision, names});
      display(data);
      dialog.close();
    } catch (error) {
      status.textContent = error.message;
    } finally {
      cancel.disabled = false;
      fields.disabled = false;
      save.disabled = false;
      saving = false;
    }
  });
  window.addEventListener('storage', event => {
    if (!endpoint && event.key === key && !dialog.open) load(false);
  });
  load(false);
})();
