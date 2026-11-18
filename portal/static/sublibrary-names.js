(() => {
  const key = location.pathname.split('/').pop();
  const report = '/r/' + key + '#library-shape';
  const api = '/api/sublibrary-names/' + key;
  const back = document.getElementById('back');
  const status = document.getElementById('editor-status');
  const form = document.getElementById('names-form');
  back.href = report;
  let token;
  try { token = sessionStorage.getItem('ont-sublibrary-editor'); } catch {}
  if (!token) {
    status.textContent = 'To edit names, open this report from your private analysis page.';
    return;
  }
  const headers = {Authorization: 'Bearer ' + token, 'X-Portal-Request': '1'};
  let revision;
  const inputs = new Map();
  async function load() {
    const response = await fetch(api, {headers});
    const data = await response.json();
    if (!response.ok) throw new Error(data.detail);
    revision = data.revision;
    status.textContent = data.sample;
    for (const [id, scope] of Object.entries(data.scopes)) {
      const row = form.querySelector('tbody').insertRow();
      const library = document.createElement('th');
      library.scope = 'row';
      library.textContent = id;
      row.append(library);
      const input = document.createElement('input');
      input.type = 'text';
      input.setAttribute('aria-label', 'Name for ' + id);
      input.value = data.names[id] || '';
      input.placeholder = scope.scope_name || id;
      row.insertCell().append(input);
      for (const index of [scope.i5, scope.i7]) {
        row.insertCell().textContent = index === '' ? 'No index' : (index ?? 'Not detected');
      }
      inputs.set(id, input);
    }
    form.hidden = false;
  }
  load().catch(error => { status.textContent = error.message; });
  form.addEventListener('submit', async event => {
    event.preventDefault();
    const button = form.querySelector('button');
    button.disabled = true;
    try {
      const names = Object.fromEntries(Array.from(inputs, ([id, input]) => [id, input.value]));
      const response = await fetch(api, {method: 'PUT', headers: {...headers, 'Content-Type': 'application/json'},
        body: JSON.stringify({revision, names})});
      if (!response.ok) throw new Error((await response.json()).detail);
      location.replace(report);
    } catch (error) { status.textContent = error.message; button.disabled = false; }
  });
})();
