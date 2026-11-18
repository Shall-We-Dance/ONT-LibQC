(() => {
  const key = location.pathname.split('/').pop();
  const report = '/r/' + key + '#library-shape';
  const api = '/api/sublibrary-names/' + key;
  const back = document.getElementById('back');
  const status = document.getElementById('editor-status');
  const form = document.getElementById('names-form');
  const unlock = document.getElementById('unlock-form');
  const linkInput = document.getElementById('analysis-link');
  back.href = report;
  let token;
  try { token = sessionStorage.getItem('ont-sublibrary-editor'); } catch {}
  let revision;
  const inputs = new Map();

  function showUnlock() {
    status.textContent = 'Paste the link from Copy link on your analysis page.';
    unlock.hidden = false;
  }
  function privateAnalysisPath(value) {
    let url;
    try { url = new URL(value.trim(), location.origin); } catch {}
    if (url?.origin === location.origin && !url.username && !url.password && url.hash &&
        (['/j/', '/b/', '/'].includes(url.pathname) || url.pathname.startsWith('/jobs/'))) {
      return url.pathname + url.hash;
    }
    throw new Error('Use the private link from Copy link on your analysis page.');
  }
  function showNames(data) {
    unlock.hidden = true;
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

  async function load() {
    if (!token) return showUnlock();
    const response = await fetch(api, {headers: {Authorization: 'Bearer ' + token}});
    if (response.status === 403) return showUnlock();
    const data = await response.json();
    if (!response.ok) throw new Error(data.detail);
    showNames(data);
  }
  unlock.addEventListener('submit', async event => {
    event.preventDefault();
    const button = unlock.querySelector('button');
    button.disabled = true;
    try {
      const link = privateAnalysisPath(linkInput.value);
      const response = await fetch(api + '/unlock', {method: 'POST', redirect: 'error',
        headers: {'X-Portal-Request': '1', 'Content-Type': 'application/json'},
        body: JSON.stringify({link})});
      const data = await response.json();
      if (!response.ok) {
        if (response.status === 404 || response.status === 410) unlock.hidden = true;
        throw new Error(data.detail);
      }
      token = data.name_edit_token;
      linkInput.value = '';
      try { sessionStorage.setItem('ont-sublibrary-editor', token); } catch {}
      showNames(data);
      inputs.values().next().value?.focus();
    } catch (error) { status.textContent = error.message; }
    finally { button.disabled = false; }
  });
  load().catch(error => { status.textContent = error.message; });
  form.addEventListener('submit', async event => {
    event.preventDefault();
    const button = form.querySelector('button');
    button.disabled = true;
    try {
      const names = Object.fromEntries(Array.from(inputs, ([id, input]) => [id, input.value]));
      const response = await fetch(api, {method: 'PUT',
        headers: {Authorization: 'Bearer ' + token, 'X-Portal-Request': '1', 'Content-Type': 'application/json'},
        body: JSON.stringify({revision, names})});
      if (!response.ok) throw new Error((await response.json()).detail);
      location.replace(report);
    } catch (error) { status.textContent = error.message; button.disabled = false; }
  });
})();
