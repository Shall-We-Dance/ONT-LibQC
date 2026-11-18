(() => {
  const labels = new Map();
  const scopes = new Map();
  for (const cell of document.querySelectorAll('.pair-table td.scope-name')) {
    const id = cell.dataset.scopeId;
    const name = cell.textContent;
    if (name === '—' || name === id || name === cell.dataset.original) continue;
    labels.set(id, name);
    scopes.set(id, id);
    scopes.set(cell.dataset.original, id);
  }
  function addLabel(element, id) {
    if (!labels.has(id)) return;
    const label = document.createElement('span');
    label.className = 'scope-display-name';
    label.textContent = ' · ' + labels.get(id);
    const count = element.querySelector('small, .mixture-count, strong');
    element.insertBefore(label, count && count.parentNode === element ? count : null);
  }
  for (const button of document.querySelectorAll('.profile-comparison-legend button[data-scope]')) {
    addLabel(button, button.dataset.scope);
    if (labels.has(button.dataset.scope)) button.setAttribute('aria-label', 'Highlight ' + button.textContent.trim());
  }
  for (const summary of document.querySelectorAll('details.profile-detail[data-scope] > summary')) {
    addLabel(summary, summary.parentElement.dataset.scope);
  }
  for (const el of document.querySelectorAll('.table-scroll th, .rna-composition-card thead th, .rna-profile-panels figcaption, .mixture-label')) {
    const text = Array.from(el.childNodes).filter(n => n.nodeType === Node.TEXT_NODE).map(n => n.textContent).join('').trim();
    addLabel(el, scopes.get(text));
  }
})();
