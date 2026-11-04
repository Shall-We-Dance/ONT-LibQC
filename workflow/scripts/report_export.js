/* Native PDF export: keep vector figures, browser-only BED plots, and offline use. */
(() => {
  'use strict';
  const button = document.getElementById('report-export-pdf');
  const status = document.getElementById('report-export-status');
  if (!button || !status) return;
  let snapshot = null;
  let preparing = false;
  let printing = false;
  let restoreTimer = null;

  function glossary() {
    const entries = new Map();
    document.querySelectorAll('[data-tip]').forEach(node => {
      const description = node.getAttribute('data-tip').trim();
      if (!description || entries.has(description)) return;
      const row = node.closest('tr');
      const heading = node.closest('h2,h3,h4,summary');
      const labelNode = row?.cells[0] || heading || node;
      const label = labelNode.textContent.replace(/\s+/g, ' ').replace(/\s*[?ⓘ]+\s*$/, '').trim();
      if (label && label !== '?') entries.set(description, label);
    });
    if (!entries.size) return;
    const section = document.createElement('section');
    section.className = 'report-print-definitions';
    section.dataset.exportTemporary = 'true';
    const title = document.createElement('h2');
    title.textContent = 'Metric notes';
    section.append(title);
    const list = document.createElement('dl');
    entries.forEach((label, description) => {
      const item = document.createElement('div');
      const term = document.createElement('dt');
      term.textContent = label;
      const meaning = document.createElement('dd');
      meaning.textContent = description;
      item.append(term, meaning);
      list.append(item);
    });
    section.append(list);
    document.body.append(section);
  }

  function prepare() {
    if (snapshot) return;
    snapshot = {
      title: document.title,
      details: Array.from(document.querySelectorAll('details'), node => ({node, open: node.open, name: node.getAttribute('name')})),
      images: Array.from(document.images, node => ({node, loading: node.getAttribute('loading')})),
      links: [],
      headers: [],
      scroll: [window.scrollX, window.scrollY]
    };
    snapshot.details.forEach(({node}) => { node.removeAttribute('name'); node.open = true; });
    snapshot.images.forEach(({node}) => { node.loading = 'eager'; });
    document.querySelectorAll('table').forEach(table => {
      const row = table.rows[0];
      if (table.tHead || !row || !row.cells.length || !Array.from(row.cells).every(cell => cell.tagName === 'TH')) return;
      const parent = row.parentNode;
      const next = row.nextSibling;
      const head = table.createTHead();
      head.append(row);
      snapshot.headers.push({row, parent, next, head});
    });
    // Relative links, including fragments, inherit the private page's access grant.
    // Remove them while printing so PDF annotations cannot expose that grant.
    document.querySelectorAll('a[href]').forEach(node => {
      const value = node.getAttribute('href');
      if (!/^https?:\/\//i.test(value)) {
        snapshot.links.push({node, value});
        node.removeAttribute('href');
      } else if (/[?&]view=/.test(value)) {
        snapshot.links.push({node, value});
        node.removeAttribute('href');
      }
    });
    const sample = document.querySelector('.report-hero h1')?.textContent.trim();
    if (sample) document.title = sample + ' - ONT-LibQC';
    glossary();
    document.documentElement.classList.add('report-exporting');
  }

  function restore() {
    clearTimeout(restoreTimer);
    printing = false;
    if (snapshot) {
      const saved = snapshot;
      snapshot = null;
      document.title = saved.title;
      saved.details.forEach(({node, open, name}) => { node.open = open; if (name !== null) node.setAttribute('name', name); });
      saved.images.forEach(({node, loading}) => { if (loading === null) node.removeAttribute('loading'); else node.setAttribute('loading', loading); });
      saved.links.forEach(({node, value}) => node.setAttribute('href', value));
      saved.headers.forEach(({row, parent, next, head}) => {
        parent.insertBefore(row, next?.parentNode === parent ? next : null);
        head.remove();
      });
      document.querySelectorAll('[data-export-temporary]').forEach(node => node.remove());
      document.documentElement.classList.remove('report-exporting');
      window.scrollTo({left: saved.scroll[0], top: saved.scroll[1], behavior: 'instant'});
    }
    button.disabled = false;
    button.removeAttribute('aria-busy');
    button.querySelector('span').textContent = 'Export PDF';
  }

  async function readyImages() {
    const images = Array.from(document.images);
    await Promise.all(images.map(image => {
      if (image.complete) return Promise.resolve();
      return new Promise(resolve => {
        const finish = () => { image.removeEventListener('load', finish); image.removeEventListener('error', finish); resolve(); };
        image.addEventListener('load', finish, {once: true});
        image.addEventListener('error', finish, {once: true});
      });
    }));
    if (images.some(image => !image.naturalWidth)) throw new Error('figures');
    await Promise.all(images.map(image => image.decode ? image.decode().catch(() => {}) : Promise.resolve()));
    if (document.fonts?.ready) await document.fonts.ready;
  }

  async function exportPDF() {
    if (preparing || printing) return;
    preparing = true;
    button.disabled = true;
    button.setAttribute('aria-busy', 'true');
    button.querySelector('span').textContent = 'Preparing…';
    status.textContent = 'Preparing all sections and figures…';
    let timeout;
    try {
      prepare();
      await Promise.race([readyImages(), new Promise((_, reject) => { timeout = setTimeout(() => reject(new Error('timeout')), 20000); })]);
      await new Promise(resolve => requestAnimationFrame(() => requestAnimationFrame(resolve)));
      status.textContent = 'Choose Save as PDF. Turn off headers and footers for a clean report without the private page URL.';
      printing = true;
      button.querySelector('span').textContent = 'Save as PDF…';
      window.print();
    } catch (error) {
      restore();
      status.textContent = 'Some figures could not be prepared. Reload the report and try Export PDF again.';
    } finally {
      clearTimeout(timeout);
      preparing = false;
    }
  }

  button.addEventListener('click', exportPDF);
  window.addEventListener('beforeprint', prepare);
  window.addEventListener('afterprint', restore);
  const media = window.matchMedia('print');
  const mediaChange = event => { if (event.matches) prepare(); else restore(); };
  if (media.addEventListener) media.addEventListener('change', mediaChange);
  else media.addListener(mediaChange);
  // Some browsers omit afterprint on cancellation; restore when the page regains focus.
  window.addEventListener('focus', () => {
    if (printing && !preparing && !media.matches) restoreTimer = setTimeout(restore, 300);
  });
  // Also used by automated print-layout checks without opening a native dialog.
  window.ONTReportExport = Object.freeze({prepare, restore});
})();
