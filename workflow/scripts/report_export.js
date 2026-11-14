/* Native PDF export: keep vector figures, browser-only BED plots, and offline use. */
(() => {
  'use strict';
  const button = document.getElementById('report-export-pdf');
  const status = document.getElementById('report-export-status');
  const contents = document.getElementById('report-export-content');
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
      if (node.matches('.metric-status') && /^(Descriptive only|Within range|Borderline|Outside range)$/.test(description)) return;
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

  // Split wide QC tables across groups of at most five scopes. Each group
  // repeats the metric label and reference, keeping every number readable.
  function printTables() {
    document.querySelectorAll('table.qc-metrics').forEach(table => {
      if (table.closest('.report-print-omit')) return;
      const columns = table.rows[0]?.cells.length || 0;
      if (columns <= 7) return;
      const wrapper = document.createElement('div');
      wrapper.className = 'report-print-tables';
      wrapper.dataset.exportTemporary = 'true';
      for (let start = 1; start < columns - 1; start += 5) {
        const indices = [0];
        for (let i = start; i < Math.min(start + 5, columns - 1); i++) indices.push(i);
        indices.push(columns - 1);
        const copy = table.cloneNode(true);
        copy.removeAttribute('id');
        copy.querySelectorAll('[id]').forEach(node => node.removeAttribute('id'));
        Array.from(copy.rows).forEach(row => {
          Array.from(row.cells).forEach((cell, index) => { if (!indices.includes(index)) cell.remove(); });
        });
        const caption = copy.createCaption();
        caption.textContent = indices.slice(1, -1).map(index => table.rows[0].cells[index].textContent.trim()).join(' · ');
        wrapper.append(copy);
      }
      table.after(wrapper);
      table.classList.add('report-print-source');
      snapshot.printSources.push(table);
    });
  }

  function omit(node) {
    node.classList.add('report-print-omit');
    snapshot.omitted.push(node);
  }

  function prepare() {
    if (snapshot) return;
    snapshot = {
      title: document.title,
      details: Array.from(document.querySelectorAll('details'), node => ({node, open: node.open, name: node.getAttribute('name')})),
      images: Array.from(document.images, node => ({node, loading: node.getAttribute('loading')})),
      links: [],
      headers: [],
      printSources: [],
      omitted: [],
      scroll: [window.scrollX, window.scrollY]
    };
    const full = contents?.value === 'full';
    document.documentElement.classList.add(full ? 'report-print-full' : 'report-print-overview');
    snapshot.details.forEach(({node, open}) => {
      if (full) { node.removeAttribute('name'); node.open = true; }
      else if (!open) omit(node);
    });
    const custom = document.getElementById('custom-regions');
    if (custom && !custom.querySelector('#custom-output svg')) omit(custom);
    if (!full) document.querySelectorAll('.library-sublibraries').forEach(omit);
    snapshot.images.forEach(({node}) => { if (!node.closest('.report-print-omit')) node.loading = 'eager'; });
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
    if (full) glossary();
    printTables();
    const note = document.createElement('p');
    note.className = 'report-print-note';
    note.dataset.exportTemporary = 'true';
    note.textContent = full ? 'Complete report · includes detailed sections and metric notes' : 'Overview · visible sections; full details are available in the HTML report';
    document.querySelector('.report-hero')?.append(note);
    if (document.querySelector('.metric-unavailable')) {
      const legend = document.createElement('p');
      legend.className = 'report-print-note';
      legend.dataset.exportTemporary = 'true';
      legend.textContent = '— Not attributable to an individual sublibrary. Coloured markers indicate evaluated QC results.';
      document.body.append(legend);
    }
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
      saved.printSources.forEach(node => node.classList.remove('report-print-source'));
      saved.omitted.forEach(node => node.classList.remove('report-print-omit'));
      document.documentElement.classList.remove('report-print-full', 'report-print-overview');
      document.documentElement.classList.remove('report-exporting');
      window.scrollTo({left: saved.scroll[0], top: saved.scroll[1], behavior: 'instant'});
    }
    button.disabled = false;
    button.removeAttribute('aria-busy');
    button.querySelector('span').textContent = 'Export PDF';
  }

  async function readyImages() {
    const images = Array.from(document.images).filter(image => !image.closest('.report-print-omit'));
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
    status.textContent = contents?.value === 'full' ? 'Preparing the complete report…' : 'Preparing the report overview…';
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
