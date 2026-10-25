'use strict';

(() => {
  const button = document.getElementById('copy-citation');
  const citation = document.getElementById('citation-text');
  const status = document.getElementById('citation-status');
  if (!button || !citation || !status) return;
  button.hidden = false;
  button.addEventListener('click', async () => {
    button.disabled = true;
    status.textContent = '';
    try {
      await navigator.clipboard.writeText(citation.textContent.replace(/\s+/g, ' ').trim());
      status.textContent = 'Citation copied.';
    } catch {
      citation.focus();
      const selection = window.getSelection();
      const range = document.createRange();
      range.selectNodeContents(citation);
      selection.removeAllRanges();
      selection.addRange(range);
      status.textContent = 'Select and copy the citation above.';
    } finally {
      button.disabled = false;
    }
  });
})();
