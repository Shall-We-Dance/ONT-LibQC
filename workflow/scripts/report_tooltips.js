
(() => {
  const tip = document.createElement('div');
  tip.id = 'metric-tooltip'; tip.setAttribute('role', 'tooltip'); tip.hidden = true;
  document.body.appendChild(tip);
  let active = null;
  function hide() {
    if (active) active.removeAttribute('aria-describedby');
    active = null; tip.hidden = true;
  }
  function show(target) {
    hide(); active = target;
    tip.textContent = target.dataset.tip;
    tip.hidden = false; tip.style.left = '8px'; tip.style.top = '8px';
    target.setAttribute('aria-describedby', tip.id);
    const r = target.getBoundingClientRect(), t = tip.getBoundingClientRect();
    const w = document.documentElement.clientWidth, h = window.innerHeight;
    const left = Math.max(8, Math.min(r.left, w - t.width - 8));
    const below = r.bottom + 8;
    const top = below + t.height <= h - 8 ? below : Math.max(8, r.top - t.height - 8);
    tip.style.left = left + 'px'; tip.style.top = Math.min(top, Math.max(8, h - t.height - 8)) + 'px';
  }
  document.addEventListener('mouseover', e => { const t = e.target.closest('.metric-tip'); if (t && t !== active) show(t); });
  document.addEventListener('mouseout', e => { if (active && active.contains(e.target) && !active.contains(e.relatedTarget)) hide(); });
  document.addEventListener('focusin', e => { if (e.target.matches('.metric-tip')) show(e.target); });
  document.addEventListener('focusout', e => { if (e.target.matches('.metric-tip')) hide(); });
  document.addEventListener('click', e => { const t = e.target.closest('.metric-tip'); if (t) show(t); else hide(); });
  document.addEventListener('keydown', e => { if (e.key === 'Escape') hide(); });
  window.addEventListener('scroll', hide, true); window.addEventListener('resize', hide);
})();
