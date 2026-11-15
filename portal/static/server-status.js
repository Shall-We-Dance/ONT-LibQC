/* Public aggregate counts; independent of private task state. */
window.ONTServerStatus = function ({ $, node }) {
function drawJobHistory(history, total) {
  const container = $('jobs-history-chart');
  if (!history.length) { container.replaceChildren(node('p', {}, 'No jobs yet.')); return; }
  const day = 86400000;
  const values = history.map(point => ({time: Date.parse(point.date + 'T00:00:00Z'), total: point.total}));
  const start = values[0].time - day;
  const end = Math.max(values.at(-1).time, Math.floor(Date.now() / day) * day);
  const left = 34, right = 306, top = 12, bottom = 136;
  const max = Math.max(1, total);
  const x = time => left + (time - start) / (end - start) * (right - left);
  const y = count => bottom - count / max * (bottom - top);
  const svgNode = (tag, attrs, text) => {
    const el = document.createElementNS('http://www.w3.org/2000/svg', tag);
    for (const [key, value] of Object.entries(attrs)) el.setAttribute(key, value);
    if (text !== undefined) el.textContent = text;
    return el;
  };
  const svg = svgNode('svg', {viewBox: '0 0 324 166', role: 'img', 'aria-label': 'Cumulative jobs by date (UTC), ' + total.toLocaleString() + ' total jobs'});
  const dateLabel = time => new Date(time).toLocaleDateString(undefined, {month: 'short', day: 'numeric', year: '2-digit', timeZone: 'UTC'});
  svg.append(svgNode('title', {}, 'Cumulative jobs · ' + dateLabel(start) + ' – ' + dateLabel(end)));
  for (const count of [...new Set([0, Math.round(max / 2), max])]) {
    svg.append(svgNode('line', {x1: left, y1: y(count), x2: right, y2: y(count), class: 'history-grid'}));
    svg.append(svgNode('text', {x: left - 6, y: y(count) + 4, 'text-anchor': 'end'}, count.toLocaleString()));
  }
  // Steps preserve zero-growth intervals instead of inventing intermediate jobs.
  let path = `M ${left} ${bottom}`;
  for (const point of values) path += ` H ${x(point.time)} V ${y(point.total)}`;
  path += ` H ${right}`;
  svg.append(svgNode('path', {d: path + ` V ${bottom} H ${left} Z`, class: 'history-area'}));
  svg.append(svgNode('path', {d: path, class: 'history-line'}));
  svg.append(svgNode('text', {x: left, y: 158}, dateLabel(start)));
  svg.append(svgNode('text', {x: right, y: 158, 'text-anchor': 'end'}, dateLabel(end)));
  container.replaceChildren(svg);
}
const serverStatus = $('server-status');
// Keep the native summary keyboard behavior while making the popover easy to dismiss.
document.addEventListener('click', event => {
  if (serverStatus.open && !serverStatus.contains(event.target)) serverStatus.open = false;
});
document.addEventListener('keydown', event => {
  if (event.key === 'Escape' && serverStatus.open) {
    serverStatus.open = false;
    serverStatus.querySelector('summary').focus();
  }
});
async function refreshServerStatus() {
  const keys = ['running', 'queued', 'total'];
  try {
    const response = await fetch('/api/statistics', {credentials: 'omit', cache: 'no-store'});
    if (!response.ok) throw new Error('Status unavailable');
    const data = await response.json();
    if (keys.some(key => !Number.isSafeInteger(data[key]) || data[key] < 0)) throw new Error('Invalid counts');
    if (!Array.isArray(data.history) || data.history.some((point, i) =>
      !/^\d{4}-\d{2}-\d{2}$/.test(point.date) || !Number.isFinite(Date.parse(point.date)) ||
      !Number.isSafeInteger(point.total) || point.total <= 0 ||
      (i && (point.date <= data.history[i - 1].date || point.total <= data.history[i - 1].total))) ||
      (data.history.at(-1)?.total || 0) !== data.total) throw new Error('Invalid history');
    for (const key of keys) $('jobs-' + key).textContent = data[key].toLocaleString();
    drawJobHistory(data.history, data.total);
  } catch {
    for (const key of keys) $('jobs-' + key).textContent = '—';
    $('jobs-history-chart').replaceChildren(node('p', {}, 'History is temporarily unavailable.'));
  }
  $('jobs-queued-mobile').textContent = $('jobs-queued').textContent;
  setTimeout(refreshServerStatus, document.hidden ? 60000 : 15000);
}
refreshServerStatus();

};
