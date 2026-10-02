"""Static research browser; sample metadata is provided by data.js."""

HTML = r"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>1000 Genomes · Local ancestry</title>
<meta name="description" content="Browse local ancestry, family relationships and Gnofix status for 3,202 high-coverage 1000 Genomes samples.">
<script>
try { document.documentElement.dataset.theme = localStorage.getItem('gnomix-theme') || (matchMedia('(prefers-color-scheme: dark)').matches ? 'dark' : 'light'); }
catch (_) { document.documentElement.dataset.theme = 'light'; }
</script>
<style>
:root { color-scheme:light; --bg:#fff; --ink:#20262d; --muted:#626c77; --rule:#dce1e6; --soft:#f6f8fa; --accent:#245777; --hover:#edf3f7; --link:#0b57d0; --head-hover:#e9edf1; }
:root[data-theme=dark] { color-scheme:dark; --bg:#15191e; --ink:#e1e6ec; --muted:#a4aeb9; --rule:#343d47; --soft:#1d232b; --accent:#93bedb; --hover:#232e39; --link:#7cacf8; --head-hover:#28313b; }
* { box-sizing:border-box; }
body { margin:0; background:var(--bg); color:var(--ink); font:14px/1.6 -apple-system,BlinkMacSystemFont,"Segoe UI",Helvetica,Arial,sans-serif; }
a { color:var(--accent); text-underline-offset:3px; }
button,input,select { font:inherit; }
button { cursor:pointer; }
:focus-visible { outline:2px solid var(--accent); outline-offset:3px; }
.wrap { width:100%; max-width:1240px; margin:auto; padding:0 32px; }
.top { display:flex; justify-content:space-between; align-items:center; gap:20px; min-height:70px; border-bottom:1px solid var(--rule); }
.project { font-size:13px; font-weight:600; letter-spacing:.02em; }
.top-right { display:flex; align-items:center; gap:24px; }
.top a { font-size:13px; text-decoration:none; }
.theme { display:flex; gap:0; border:1px solid var(--rule); border-radius:3px; padding:2px; }
.theme button { color:var(--muted); background:transparent; border:0; border-radius:2px; padding:3px 10px; font-size:12px; }
.theme button[aria-pressed=true] { color:var(--ink); background:var(--soft); font-weight:600; }
.intro { padding:38px 0 28px; }
h1 { font-size:32px; font-weight:600; letter-spacing:-.035em; line-height:1.2; margin:0 0 12px; }
.intro p { color:var(--muted); margin:0; max-width:76ch; }
.filters { display:grid; grid-template-columns:1fr 1.65fr 1fr; gap:20px; padding:24px 0; border-top:1px solid var(--rule); border-bottom:1px solid var(--rule); }
.field { display:flex; flex-direction:column; gap:7px; min-width:0; }
.field label { font-size:12px; font-weight:600; }
select,input[type=search] { width:100%; min-width:0; color:var(--ink); background:var(--bg); border:1px solid var(--rule); border-radius:3px; height:40px; padding:0 11px; }
select { padding-right:28px; text-overflow:ellipsis; }
input::placeholder { color:var(--muted); opacity:1; }
:disabled { opacity:.55; cursor:default; }
.filter-bottom { display:flex; justify-content:space-between; align-items:center; gap:16px; padding:15px 0; }
.checkbox { display:flex; align-items:center; gap:8px; cursor:pointer; font-size:13px; }
.checkbox input { width:15px; height:15px; margin:0; accent-color:var(--accent); }
.count { font-size:13px; color:var(--muted); font-variant-numeric:tabular-nums; }
main { min-height:390px; padding-bottom:48px; }
.empty { padding:60px 0 80px; text-align:center; }
.empty h2 { margin:0 0 6px; font-size:18px; font-weight:500; }
.empty p { color:var(--muted); margin:0; }
.pophead { padding:20px 0 16px; }
.pophead h2 { font-size:21px; font-weight:500; line-height:1.4; letter-spacing:-.015em; margin:0 0 5px; }
.pophead p { margin:0; color:var(--muted); font-size:13px; }
.phase-note { margin:0 0 20px; padding-left:13px; border-left:2px solid var(--rule); color:var(--muted); font-size:13px; max-width:100ch; }
.table-wrap { overflow-x:auto; border-top:1px solid var(--rule); }
table { border-collapse:collapse; width:100%; font-size:13px; font-variant-numeric:tabular-nums; white-space:nowrap; }
caption { text-align:left; color:var(--muted); font-size:12px; padding:0 0 12px; caption-side:bottom; padding-top:12px; }
th,td { text-align:left; padding:12px 13px; border-bottom:1px solid var(--rule); }
th:first-child,td:first-child { padding-left:12px; }
th:last-child,td:last-child { padding-right:12px; }
th { background:var(--soft); font-weight:600; font-size:12px; line-height:1.35; color:var(--muted); vertical-align:bottom; white-space:normal; }
th.num,td.num { text-align:right; }
th.num { min-width:96px; }
th.sortable { padding:0; }
.sort-head { display:flex; align-items:center; justify-content:flex-start; gap:7px; width:100%; color:inherit; background:none; border:0; padding:12px 13px; margin:0; font:inherit; text-align:left; }
th.sortable:first-child .sort-head { padding-left:12px; }
th.num .sort-head { justify-content:flex-end; text-align:right; }
.sort-head:hover { background:var(--head-hover); color:var(--ink); }
th[aria-sort=ascending] .sort-head, th[aria-sort=descending] .sort-head { color:var(--ink); }
.plain-head { display:inline-flex; align-items:center; min-height:20px; }
.sort-icon { flex:none; width:14px; height:20px; fill:currentColor; }
.sort-icon path { opacity:.35; }
.sort-head:hover .sort-icon path { opacity:.6; }
th[aria-sort=ascending] .sort-icon .up, th[aria-sort=descending] .sort-icon .down { opacity:1; }
th[aria-sort=ascending] .sort-icon .down, th[aria-sort=descending] .sort-icon .up { opacity:.12; }
tbody tr:hover { background:var(--hover); }
.sample-link { font:600 13px/1.5 ui-monospace,SFMono-Regular,Consolas,monospace; color:var(--link); text-decoration:underline; text-decoration-thickness:1px; text-underline-offset:3px; }
.sample-link:hover { text-decoration-thickness:2px; }
.below { color:var(--muted); }
.relationship { line-height:1.65; }
.phase { color:var(--muted); }
.phase.applied { color:var(--ink); }
footer { border-top:1px solid var(--rule); padding:22px 0 32px; color:var(--muted); font-size:12px; }
.footer-row { display:flex; justify-content:space-between; gap:24px; }
footer p { margin:0; max-width:88ch; }
footer a { white-space:nowrap; }
dialog { width:calc(100% - 48px); max-width:1400px; max-height:92vh; padding:0; border:1px solid var(--rule); border-radius:4px; background:var(--bg); color:var(--ink); }
dialog::backdrop { background:rgba(0,0,0,.65); }
.dlg-head { position:sticky; top:0; background:var(--bg); display:flex; justify-content:space-between; align-items:start; gap:20px; padding:20px 24px; border-bottom:1px solid var(--rule); z-index:1; }
.dlg-head h2 { margin:0 0 3px; font:600 18px/1.4 ui-monospace,SFMono-Regular,Consolas,monospace; }
.dlg-head p { font-size:13px; color:var(--muted); margin:0; }
.close { border:1px solid var(--rule); border-radius:3px; background:var(--bg); color:var(--ink); padding:5px 12px; }
.dlg-body { padding:20px 24px; }
.dlg-body img { width:100%; height:auto; display:block; background:#fff; }
.kname { font-size:13px; color:var(--muted); margin:24px 0 8px; }
.kname:first-child { margin-top:0; }
.dlg-related { border-top:1px solid var(--rule); margin-top:24px; padding-top:16px; }
.dlg-related h3 { margin:0 0 8px; font-size:14px; font-weight:600; }
.martin-link { display:inline-block; margin-top:16px; font-size:13px; }
[hidden] { display:none !important; }
@media(max-width:760px) {
 .wrap { padding:0 20px; } .top { min-height:64px; } .top-right { gap:12px; } .top-right>a { display:none; }
 .intro { padding:28px 0 24px; } h1 { font-size:27px; }
 .filters { grid-template-columns:1fr; gap:16px; padding:20px 0; } .filter-bottom { align-items:start; }
 .pophead h2 { font-size:18px; } .footer-row { flex-direction:column; gap:12px; }
 dialog { width:calc(100% - 20px); } .dlg-head,.dlg-body { padding:16px; }
}
</style>
</head>
<body>
<div class="wrap">
<header>
 <div class="top">
  <span class="project">1000 Genomes / Gnomix</span>
  <div class="top-right"><a href="REPO_URL">Code &amp; methods</a>
   <div class="theme" role="group" aria-label="Color theme">
    <button id="light" aria-pressed="false">Light</button><button id="dark" aria-pressed="false">Dark</button>
   </div>
  </div>
 </div>
 <div class="intro">
  <h1>Local ancestry of the 1000 Genomes panel</h1>
  <p>Global ancestry and karyograms for 3,202 high-coverage samples. Select a superpopulation and population to browse the results. Select a sample ID to view its karyogram.</p>
 </div>
</header>
<section aria-label="Sample filters">
 <div class="filters">
  <div class="field"><label for="sup">Superpopulation</label><select id="sup"><option value="">Select superpopulation</option></select></div>
  <div class="field"><label for="pop">Population</label><select id="pop" disabled><option value="">Select population</option></select></div>
  <div class="field"><label for="q">Sample ID</label><input type="search" id="q" placeholder="Search within population" disabled autocomplete="off" spellcheck="false"></div>
 </div>
 <div class="filter-bottom">
  <label class="checkbox"><input type="checkbox" id="hide-trios" aria-describedby="trio-help">Hide trios</label>
  <span class="count" id="count" role="status" aria-live="polite"></span>
 </div>
 <span id="trio-help" hidden>Hides both parents and children in trios. Duo members remain visible.</span>
</section>
<main id="list"></main>
<footer><div class="footer-row">
 <p>Pretrained <a href="https://github.com/AI-sandbox/gnomix">Gnomix</a> calls. Gnofix is used for every sample of the six admixed populations except trio children, who retain their published phase (made with both parents).</p>
 <a href="RELEASE_URL">Download all outputs ↗</a>
</div></footer>
</div>
<dialog id="dlg" aria-labelledby="dlg-id" aria-describedby="dlg-sub">
 <div class="dlg-head"><div><h2 id="dlg-id"></h2><p id="dlg-sub"></p></div><button class="close" id="dlg-close" autofocus>Close</button></div>
 <div class="dlg-body" id="dlg-body"></div>
</dialog>
<script src="data.js"></script>
<script>
const D = window.DATA, byId = Object.fromEntries(D.samples.map(s => [s.id, s]));
// European and West Asian are shown separately for the EUR and SAS superpopulations, combined elsewhere.
const SPLIT_EUR_WAS = new Set(['EUR','SAS']);
const COLUMN_ORDER = ['AFR','AHG','EUR','WAS','SAS','EAS','NAT','OCE'];
function columnsFor(sup) {
 const was = D.anc.findIndex(a => a.code === 'WAS');
 return COLUMN_ORDER.map(code => [D.anc.find(a => a.code === code), D.anc.findIndex(a => a.code === code)]).flatMap(([a,i]) => {
  if (SPLIT_EUR_WAS.has(sup)) return [{label:a.name, codes:a.code, indices:[i]}];
  if (a.code === 'WAS') return [];
  if (a.code === 'EUR') return [{label:'European + West Asian', codes:'EUR + WAS', indices:[i,was]}];
  return [{label:a.name, codes:a.code, indices:[i]}];
 });
}
function fraction(s,col) { return col.indices.reduce((total,index) => total + s.a[index], 0); }
const $ = id => document.getElementById(id);
const supSel = $('sup'), popSel = $('pop');
const admixed = new Set(['ACB','ASW','CLM','MXL','PEL','PUR']);
const state = { sort:'id', desc:false };
const escapeHTML = value => String(value).replace(/[&<>"']/g, c => ({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c]));
function family(s) {
 if (s.family) return s.family;
 const roles = [];
 if (s.par.length) roles.push(s.par.length === 2 ? 'Trio child' : 'Duo child');
 if (s.kids.some(k => byId[k]?.par.length === 2)) roles.push('Trio parent');
 if (s.kids.some(k => byId[k]?.par.length === 1)) roles.push('Duo parent');
 return roles.length ? roles : ['Unrelated'];
}
function phase(s) { return s.fix ? 'Applied' : 'Not applied'; }
function phaseDetail(s) {
 if (s.fix) return 'Gnofix applied';
 return s.par.length === 2 ? 'Gnofix not applied · trio child, phased with both parents' : 'Gnofix not applied · population policy';
}
function setTheme(theme) {
 document.documentElement.dataset.theme = theme;
 ['light','dark'].forEach(t => $(t).setAttribute('aria-pressed', String(t === theme)));
 try { localStorage.setItem('gnomix-theme', theme); } catch (_) {}
}
['light','dark'].forEach(t => $(t).onclick = () => setTheme(t));
setTheme(document.documentElement.dataset.theme);
Object.entries(D.sup).forEach(([code,name]) => supSel.add(new Option(name + ' (' + code + ')', code)));
function fillPops() {
 popSel.replaceChildren(new Option('Select population', ''));
 D.pops.filter(p => p.sup === supSel.value).forEach(p => popSel.add(new Option(p.pop + ' — ' + p.name, p.pop)));
 popSel.disabled = !supSel.value;
}
const SORT_ICON = '<svg class="sort-icon" viewBox="0 0 14 20" width="14" height="20" aria-hidden="true"><path class="up" d="M7 1.5 13.5 8.5H.5z"/><path class="down" d="M7 18.5.5 11.5h13z"/></svg>';
function heading(key, label, numeric=false, title='') {
 const sorted = state.sort === key ? (state.desc ? 'descending' : 'ascending') : 'none';
 return '<th scope="col" class="sortable' + (numeric ? ' num' : '') + '" aria-sort="' + sorted + '"><button class="sort-head" data-key="' + key + '" title="' + escapeHTML(title || 'Sort by ' + label) + '"><span>' + escapeHTML(label) + '</span>' + SORT_ICON + '</button></th>';
}
function render() {
 const pop = popSel.value, meta = D.pops.find(p => p.pop === pop);
 $('q').disabled = !pop;
 if (!meta) {
  $('count').textContent = '';
  $('list').innerHTML = '<div class="empty"><h2>Select a population</h2><p>Choose a superpopulation, then a population, to view its samples.</p></div>';
  return;
 }
 const population = D.samples.filter(s => s.pop === pop);
 const q = $('q').value.trim().toUpperCase();
 const cols = columnsFor(meta.sup);
 const rows = population.filter(s => (!q || s.id.toUpperCase().includes(q)) && (!$('hide-trios').checked || !family(s).some(r => r.startsWith('Trio'))));
 rows.sort((x,y) => {
  const c = state.sort === 'id' ? x.id.localeCompare(y.id) : fraction(x,cols[+state.sort.slice(1)]) - fraction(y,cols[+state.sort.slice(1)]);
  return (state.desc ? -c : c) || x.id.localeCompare(y.id);
 });
 $('count').textContent = rows.length + ' of ' + population.length + ' samples';
 const shown = cols.map((_,i) => i).filter(i => population.some(s => fraction(s,cols[i]) >= .002));
 const applied = population.filter(s => s.fix).length;
 const protectedCount = population.filter(s => !s.fix && s.par.length === 2).length;
 let html = '<section aria-labelledby="pop-title"><div class="pophead"><h2 id="pop-title">' + escapeHTML(meta.name) + '</h2><p>' + escapeHTML(meta.pop + ' · ' + D.sup[meta.sup]) + ' · ' + population.length + ' samples' + (admixed.has(pop) ? ' · Gnofix applied to ' + applied + '; ' + protectedCount + ' trio children excluded' : '') + '</p></div>';
 if (admixed.has(pop)) html += '<p class="phase-note">Trio children retain their published phase: it was made with both parents and is near-exact, so Gnofix could only add errors. All other samples get Gnofix, including trio parents and duo members, whose published phase is no better than that of unrelated samples.</p>';
 if (!rows.length) html += '<div class="empty"><h2>No samples match</h2><p>Change the sample ID search or uncheck Hide trios.</p></div>';
 else {
  html += '<div class="table-wrap" role="region" aria-label="Sample ancestry table" tabindex="0"><table><caption>Global ancestry (% of genetic length). Values below 0.2% are omitted (—). Select a sample ID to open its karyogram; select a column heading to sort.</caption><thead><tr>' + heading('id','Sample ID') + shown.map(i => heading('a'+i, cols[i].label, true, 'Sort by ' + cols[i].label + ' (' + cols[i].codes + ')')).join('') + '<th scope="col"><span class="plain-head">Relationship</span></th><th scope="col"><span class="plain-head">Gnofix</span></th></tr></thead><tbody>';
  rows.forEach(s => {
   html += '<tr><td><a class="sample-link" href="#sample=' + encodeURIComponent(s.id) + '" data-id="' + escapeHTML(s.id) + '">' + escapeHTML(s.id) + '</a></td>' + shown.map(i => {
    const value = fraction(s,cols[i]);
    return '<td class="num' + (value < .002 ? ' below' : '') + '">' + (value >= .002 ? (100*value).toFixed(1) + '%' : '—') + '</td>';
   }).join('') + '<td class="relationship">' + family(s).map(escapeHTML).join('<br>') + '</td><td class="phase' + (s.fix ? ' applied' : '') + '" title="' + escapeHTML(phaseDetail(s)) + '">' + phase(s) + '</td></tr>';
  });
  html += '</tbody></table></div>';
 }
 $('list').innerHTML = html + '</section>';
 document.querySelectorAll('[data-key]').forEach(b => b.onclick = () => {
  const key = b.dataset.key;
  state.desc = state.sort === key ? !state.desc : key !== 'id';
  state.sort = key; render();
  document.querySelector('[data-key="' + key + '"]').focus();
 });
 document.querySelectorAll('[data-id]').forEach(a => a.onclick = e => { e.preventDefault(); show(a.dataset.id); });
}
function karyoPath(s) { return 'karyograms/' + encodeURIComponent(s.pop) + '/' + encodeURIComponent(s.id) + '.png'; }
function imageHTML(s) { return '<img src="' + karyoPath(s) + '" alt="Local ancestry karyogram of ' + escapeHTML(s.id) + '" loading="lazy">'; }
function show(id) {
 const s = byId[id];
 if (!s) return;
 const meta = D.pops.find(p => p.pop === s.pop);
 $('dlg-id').textContent = id;
 $('dlg-sub').textContent = meta.pop + ' · ' + family(s).join(', ') + ' · ' + phaseDetail(s);
 let body = imageHTML(s);
 if (s.m) body += '<a class="martin-link" target="_blank" rel="noopener" href="' + escapeHTML(D.martin_url.replace('{}',id)) + '">Compare with Martin et al. (2017) karyogram (PDF) ↗</a>';
 const coParents = [...new Set(s.kids.flatMap(k => byId[k] ? byId[k].par : []))].filter(p => p !== s.id);
 const relatives = [...s.par.map(p => [p,'Parent']), ...coParents.map(p => [p,'Other parent']), ...s.kids.map(k => [k,'Child'])];
 if (relatives.length) {
  body += '<section class="dlg-related"><h3>Family karyograms</h3>';
  relatives.forEach(([rid,rel]) => {
   const r = byId[rid];
   if (r) body += '<p class="kname">' + escapeHTML(rid + ' · ' + rel + ' · ' + family(r).join(', ') + ' · ' + phaseDetail(r)) + '</p>' + imageHTML(r);
  });
  body += '</section>';
 }
 $('dlg-body').innerHTML = body;
 $('dlg').showModal();
 $('dlg').scrollTop = 0;
}
supSel.onchange = () => { fillPops(); $('q').value = ''; state.sort = 'id'; state.desc = false; render(); };
popSel.onchange = () => { $('q').value = ''; state.sort = 'id'; state.desc = false; render(); };
$('q').oninput = render;
$('hide-trios').onchange = render;
$('dlg-close').onclick = () => $('dlg').close();
$('dlg').onclick = e => {
 if (e.target !== $('dlg')) return;
 const rect = $('dlg').getBoundingClientRect();
 if (e.clientX < rect.left || e.clientX > rect.right || e.clientY < rect.top || e.clientY > rect.bottom) $('dlg').close();
};
fillPops(); render();
</script>
</body>
</html>
""".replace("REPO_URL", "https://github.com/human-genomics/gnomix-1000g").replace(
    "RELEASE_URL", "https://github.com/human-genomics/gnomix-1000g/releases/tag/v1.1.0")
