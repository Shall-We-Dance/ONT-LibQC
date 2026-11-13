/* The same exact base-coverage integration as custom_regions.py; no network calls. */
const CustomRegions = (() => {
// Do not infer file type from a .bed suffix: directories may have that suffix too.
async function checkedBedFile(input, maxBytes = 10 * 1024 ** 2) {
  const file = input.files?.[0];
  if ([...(input.webkitEntries || [])].some(entry => entry.isDirectory) || ['inode/directory', 'application/x-directory'].includes(file?.type)) {
    throw Error('Choose a BED file inside the folder, not the folder itself.');
  }
  if (!file) throw Error('Choose a BED file first.');
  if (file.size > maxBytes) throw Error('Use a BED file of at most ' + maxBytes / 1024 ** 2 + ' MiB.');
  // Dialog selections often have no filesystem entry. A read probe rejects
  // directory-like or unreadable selections before parsing or upload.
  try { await file.slice(0, 1).arrayBuffer(); }
  catch { throw Error('Cannot read this selection. Open the folder and choose a readable BED file inside it.'); }
  if (!file.size) throw Error('The BED file is empty. Choose a file containing regions.');
  return file;
}
function bedDropIsDirectory(event) {
  return [...(event.dataTransfer?.items || [])].some(item => (item.getAsEntry?.() || item.webkitGetAsEntry?.())?.isDirectory);
}
function parseBed(text, sizes, primary) {
  const rows=[],seen=new Set(),keep=new Set(primary);const renamed=new Map();let skipped=0,duplicates=0,renamed_rows=0,unknown=[];
  for (const [i,line] of text.split(/\r?\n/).entries()) {
    if (!line.trim() || /^(#|track|browser)/.test(line)) continue;
    const v=line.trim().split(/\s+/),s=Number(v[1]),e=Number(v[2]);
    if(v.length<3 || !/^\d+$/.test(v[1]) || !/^\d+$/.test(v[2]) || !Number.isSafeInteger(s) || !Number.isSafeInteger(e)) throw Error(`BED line ${i+1}: integer coordinates required`);
    if(s<0 || e<=s) throw Error(`BED line ${i+1}: invalid interval`);
    if(!Object.prototype.hasOwnProperty.call(sizes,v[0])){
      const alternate=v[0].startsWith('chr')?v[0].slice(3):'chr'+v[0];
      if(Object.prototype.hasOwnProperty.call(sizes,alternate)){renamed.set(v[0],alternate);renamed_rows++;v[0]=alternate;}
    }
    if(!Object.prototype.hasOwnProperty.call(sizes,v[0])){unknown.push({line:i+1,chrom:v[0],start:s,end:e,reason:'chromosome_absent_from_reference'});continue;}
    if(s<0 || e<=s || e>sizes[v[0]]) throw Error(`BED line ${i+1}: invalid or out-of-reference interval`);
    const strand=v.length>=6?v[5]:'.';if(!['+','-','.'].includes(strand)) throw Error(`BED line ${i+1}: invalid strand`);
    if(!keep.has(v[0])) {skipped++;continue;}
    const key=JSON.stringify([v[0],s,e,strand]);if(seen.has(key)){duplicates++;continue;}seen.add(key);
    rows.push({chrom:v[0],start:s,end:e,strand});if(rows.length>100000) throw Error('Use at most 100,000 regions per browser profile');
  }
  if(!rows.length) throw Error('No regions on configured primary chromosomes');
  return {rows,skipped,duplicates,unknown,renamed_rows,renamed_chromosomes:Object.fromEntries(renamed)};
}
function makeTrack(records) {
  const events=new Map();for(const [s,e] of records){events.set(s,(events.get(s)||0)+1);events.set(e,(events.get(e)||0)-1);}
  const x=[...events.keys()].sort((a,b)=>a-b),y=[],area=[];let h=0,a=0;
  for(let i=0;i<x.length;i++){if(i)a+=(x[i]-x[i-1])*h;h+=events.get(x[i]);area.push(a);y.push(h);}
  return {x,y,area};
}
function integral(t,v){if(!t || !t.x.length || v<t.x[0])return 0;let lo=0,hi=t.x.length;while(lo<hi){let m=(lo+hi)>>1;if(t.x[m]<=v)lo=m+1;else hi=m;}let i=lo-1;return t.area[i]+(v-t.x[i])*t.y[i];}
function resolveParameters(rows,params={}) {
  const p={mode:'auto',upstream_bp:5000,downstream_bp:5000,flank_bins:50,body_bins:100,reference_bins:100,...params};
  if(!['auto','reference_point','region'].includes(p.mode))throw Error('Custom region mode must be auto, reference_point or region');
  p.requested_mode=p.mode==='auto'?'auto':(p.requested_mode||p.mode);
  if(p.mode==='auto')p.mode=rows.every(r=>r.end-r.start===1)?'reference_point':'region';
  return p;
}
function compute(records,rows,sizes,p){
  p=resolveParameters(rows,p);
  const grouped=new Map();for(const [c,s,e,blocks] of records){if(!grouped.has(c))grouped.set(c,[]);for(const block of (blocks||[[s,e]]))grouped.get(c).push(block);}
  const tracks=new Map([...grouped].map(([c,rs])=>[c,makeTrack(rs)]));
  const nbin=p.mode==='reference_point'?p.reference_bins:2*p.flank_bins+p.body_bins,sum=Array(nbin).fill(0),square=Array(nbin).fill(0),counts=Array(nbin).fill(0),factor=records.length?1e6/records.length:0;
  for(const r of rows){
    const left=r.strand==='-'?p.downstream_bp:p.upstream_bp,right=r.strand==='-'?p.upstream_bp:p.downstream_bp,edges=[];
    if(p.mode==='reference_point'){
      const center=(r.start+r.end)/2;
      for(let i=0;i<=p.reference_bins;i++)edges.push(center-left+i*(left+right)/p.reference_bins);
    }else{
    for(let i=0;i<p.flank_bins;i++)edges.push(r.start-left+i*left/p.flank_bins);
    for(let i=0;i<p.body_bins;i++)edges.push(r.start+i*(r.end-r.start)/p.body_bins);
    for(let i=0;i<=p.flank_bins;i++)edges.push(r.end+i*right/p.flank_bins);
    }
    for(let i=0;i<nbin;i++){
      const s=Math.max(0,Math.min(sizes[r.chrom],edges[i])),e=Math.max(0,Math.min(sizes[r.chrom],edges[i+1]));if(e<=s)continue;
      const j=r.strand==='-'?nbin-1-i:i,v=(integral(tracks.get(r.chrom),e)-integral(tracks.get(r.chrom),s))/(e-s)*factor;
      sum[j]+=v;square[j]+=v*v;counts[j]++;
    }
  }
  return {mean:counts.map((n,i)=>n&&records.length?sum[i]/n:null),sem:counts.map((n,i)=>n>1&&records.length?Math.sqrt(Math.max(0,square[i]-sum[i]*sum[i]/n)/(n*(n-1))):null),regions_per_bin:counts,eligible_alignments:records.length};
}
return {parseBed,compute,resolveParameters,checkedBedFile,bedDropIsDirectory};
})();
if(typeof module!=='undefined')module.exports=CustomRegions;
if(typeof document!=='undefined')(() => {
const payload=JSON.parse(document.getElementById('custom-payload').textContent),message=document.getElementById('custom-message'),output=document.getElementById('custom-output'),downloads=document.getElementById('custom-downloads');let urls=[];
const localResult=document.getElementById('custom-local'),currentSource=document.getElementById('custom-current-source');
function clearResult(){output.replaceChildren();downloads.replaceChildren();urls.forEach(URL.revokeObjectURL);urls=[];if(localResult)localResult.hidden=true;if(currentSource)currentSource.textContent='';}
const bedInput=document.getElementById('custom-bed'),bedName=document.getElementById('custom-bed-name'),runButton=document.getElementById('custom-run');
let selectionVersion=0;
// Let Safari navigate folders normally; validate the selected file afterward.
bedInput.removeAttribute('accept');
function showBedName(){const name=bedInput.files[0]?.name||'No new BED selected';if(bedName)bedName.textContent=name;bedInput.closest('.custom-file-picker').title=name;}
bedInput.addEventListener('change',async()=>{
  const version=++selectionVersion;
  document.getElementById('custom-assembly').checked=false;
  clearResult();showBedName();
  if(!bedInput.files.length){message.textContent='';runButton.disabled=false;return;}
  message.textContent='Checking BED file…';runButton.disabled=true;
  try{
    await CustomRegions.checkedBedFile(bedInput);
    if(version===selectionVersion)message.textContent='Confirm coordinates and click Draw metaplot to analyze the selected BED.';
  }catch(e){
    if(version===selectionVersion){bedInput.value='';showBedName();message.textContent=e.message;}
  }finally{if(version===selectionVersion)runButton.disabled=false;}
});
bedInput.addEventListener('drop',event=>{
  if(!CustomRegions.bedDropIsDirectory(event))return;
  event.preventDefault();++selectionVersion;bedInput.value='';showBedName();
  document.getElementById('custom-assembly').checked=false;clearResult();runButton.disabled=false;
  message.textContent='Choose a BED file inside the folder, not the folder itself.';
});
const escape=s=>String(s).replace(/[&<>"']/g,c=>({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c]));
runButton.addEventListener('click',async()=>{
  const version=++selectionVersion;runButton.disabled=true;
  clearResult();
  try {
    if(!document.getElementById('custom-assembly').checked)throw Error(`Confirm that the BED uses ${payload.reference} coordinates`);
    const file=await CustomRegions.checkedBedFile(bedInput);if(version!==selectionVersion)return;
    if(localResult)localResult.hidden=false;if(currentSource)currentSource.textContent='BED: '+file.name;
    message.textContent='Calculating…';await new Promise(resolve=>setTimeout(resolve,0));
    const text=await file.text();if(version!==selectionVersion)return;
    const bed=CustomRegions.parseBed(text,payload.sizes,payload.primary),scope=payload.scopes[document.getElementById('custom-scope').value],params=CustomRegions.resolveParameters(bed.rows,{...payload.parameters,mode:document.getElementById('custom-mode').value}),p=CustomRegions.compute(scope.records,bed.rows,payload.sizes,params);
    const point=params.mode==='reference_point',modeLabel=point?'Reference point (interval midpoint)':'Region (scaled body)';
    message.textContent=`${params.requested_mode==='auto'?'Auto → ':''}${modeLabel} · ${bed.rows.length} regions · ${p.eligible_alignments.toLocaleString()} eligible alignments · ${bed.skipped} non-primary rows skipped · ${bed.duplicates} duplicates removed · ${bed.unknown.length} regions on absent contigs excluded${bed.unknown.length ? ": "+[...new Set(bed.unknown.map(x=>x.chrom))].join(", ") : ""}`;
    if(!p.eligible_alignments){message.textContent+=' · No eligible nuclear alignments; profile unavailable.';return;}
    // Match the report's TSS–TES plots: boxed axes, round ticks and a quiet SEM band.
    const w=820,h=440,l=80,t=52,pw=712,ph=298,n=p.mean.length,color='#000080',ink='#202b36';
    const maximum=Math.max(0,...p.mean.map((v,i)=>v==null?0:v+(p.sem[i]||0)))||1;
    const magnitude=10**Math.floor(Math.log10(maximum/4));
    const step=[1,2,2.5,5,10].map(v=>v*magnitude).reduce((a,b)=>Math.abs(Math.log(a/(maximum/4)))<Math.abs(Math.log(b/(maximum/4)))?a:b);
    const ylim=Math.ceil(maximum*1.06/step)*step,x=i=>l+(i+.5)/n*pw,y=v=>t+ph-v/ylim*ph;
    let parts=[`<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 ${w} ${h}" role="img" aria-label="Custom region metaplot" style="font-family:Arial,Helvetica,sans-serif"><title>${escape(scope.name)} · Custom region metaplot</title><desc>${escape(payload.sample+' · '+file.name)} · ${bed.rows.length} regions · ${escape(payload.reference)} · ${modeLabel}. Equal region weights; shaded band = SEM across regions.</desc><rect width="100%" height="100%" fill="white"/><defs><clipPath id="custom-profile-frame"><rect x="${l}" y="${t}" width="${pw}" height="${ph}"/></clipPath></defs><text x="${l+pw/2}" y="29" text-anchor="middle" font-size="20" fill="${ink}">${escape(scope.name)}</text>`];
    for(let i=0;i<=Math.round(ylim/step);i++){const v=i*step;parts.push(`<line x1="${l-5}" x2="${l}" y1="${y(v)}" y2="${y(v)}" stroke="${ink}"/><text x="${l-11}" y="${y(v)+5}" text-anchor="end" font-size="14" fill="${ink}">${Number(v.toPrecision(5))}</text>`);}
    let runs=[],run=[];for(let i=0;i<n;i++){if(p.mean[i]==null){if(run.length)runs.push(run);run=[];}else run.push(i);}if(run.length)runs.push(run);
    for(const indices of runs){const upper=indices.map(i=>`${x(i)},${y(p.mean[i]+(p.sem[i]||0))}`),lower=[...indices].reverse().map(i=>`${x(i)},${y(Math.max(0,p.mean[i]-(p.sem[i]||0)))}`);parts.push(`<polygon points="${upper.concat(lower).join(' ')}" fill="${color}" fill-opacity=".20" clip-path="url(#custom-profile-frame)"/><polyline points="${indices.map(i=>`${x(i)},${y(p.mean[i])}`).join(' ')}" stroke="${color}" stroke-width="2.4" fill="none" clip-path="url(#custom-profile-frame)"/>`);}
    const center=params.upstream_bp/(params.upstream_bp+params.downstream_bp);
    const ticks=[[0,`−${params.upstream_bp/1000} kb`],...(point?[[center,'Reference point']]:[[params.flank_bins/n,'Region start'],[(params.flank_bins+params.body_bins)/n,'Region end']]),[1,`+${params.downstream_bp/1000} kb`]];
    for(const [f,label] of ticks)parts.push(`<line x1="${l+pw*f}" x2="${l+pw*f}" y1="${t+ph}" y2="${t+ph+5}" stroke="${ink}"/><text x="${l+pw*f}" y="${t+ph+28}" text-anchor="${f===1?'end':'middle'}" font-size="15" fill="${ink}">${label}</text>`);
    parts.push(`<rect x="${l}" y="${t}" width="${pw}" height="${ph}" fill="none" stroke="${ink}" stroke-width="1.2"/><text x="22" y="${t+ph/2}" transform="rotate(-90 22 ${t+ph/2})" text-anchor="middle" font-size="15" fill="${ink}">Mean coverage (CPM)</text><line x1="365" x2="389" y1="417" y2="417" stroke="${color}" stroke-width="2.4"/><text x="397" y="421" font-size="13" fill="${ink}">regions</text></svg>`);
    const svg=parts.join('');output.innerHTML=svg;
    const tsv='bin\tmean_cpm\tsem_across_regions\tregions_in_bin\n'+p.mean.map((v,i)=>[i,v??'NA',p.sem[i]??'NA',p.regions_per_bin[i]].join('\t')).join('\n')+'\n';
    for(const [label,content,mime,ext] of [['Download SVG',svg,'image/svg+xml','svg'],['Download data',tsv,'text/tab-separated-values','tsv'],['Download provenance',JSON.stringify({sample:payload.sample,scope:scope.name,assembly:payload.reference,bed_name:file.name,parameters:params,retained_regions:bed.rows.length,renamed_chromosomes:bed.renamed_chromosomes,renamed_rows:bed.renamed_rows,excluded_unknown_contigs:bed.unknown,skipped_nonprimary:bed.skipped,duplicate_intervals_removed:bed.duplicates,eligible_alignments:p.eligible_alignments,normalization:'exact mean base coverage; CPM by scope primary nuclear alignments; equally weighted regions'},null,2),'application/json','json']]){
      const a=document.createElement('a'),url=URL.createObjectURL(new Blob([content],{type:mime}));urls.push(url);a.href=url;a.download=payload.sample+'.'+scope.name+'.custom_regions.'+ext;a.textContent=label;downloads.appendChild(a);
    }
  }catch(e){if(version===selectionVersion)message.textContent=e.message;}
  finally{if(version===selectionVersion)runButton.disabled=false;}
});
})();
