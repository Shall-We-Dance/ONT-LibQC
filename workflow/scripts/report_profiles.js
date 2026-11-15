
(()=>{
 const root=document.currentScript.closest('.profile-comparison');
 if(!root)return;
 const buttons=[...root.querySelectorAll('.profile-comparison-legend button[data-scope]')];
 const traces=[...root.querySelectorAll('.profile-size-chart [data-scope]')];
 let pinned=null, pointerType=null;
 function show(scope){
  for(const node of [...buttons,...traces])node.classList.toggle('is-dimmed',!!scope&&node.dataset.scope!==scope);
  for(const b of buttons){
   b.setAttribute('aria-pressed',String(b.dataset.scope===pinned));
   b.classList.toggle('is-highlighted',b.dataset.scope===scope);
  }
 }
 for(const b of buttons){
  b.addEventListener('pointerenter',e=>{if(e.pointerType!=='touch'&&matchMedia('(hover:hover)').matches){pinned=null;show(b.dataset.scope)}});
  b.addEventListener('pointerleave',e=>{if(e.pointerType!=='touch'&&matchMedia('(hover:hover)').matches){pinned=null;show(null)}});
  b.addEventListener('focus',()=>{if(b.matches(':focus-visible'))show(b.dataset.scope)});
  b.addEventListener('blur',()=>show(pinned));
  b.addEventListener('pointerdown',e=>{pointerType=e.pointerType});
  b.addEventListener('keydown',()=>{pointerType=null});
  b.addEventListener('click',e=>{
   const touch=pointerType==='touch'||e.pointerType==='touch'||(e.detail>0&&matchMedia('(hover:none)').matches);
   pointerType=null;
   if(touch){pinned=pinned===b.dataset.scope?null:b.dataset.scope;show(pinned)}
   else{pinned=null;show(b.dataset.scope)}
  });
 }
 for(const line of root.querySelectorAll('.profile-size-chart polyline[data-scope]')){
  line.addEventListener('pointerenter',e=>{if(e.pointerType!=='touch'&&matchMedia('(hover:hover)').matches){pinned=null;show(line.dataset.scope)}});
  line.addEventListener('pointerleave',e=>{if(e.pointerType!=='touch'&&matchMedia('(hover:hover)').matches){pinned=null;show(null)}});
 }
 function reset(){pinned=null;pointerType=null;show(null)}
 root.addEventListener('pointerleave',e=>{if(e.pointerType!=='touch'&&matchMedia('(hover:hover)').matches)reset()});
 root.addEventListener('pointercancel',reset);
 root.addEventListener('keydown',e=>{if(e.key==='Escape')reset()});
 document.addEventListener('pointerdown',e=>{if(!root.contains(e.target))reset()});
 window.addEventListener('blur',reset);
 document.addEventListener('visibilitychange',()=>{if(document.hidden)reset()});
})();
