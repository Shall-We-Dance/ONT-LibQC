/* Per-file editor for the physical P5-to-P7 library strand. */
const LibraryStructure = (() => {
  const positions=['before_i5','after_i5','before_insert','after_insert','before_i7','after_i7'];
  const labels=['P5 / i5','i5 / Adapter','Adapter / Insert','Insert / Adapter','Adapter / i7','i7 / P7'];
  const blocks=[['P5','handle'],['i5','index'],['Adapter','anchor','i5 side'],['Insert','insert'],['Adapter','anchor','i7 side'],['i7','index'],['P7','handle']];
  const defaults=()=>({mode:'auto',adapter_family:'auto',p5_sequence:null,p7_sequence:null,i5_length:null,i7_length:null,adapter_i5_sequence:null,adapter_i7_sequence:null,umis:[]});
  function create(catalog,onApply) {
    let state=defaults(),selected=null,dragging=null,nextKey=0;
    const keys=new WeakMap(),rows=new Map();
    const key=umi=>{if(!keys.has(umi))keys.set(umi,String(++nextKey));return keys.get(umi);};
    const sorted=()=>[...state.umis].sort((a,b)=>positions.indexOf(a.position)-positions.indexOf(b.position));
    const motion=()=>!matchMedia('(prefers-reduced-motion: reduce)').matches;
    const limits=catalog||{p5_sequence:'AATGATACGGCGACCACCGAGATCTACAC',p7_sequence:'ATCTCGTATGCCGTCTTCTGCTTG'};
    const root=node('section',{class:'library-structure','aria-label':'Library structure'});
    const auto=node('button',{type:'button','aria-pressed':'true'},'Auto');
    const custom=node('button',{type:'button','aria-pressed':'false'},'Customize');
    const description=node('p',{class:'structure-description'});
    const track=node('div',{class:'structure-track','aria-label':'Library sequence from P5 to P7'});
    const scroller=node('div',{class:'structure-scroll',tabindex:'0','aria-label':'Library structure diagram; scroll horizontally to view all blocks'},track);
    const hint=node('p',{class:'structure-hint'});
    const umiTitle=node('h4',{},'UMIs');
    const list=node('div',{class:'structure-umi-list'});
    const umiPanel=node('div',{class:'structure-umi-panel',hidden:true},node('div',{class:'structure-umi-heading'},umiTitle,node('span',{},'Fixed length · each position used once')),list);
    const live=node('span',{class:'visually-hidden',role:'status','aria-live':'polite'});
    const notice=node('p',{class:'structure-copy-note',role:'status',hidden:true});
    const error=node('p',{class:'error',role:'alert',hidden:true});
    const add=node('button',{type:'button',class:'structure-add'},'+ Add UMI');
    const apply=node('button',{type:'button',class:'structure-apply',hidden:true},'Apply structure to all');
    const advanced=node('details',{class:'structure-settings'},node('summary',{},'Adapter & index settings'));
    const family=select([['auto','Auto · TruSeq / Nextera'],['truseq','TruSeq'],['nextera','Nextera']],'auto');
    const i5=node('input',{type:'number',min:'0',max:'24',step:'1',placeholder:'Auto'});
    const i7=node('input',{type:'number',min:'0',max:'24',step:'1',placeholder:'Auto'});
    const p5=node('input',{type:'text',placeholder:limits.p5_sequence,pattern:'[ACGTacgt]{16,80}',maxlength:'80',spellcheck:'false',autocomplete:'off'});
    const p7=node('input',{type:'text',placeholder:limits.p7_sequence,pattern:'[ACGTacgt]{16,80}',maxlength:'80',spellcheck:'false',autocomplete:'off'});
    const adapter5=node('input',{type:'text',placeholder:'Use preset / auto-detect',pattern:'[ACGTacgt]{16,80}',maxlength:'80',spellcheck:'false',autocomplete:'off'});
    const adapter7=adapter5.cloneNode();
    p5.className=p7.className=adapter5.className=adapter7.className='structure-sequence';
    const fields={adapter_family:family,i5_length:i5,i7_length:i7,p5_sequence:p5,p7_sequence:p7,adapter_i5_sequence:adapter5,adapter_i7_sequence:adapter7};
    advanced.append(node('div',{class:'structure-settings-grid'},control('Adapter preset',family),control('i5 length (bp)',i5),control('i7 length (bp)',i7),node('div',{class:'structure-setting-wide'},control('Adapter · i5 side · 5′ → 3′',adapter5)),node('div',{class:'structure-setting-wide'},control('Adapter · i7 side · 5′ → 3′',adapter7)),node('div',{class:'structure-setting-wide'},control('P5 sequence · 5′ → 3′',p5)),node('div',{class:'structure-setting-wide'},control('P7 sequence · 5′ → 3′',p7))),node('p',{class:'field-note'},'Enter 16–80 A/C/G/T bases per sequence, following the P5 → P7 strand above. Enter the i7-side adapter as it appears on that strand, not the reverse-complement primer. A custom adapter overrides its side of the preset; leave a side blank to use the preset / auto-detect. Blank index length means Auto; 0 means no index.'));
    root.append(node('div',{class:'structure-heading'},node('h3',{},'Library structure'),node('div',{class:'structure-modes',role:'group','aria-label':'Library structure mode'},auto,custom)),description,scroller,hint,umiPanel,node('div',{class:'structure-footer'},add,apply),advanced,notice,error,live);
    function announce(text){live.textContent=text;}
    function setError(text){error.hidden=!text;error.textContent=text;}
    function refreshFields(){for(const [name,input] of Object.entries(fields))input.value=state[name]??'';}
    function animate(element,frames){if(motion()&&element.animate)element.animate(frames,{duration:180,easing:'ease-out'});}
    function clearDrag(){dragging=null;root.classList.remove('is-dragging');root.querySelectorAll('.is-drag-source,.is-drop-target').forEach(e=>e.classList.remove('is-drag-source','is-drop-target'));}
    function dragSource(element,umi){
      element.draggable=true;
      element.addEventListener('dragstart',event=>{dragging=key(umi);selected=dragging;event.dataTransfer.effectAllowed='move';event.dataTransfer.setData('text/plain',dragging);root.classList.add('is-dragging');element.classList.add('is-drag-source');});
      element.addEventListener('dragend',clearDrag);
    }
    function dropTarget(element,position){
      element.addEventListener('dragover',event=>{if(!dragging)return;event.preventDefault();event.dataTransfer.dropEffect='move';element.classList.add('is-drop-target');});
      element.addEventListener('dragleave',event=>{if(!element.contains(event.relatedTarget))element.classList.remove('is-drop-target');});
      element.addEventListener('drop',event=>{if(!dragging)return;event.preventDefault();const umi=state.umis.find(u=>key(u)===dragging);clearDrag();if(umi)moveUMI(umi,position);});
    }
    function moveUMI(umi,position){
      if(umi.position===position)return;
      const other=state.umis.find(u=>u.position===position);
      if(other)other.position=umi.position;
      umi.position=position;selected=key(umi);setError('');render(true);
      announce(other?'UMI positions exchanged.':'UMI moved between '+labels[positions.indexOf(position)]+'.');
    }
    function selectUMI(umi,focus=false){
      selected=key(umi);
      for(const [id,row] of rows)row.root.classList.toggle('is-selected',id===selected);
      track.querySelectorAll('.umi').forEach(chip=>{const active=chip.dataset.umiKey===selected;chip.classList.toggle('selected',active);chip.setAttribute('aria-pressed',String(active));});
      if(focus){const row=rows.get(selected);row?.length.focus();row?.length.select();}
    }
    function addUMI(position){
      state.mode='custom';let umi=state.umis.find(u=>u.position===position);
      if(!umi){umi={position,length:8};state.umis.push(umi);}
      selected=key(umi);setError('');render(true);selectUMI(umi,true);announce('Added an 8 bp UMI between '+labels[positions.indexOf(position)]+'.');
    }
    function removeUMI(umi){
      const ordered=sorted(),index=ordered.indexOf(umi);state.umis=state.umis.filter(u=>u!==umi);
      const next=state.umis.length?sorted()[Math.min(index,state.umis.length-1)]:null;
      selected=next?key(next):null;setError('');render(true);
      if(next)rows.get(selected).remove.focus();else add.focus();announce('UMI removed.');
    }
    function makeRow(umi){
      const id=key(umi),title=node('strong',{class:'structure-umi-title'});
      const grip=node('span',{class:'structure-drag-handle',title:'Drag to move this UMI','aria-hidden':'true'},'⠿');dragSource(grip,umi);
      const position=select(positions.map((p,i)=>[p,labels[i]]),umi.position);
      const length=node('input',{type:'number',min:'1',max:'64',step:'1'});
      const remove=node('button',{type:'button',class:'structure-remove'},'Remove');
      const row=node('div',{class:'structure-umi-row','data-umi-key':id},node('div',{class:'structure-umi-name'},grip,title),node('label',{class:'structure-position'},node('span',{},'Position'),position),node('label',{class:'structure-length'},node('span',{},'Length'),node('span',{},length,node('span',{class:'structure-bp'},'bp'))),remove);
      row.addEventListener('focusin',()=>selectUMI(umi));
      position.addEventListener('change',()=>moveUMI(umi,position.value));
      length.addEventListener('input',()=>{umi.length=length.value===''?null:Number(length.value);const chip=track.querySelector('.umi[data-umi-key="'+id+'"] small');if(chip)chip.textContent=umi.length===null?'Set length':umi.length+' bp';setError('');});
      remove.addEventListener('click',()=>removeUMI(umi));
      return {root:row,title,position,length,remove};
    }
    function renderRows(withMotion){
      const current=sorted(),valid=new Set(current.map(key));
      for(const [id,row] of rows)if(!valid.has(id)){row.root.remove();rows.delete(id);}
      current.forEach((umi,index)=>{
        const id=key(umi);let row=rows.get(id);const fresh=!row;
        if(!row){row=makeRow(umi);rows.set(id,row);}
        row.root.dataset.position=umi.position;row.root.classList.toggle('is-selected',id===selected);
        row.title.textContent='UMI '+(index+1);row.position.value=umi.position;
        if(document.activeElement!==row.length)row.length.value=umi.length??'';
        row.position.setAttribute('aria-label','UMI '+(index+1)+' position');row.length.setAttribute('aria-label','UMI '+(index+1)+' length in base pairs');row.remove.setAttribute('aria-label','Remove UMI '+(index+1));
        for(const option of row.position.options){const other=current.find(u=>u!==umi&&u.position===option.value);option.textContent=labels[positions.indexOf(option.value)]+(other?' · swap with UMI '+(current.indexOf(other)+1):'');}
        for(const input of [row.position,row.length,row.remove])input.disabled=state.mode!=='custom';
        // Keep existing inputs attached while typing; move only when order changes.
        if(list.children[index]!==row.root)list.insertBefore(row.root,list.children[index]||null);
        if(withMotion&&fresh)animate(row.root,[{opacity:0,transform:'translateY(-5px)'},{opacity:1,transform:'translateY(0)'}]);
      });
      umiTitle.textContent='UMIs · '+current.length;umiPanel.hidden=state.mode!=='custom'||!current.length;
    }
    function render(withMotion=false){
      const before=new Map();if(withMotion)root.querySelectorAll('[data-umi-key]').forEach(el=>before.set(el.classList.contains('umi')?'chip-'+el.dataset.umiKey:'row-'+el.dataset.umiKey,el.getBoundingClientRect()));
      const customizing=state.mode==='custom',ordered=sorted();notice.hidden=true;
      root.classList.toggle('is-custom',customizing);auto.setAttribute('aria-pressed',String(!customizing));custom.setAttribute('aria-pressed',String(customizing));
      description.textContent=customizing?'Customize adapter sequences or add UMIs with fixed lengths.':'Detect P5, P7 and i5/i7 indexes automatically during analysis.';
      hint.textContent=customizing?'Click an Adapter to edit its sequence. Click + to add a UMI; drag a UMI to move or swap it.':'Short arrows mark the NGS read starts and directions, not the read lengths.';
      track.replaceChildren(node('span',{class:'structure-end','aria-hidden':'true'},'5′'));
      blocks.forEach(([label,kind,side],index)=>{
        const adapterField=index===2?'adapter_i5_sequence':'adapter_i7_sequence';
        const editable=customizing&&(['P5','P7','i5','i7'].includes(label)||kind==='anchor');
        const block=node(editable?'button':'span',{class:'structure-block '+kind,...(editable?{type:'button','aria-label':'Edit '+(kind==='anchor'?side+' Adapter':label)+' settings'}:{})},node('strong',{},label));
        if(kind==='index')block.append(node('small',{},customizing&&state[label+'_length']!==null?state[label+'_length']+' bp':'Auto'));
        if(kind==='handle')block.append(node('small',{},customizing&&state[label.toLowerCase()+'_sequence']?'Custom':'Auto'));
        if(kind==='insert')block.append(node('small',{},'DNA'));
        if(kind==='anchor'){block.append(node('small',{},side));if(customizing&&state[adapterField])block.append(node('small',{},'Custom'));}
        if(editable)block.addEventListener('click',()=>{advanced.open=true;fields[kind==='anchor'?adapterField:kind==='index'?label+'_length':label.toLowerCase()+'_sequence'].focus();});
        track.append(block);if(index>=positions.length)return;
        const position=positions[index],umi=state.umis.find(u=>u.position===position);
        if(customizing&&umi){
          const id=key(umi),number=ordered.indexOf(umi)+1;
          const chip=node('button',{type:'button',class:'structure-block umi'+(selected===id?' selected':''),'aria-pressed':String(selected===id),'aria-label':'Edit UMI '+number+' between '+labels[index]+', '+umi.length+' bp',title:'Drag to move or swap; click to edit','data-position':position,'data-umi-key':id},node('strong',{},'UMI '+number),node('small',{},umi.length===null?'Set length':umi.length+' bp'));
          chip.addEventListener('click',()=>selectUMI(umi,true));dragSource(chip,umi);dropTarget(chip,position);track.append(chip);
        }else if(customizing){const gap=node('button',{type:'button',class:'structure-gap','aria-label':'Add UMI between '+labels[index],title:'Add or drop a UMI here','data-position':position,onclick:()=>addUMI(position)},'+');dropTarget(gap,position);track.append(gap);}
        else track.append(node('span',{class:'structure-connector','aria-hidden':'true'},'—'));
      });
      track.append(node('span',{class:'structure-end','aria-hidden':'true'},'3′'));
      [...track.children].forEach((part,index)=>{part.style.gridColumn=String(index+1);});
      const left=customizing&&state.umis.some(u=>u.position==='before_insert'),right=customizing&&state.umis.some(u=>u.position==='after_insert');
      const r1=left?'Read 1 starts at the i5-side UMI':'Read 1 starts at the left end of the insert';
      const r2=right?'Read 2 starts at the i7-side UMI':'Read 2 starts at the right end of the insert on the opposite strand';
      const directions=node('div',{class:'structure-directions',role:'group','aria-label':'NGS sequencing directions'},node('span',{class:'structure-read forward','aria-label':r1,title:r1},node('span',{},'R1'),node('span',{class:'structure-read-arrow','aria-hidden':'true'})),node('span',{class:'structure-read reverse','aria-label':r2,title:r2},node('span',{},'R2'),node('span',{class:'structure-read-arrow','aria-hidden':'true'})));
      directions.style.gridColumn=(left?7:8)+' / '+(right?10:9);track.append(directions);
      advanced.hidden=add.hidden=!customizing;add.disabled=state.umis.length===positions.length;
      for(const input of Object.values(fields))input.disabled=!customizing;
      renderRows(withMotion);
      if(withMotion)root.querySelectorAll('[data-umi-key]').forEach(el=>{const id=(el.classList.contains('umi')?'chip-':'row-')+el.dataset.umiKey,old=before.get(id),now=el.getBoundingClientRect();if(old&&old.width&&now.width)animate(el,[{transform:'translate('+(old.left-now.left)+'px,'+(old.top-now.top)+'px)'},{transform:'translate(0,0)'}]);else if(el.classList.contains('umi'))animate(el,[{opacity:0,transform:'scale(.88)'},{opacity:1,transform:'scale(1)'}]);});
    }
    auto.addEventListener('click',()=>{state.mode='auto';clearDrag();setError('');render(true);announce('Auto detection selected. Custom settings are kept as a draft.');});
    custom.addEventListener('click',()=>{state.mode='custom';setError('');render(true);});
    add.addEventListener('click',()=>{const p=['before_insert','after_insert',...positions].find(p=>!state.umis.some(u=>u.position===p));if(p)addUMI(p);});
    for(const [name,input] of Object.entries(fields))input.addEventListener('input',()=>{state[name]=name.endsWith('_length')?(input.value===''?null:Number(input.value)):name.endsWith('_sequence')?(input.value.replace(/\s+/g,'').toUpperCase()||null):input.value;if(name.endsWith('_sequence'))input.value=state[name]||'';setError('');render();});
    function value(){
      if(state.mode==='auto')return {mode:'auto'};
      try{
        for(const umi of state.umis)if(!Number.isInteger(umi.length)||umi.length<1||umi.length>64)throw Error('Set each UMI length to a whole number from 1 to 64 bp.');
        for(const name of ['i5_length','i7_length'])if(state[name]!==null&&(!Number.isInteger(state[name])||state[name]<0||state[name]>24))throw Error('Index lengths must be Auto or a whole number from 0 to 24.');
        for(const name of ['p5_sequence','p7_sequence','adapter_i5_sequence','adapter_i7_sequence'])if(state[name]!==null&&!/^[ACGT]{16,80}$/.test(state[name]))throw Error('Custom handle and adapter sequences need 16–80 A/C/G/T bases.');
        for(const input of Object.values(fields))if(!input.checkValidity())throw Error('Check the adapter and index settings.');
        const result=structuredClone(state);result.umis.sort((a,b)=>positions.indexOf(a.position)-positions.indexOf(b.position));setError('');return result;
      }catch(e){setError(e.message);throw e;}
    }
    apply.addEventListener('click',()=>{try{onApply(value());}catch(e){setError(e.message);}});
    refreshFields();render();
    return {root,value,applyButton:apply,setValue(value){clearDrag();state={...defaults(),...structuredClone(value)};selected=state.umis[0]?key(state.umis[0]):null;refreshFields();setError('');render();},feedback(text){notice.hidden=false;notice.textContent=text;}};
  }
  return {create};
})();
