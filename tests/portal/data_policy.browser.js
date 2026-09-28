(async()=>{
  const checks={};
  const check=(k,v)=>{checks[k]=Boolean(v);if(!v)throw Error(k);};
  const change=e=>e.dispatchEvent(new Event('change',{bubbles:true}));
  const pause=()=>new Promise(r=>setTimeout(r,25));
  try{
    for(let i=0;i<200&&!catalog;i++)await pause();
    const acknowledgment=$('data-policy');
    check('unchecked_and_required',!acknowledgment.checked&&acknowledgment.required);
    addFiles([new File(['@test\nACGT\n+\nIIII\n'],'example.fastq'),new File(['@test2\nTGCA\n+\nIIII\n'],'example2.fastq')]);
    for(const entry of files){entry.ref.value='dm6';change(entry.ref);entry.assay.value='chromatin';change(entry.assay);}
    $('email').value='test@example.org';$('consent').checked=true;
    let calls=0,payload;const original=api;
    api=async(url,method,body)=>{calls++;payload=body;throw Error('Synthetic test stops before upload');};
    try{
      check('start_disabled_without_acknowledgment',$('submit-button').disabled);
      $('submission-form').dispatchEvent(new Event('submit',{bubbles:true,cancelable:true}));await pause();
      check('manual_submit_cannot_bypass_acknowledgment',calls===0&&$('form-error').textContent.includes('human patient data')&&document.activeElement===acknowledgment);
      acknowledgment.checked=true;change(acknowledgment);
      check('enabled_after_active_confirmation',!$('submit-button').disabled);
      $('submission-form').dispatchEvent(new Event('submit',{bubbles:true,cancelable:true}));await pause();
      check('acknowledgment_in_every_file_payload',calls===1&&payload.files.length===2&&payload.files.every(f=>f.data_policy_acknowledged===true));
      acknowledgment.checked=false;change(acknowledgment);
      check('unchecking_blocks_further_submission',$('submit-button').disabled);
      $('submission-form').dispatchEvent(new Event('submit',{bubbles:true,cancelable:true}));await pause();
      check('no_request_after_unchecking',calls===1);
      acknowledgment.checked=true;change(acknowledgment);
      window.dispatchEvent(new PageTransitionEvent('pageshow',{persisted:true}));
      check('history_restore_requires_new_confirmation',!acknowledgment.checked&&$('submit-button').disabled);
      check('no_horizontal_overflow',document.documentElement.scrollWidth<=innerWidth);
    }finally{api=original;}
    $('form-error').hidden=true;document.querySelector('.submission-consents').scrollIntoView({block:'center',behavior:'instant'});
    window.uiVerification={ok:true,checks,viewport:innerWidth};
  }catch(e){window.uiVerification={ok:false,checks,error:String(e.stack||e)};}
})();
