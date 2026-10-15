(async () => {
  const checks = {};
  const check = (name, value) => {
    checks[name] = Boolean(value);
    if (!value) throw Error(name);
  };
  const change = (field, value) => {
    field.value = value;
    field.dispatchEvent(new Event('change', { bubbles: true }));
  };
  try {
    for (let i = 0; i < 200 && !catalog; i++) await new Promise(resolve => setTimeout(resolve, 25));
    const note = 'BWA-MEM2 is unavailable. BWA-MEM is available as an alternative.';
    catalog.references.push({ ...catalog.references[0], id: 'classic', label: 'Classic BWA reference',
      mappers: ['bowtie2', 'bwa-mem'], microc_mapper: 'bwa-mem', availability_note: note });
    catalog.references.push({ ...catalog.references[0], id: 'both', label: 'Both BWA variants',
      mappers: ['bowtie2', 'bwa-mem2', 'bwa-mem'] });
    addFiles([new File(['@synthetic\nACGT\n+\nIIII\n'], 'example.fastq')]);
    const entry = files[0];
    check('initial_note_hidden', entry.availabilityNote.hidden && entry.mapper.disabled);
    change(entry.ref, 'both');
    change(entry.assay, 'microc');
    check('mem2_default_with_classic_available', entry.mapper.value === 'auto' &&
      entry.mapper.options[0].textContent === 'Default · BWA-MEM2' &&
      [...entry.mapper.options].some(option => option.value === 'bwa-mem'));
    change(entry.assay, '');
    change(entry.ref, 'classic');
    check('reference_alone_does_not_show_note', entry.availabilityNote.hidden);
    change(entry.assay, 'chromatin');
    check('chromatin_uses_bowtie2_without_microc_note', entry.availabilityNote.hidden &&
      entry.mapper.options[0].textContent === 'Default · Bowtie2');
    change(entry.assay, 'microc');
    check('microc_displays_classic_default_and_note', !entry.availabilityNote.hidden &&
      entry.availabilityNote.textContent === note && entry.mapper.options[0].textContent === 'Default · BWA-MEM');
    check('unavailable_mappers_not_offered', [...entry.mapper.options].map(option => option.value).join(',') === 'auto,bwa-mem');
    change(entry.ref, 'dm6');
    check('switching_reference_restores_mem2_and_hides_note', entry.availabilityNote.hidden &&
      entry.mapper.options[0].textContent === 'Default · BWA-MEM2');
    change(entry.ref, 'classic');
    change(entry.assay, 'chromatin');
    check('switching_back_to_chromatin_hides_note', entry.availabilityNote.hidden);
    change(entry.assay, '');
    check('clearing_assay_disables_mapper', entry.availabilityNote.hidden && entry.mapper.disabled);
    check('no_horizontal_overflow', document.documentElement.scrollWidth <= innerWidth);
    window.uiVerification = { ok: true, checks, viewport: innerWidth };
  } catch (error) {
    window.uiVerification = { ok: false, checks, error: String(error.stack || error) };
  }
})();
