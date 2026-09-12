(() => {
  const result = {};
  const event = () => new Event('change', { bubbles: true });
  const bedFile = (entry, name = 'synthetic_regions.bed') => {
    const transfer = new DataTransfer();
    transfer.items.add(new File(['chr2L\t100\t101\n'], name, { type: 'text/plain' }));
    entry.bed.files = transfer.files;
    entry.bed.dispatchEvent(event());
  };
  try {
    addFiles([new File(['@synthetic\nACGT\n+\nIIII\n'], 'synthetic.fastq')]);
    const entry = files[0];
    entry.ref.value = 'dm6'; entry.ref.dispatchEvent(event());
    result.optional_controls_disabled = entry.bedBox.hidden && entry.bedConfirm.disabled && !entry.bedConfirm.required;
    entry.bedCheck.input.click();
    result.controls_visible = !entry.bedBox.hidden && entry.bedConfirm.required && !entry.bedMode.disabled;
    result.plot_modes_match_report = JSON.stringify([...entry.bedMode.options].map(o => [o.value, o.text])) === JSON.stringify([['auto','Auto'],['reference_point','Reference point'],['region','Region (scaled)']]);
    bedFile(entry);
    result.filename_display = entry.bedName.textContent === 'synthetic_regions.bed';
    result.confirmation_required = !entry.bedConfirm.checked && !entry.bedConfirm.checkValidity();
    try { optionsFor(entry); result.unconfirmed_submission_blocked = false; }
    catch (e) { result.unconfirmed_submission_blocked = e.message.includes('dm6 coordinates'); }
    entry.bedConfirm.click();
    result.mode_payloads = ['auto', 'reference_point', 'region'].every(mode => {
      entry.bedMode.value = mode;
      const options = optionsFor(entry);
      return options.regions_mode === mode && options.regions_assembly === 'dm6';
    });
    entry.ref.value = 'hg38'; entry.ref.dispatchEvent(event());
    result.reference_change_resets_confirmation = !entry.bedConfirm.checked && entry.bedConfirmText.textContent === 'I confirm this BED uses hg38 coordinates';
    entry.bedConfirm.click(); bedFile(entry, 'replacement.bed');
    result.file_change_resets_confirmation = !entry.bedConfirm.checked;
    entry.bedCheck.input.click();
    entry.bedLabel.value = '<invalid label>';
    result.optional_fields_do_not_block_form = !entry.bedConfirm.required && entry.bedConfirm.disabled && entry.bedLabel.disabled && entry.bedLabel.checkValidity();
    entry.bedLabel.value = 'Uploaded regions';
    entry.bedCheck.input.click();
    entry.assay.value = 'microc'; entry.assay.dispatchEvent(event());
    result.microc_disables_bed = entry.bedCheck.input.disabled && !entry.bedCheck.input.checked && entry.bedBox.hidden && entry.bedConfirm.disabled;
    entry.assay.value = 'chromatin'; entry.assay.dispatchEvent(event());
    entry.bedCheck.input.click();
    entry.ref.value = 'dm6'; entry.ref.dispatchEvent(event());
    entry.bedMode.value = 'reference_point';
    bedFile(entry, 'synthetic_long_region_filename_for_layout_testing_summits.bed');
    entry.bedConfirm.click();
    addFiles([new File(['@second\nACGT\n+\nIIII\n'], 'second.fastq')]);
    const second = files[1];
    second.ref.value = 'hg38'; second.ref.dispatchEvent(event());
    second.bedCheck.input.click(); bedFile(second); second.bedMode.value = 'region'; second.bedConfirm.click();
    result.batch_settings_independent = optionsFor(entry).regions_mode === 'reference_point' && optionsFor(entry).regions_assembly === 'dm6' && optionsFor(second).regions_mode === 'region' && optionsFor(second).regions_assembly === 'hg38';
    second.card.querySelector('.remove').click();
    entry.bedConfirm.checked = false;
    entry.bedBox.scrollIntoView({ block: 'start' }); window.scrollBy(0, -40);
    const picker = entry.bed.closest('.custom-file-picker').getBoundingClientRect();
    const mode = entry.bedMode.getBoundingClientRect();
    const box = entry.bedBox.querySelector('.custom-controls').getBoundingClientRect();
    result.matching_control_heights = Math.abs(picker.height - 44) < 1 && Math.abs(mode.height - 44) < 1;
    result.controls_within_panel = [picker, mode].every(r => r.left >= box.left && r.right <= box.right + 1);
    result.no_document_overflow = document.documentElement.scrollWidth <= innerWidth;
    result.viewport = innerWidth;
    return result;
  } catch (error) { return { ...result, error: String(error.stack || error) }; }
})();
