"""Offline similarity panel: controls, downloads, mobile layout and print visibility."""
from pathlib import Path
import os
import sys
import tempfile
from playwright.sync_api import sync_playwright
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'workflow/scripts'))
from report_similarity import write_assets,report_html

with tempfile.TemporaryDirectory() as tmp:
    root=Path(tmp);assets=root/'assets'
    scopes=[dict(scope_id=f'lib{i:02d}') for i in range(1,4)]
    records={s['scope_id']:{str(j):dict(chrom='chr1',start=10000*j*(i+1),end=10000*j*(i+1)+70,mapq=60)
                          for j in range(15)} for i,s in enumerate(scopes)}
    write_assets('example',scopes,records,{'chr1':5000000},['chr1'],20,assets)
    (root/'index.html').write_text('<!doctype html><html><head><meta name="viewport" content="width=device-width,initial-scale=1"><style>body{font-family:Arial;margin:16px}</style></head><body>'+report_html('example',assets)+'</body></html>')
    with sync_playwright() as p:
        browser_name=sys.argv[1] if len(sys.argv)>1 else 'chromium'
        options={'headless':True}
        if os.environ.get('ONT_BROWSER_EXECUTABLE'):options['executable_path']=os.environ['ONT_BROWSER_EXECUTABLE']
        browser=getattr(p,browser_name).launch(**options)
        page=browser.new_page(viewport={'width':1280,'height':900});errors=[]
        page.on('pageerror',lambda e:errors.append(str(e)))
        page.goto((root/'index.html').as_uri(),wait_until='networkidle')
        for method in ('spearman','pearson'):
            page.get_by_label('Similarity correlation method').select_option(method)
            for size in ('10000','100000','1000000'):
                page.get_by_label('Similarity genomic window').select_option(size)
                for policy in ('all','nonzero'):
                    page.get_by_label('Similarity zero-bin policy').select_option(policy)
                    panel=page.locator('.similarity-panel:visible');assert panel.count()==1
                    assert panel.get_attribute('data-method')==method and panel.get_attribute('data-bin')==size
                    assert panel.get_attribute('data-zeros')==policy
                    img=panel.locator('img');img.scroll_into_view_if_needed()
                    page.wait_for_function('(img)=>img.complete&&img.naturalWidth>0',arg=img.element_handle())
                    for link in panel.locator('a').all():assert (root/link.get_attribute('href')).exists()
        page.emulate_media(media='print');assert page.locator('.similarity-panel:visible').count()==1
        assert not page.locator('.similarity-controls').is_visible()
        page.emulate_media(media='screen');page.set_viewport_size({'width':390,'height':844})
        assert page.evaluate('document.documentElement.scrollWidth<=innerWidth+1')
        assert not errors,errors
        browser.close()
print('Similarity controls, downloads, print and mobile checks passed')
