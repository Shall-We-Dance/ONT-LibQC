"""Sublibrary highlighting works offline, with keyboard/touch and in print."""
from pathlib import Path
import sys,tempfile,os
from playwright.sync_api import sync_playwright
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'workflow/scripts'))
from report_profiles import render_profiles_card
with tempfile.TemporaryDirectory() as tmp:
 root=Path(tmp)
 for sid,signal in [('lib01',100),('lib02',60)]:
  (root/f's.{sid}.full_library_profiles.tsv').write_text(f'size_bp\tlength_weighted_relative_signal\n0\t0\n100\t{signal}\n2000\t0\n')
 scopes=[{'scope_id':sid} for sid in ('total','lib01','lib02')]
 doc=render_profiles_card('s',scopes,['<section class="profile-scope">Details</section>']*3,root,root,asset_prefix='')
 (root/'index.html').write_text('<html><head><meta name="viewport" content="width=device-width"><meta http-equiv="Content-Security-Policy" content="default-src \'none\'; script-src \'unsafe-inline\'; style-src \'unsafe-inline\'; img-src \'self\' data:"></head><body>'+doc+'</body></html>')
 with sync_playwright() as p:
  for engine in sys.argv[1:] or ['chromium']:
   browser=getattr(p,engine).launch(**({'executable_path':os.environ['ONT_BROWSER_EXECUTABLE']} if os.environ.get('ONT_BROWSER_EXECUTABLE') else {}))
   for mobile in (False,True):
    page=browser.new_page(viewport={'width':390 if mobile else 1440,'height':900},has_touch=mobile,is_mobile=mobile)
    errors=[];page.on('pageerror',lambda e:errors.append(str(e)));page.goto((root/'index.html').as_uri())
    buttons=page.locator('.profile-comparison-legend button');lines=page.locator('.profile-size-chart polyline')
    assert buttons.count()==2 and lines.count()==2
    original=lines.evaluate_all('(es)=>es.map(e=>e.getAttribute("points"))')
    def opacities():
     page.wait_for_timeout(200)
     return lines.evaluate_all('(es)=>es.map(e=>Number(getComputedStyle(e).opacity))')
    if mobile:
     buttons.nth(1).tap();assert opacities()==[.12,1]
     buttons.nth(1).tap();assert opacities()==[1,1]
    else:
     buttons.first.hover();assert opacities()==[1,.12]
     page.mouse.move(0,0);assert opacities()==[1,1]
     page.keyboard.press('Tab');assert buttons.first.evaluate('(e)=>e===document.activeElement')
     assert opacities()==[1,.12]
     page.keyboard.press('Tab');assert opacities()==[.12,1]
     page.keyboard.press('Escape');assert opacities()==[1,1]
    if mobile:
     buttons.first.tap();assert opacities()==[1,.12]
     assert buttons.first.get_attribute('aria-pressed')=='true'
     page.touchscreen.tap(0,0);assert opacities()==[1,1]
     buttons.first.tap()
    else:
     # A mouse click must never pin a series after the pointer leaves.
     buttons.first.click();assert opacities()==[1,.12]
     assert buttons.first.get_attribute('aria-pressed')=='false'
     page.mouse.move(0,0);assert opacities()==[1,1]
     buttons.nth(1).click();page.mouse.move(0,0);assert opacities()==[1,1]
     buttons.first.hover()
    backgrounds=buttons.evaluate_all('(es)=>es.map(e=>getComputedStyle(e).backgroundColor)')
    assert backgrounds[0]!=backgrounds[1]
    assert buttons.first.evaluate('(e)=>e.style.getPropertyValue("--legend-tint")')=='#2a6cbd0f'
    page.emulate_media(media='print');assert opacities()==[1,1]
    page.emulate_media(media='screen')
    # Print layout can move the legend away from a stationary mouse pointer.
    if not mobile: buttons.first.hover()
    assert opacities()==[1,.12]
    assert lines.evaluate_all('(es)=>es.map(e=>e.getAttribute("points"))')==original
    assert page.locator('.comparison-size a').get_attribute('href')=='s.sublibraries.size_profiles.svg'
    assert not errors,errors
    page.close()
   browser.close();print(engine,'desktop, keyboard, touch, offline, print: OK')
