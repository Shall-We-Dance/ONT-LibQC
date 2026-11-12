"""Public counts, private short links, and legacy links on desktop and mobile."""
import json,os,socket,sys,tempfile,threading,time
from pathlib import Path
import uvicorn
from playwright.sync_api import sync_playwright
from ont_portal.app import create_app
from ont_portal.settings import Settings
from ont_portal.store import token_for
from ont_portal.short_links import link_key,key_hash
OUTPUT=Path(sys.argv[1]) if len(sys.argv)>1 else Path('test-results/status-links')
OUTPUT.mkdir(parents=True,exist_ok=True)
results=[]
def wait_tasks(page,count):
 deadline=time.monotonic()+10
 while time.monotonic()<deadline:
  if page.evaluate('tasks.length === '+str(count)+' && tasks.every(t => t.job && !t.fetchError)'):return
  page.wait_for_timeout(50)
 raise AssertionError('Private page did not load')
with tempfile.TemporaryDirectory() as tmp:
 sock=socket.socket();sock.bind(('127.0.0.1',0));origin='http://127.0.0.1:'+str(sock.getsockname()[1])
 cfg=Settings(state=Path(tmp),public_url=origin,allowed_hosts=['127.0.0.1'],secret='browser-fixture-'*5,submission_key='',public_catalog={'references':[]},require_https=False,min_free_bytes=0)
 app=create_app(cfg);store=app.state.store;store.heartbeat();now=time.time();batch='batch-fixture';keys={}
 with store.connect() as db:
  for n,status in enumerate(('running','queued','succeeded','failed','expired'),1):
   jid=str(n)*32;key=link_key(cfg.secret,'job',jid);keys[jid]=key
   db.execute('INSERT INTO jobs (id,created,updated,status,email,filename,expected_bytes,options,client,expires,batch_id,short_link_hash) VALUES (?,?,?,?,?,?,?,?,?,?,?,?)',
    (jid,now-(5-n)*86400,now,status,'private@example.org','private.fastq',20,json.dumps({'sample':'Sample'+str(n)}),'private-address',now+86400,batch if n<=2 else None,key_hash(key) if n!=3 else None))
  db.execute('INSERT INTO batch_links VALUES (?,?)',(batch,key_hash(link_key(cfg.secret,'batch',batch))))
 server=uvicorn.Server(uvicorn.Config(app,log_level='error',access_log=False));thread=threading.Thread(target=server.run,kwargs={'sockets':[sock]},daemon=True);thread.start()
 while not server.started:time.sleep(.02)
 try:
  with sync_playwright() as p:
   for engine in sys.argv[2:] or ['chromium','webkit']:
    browser=getattr(p,engine).launch(**({'executable_path':os.environ['ONT_WEBKIT_EXECUTABLE']} if engine=='webkit' and os.environ.get('ONT_WEBKIT_EXECUTABLE') else {}))
    try:
     for width in (1280,390):
      context=browser.new_context(viewport={'width':width,'height':900});urls=[];errors=[]
      context.on('request',lambda request:urls.append(request.url));page=context.new_page();page.on('pageerror',lambda error:errors.append(str(error)))
      page.goto(origin)
      page.get_by_text('Jobs running: 1',exact=True).wait_for();page.get_by_text('Jobs waiting to run: 1',exact=True).wait_for()
      page.locator('#server-status summary').click()
      assert page.locator('#jobs-total').inner_text()=='5'
      assert page.locator('#jobs-history-chart svg').count()==1
      assert page.locator('#jobs-history-chart .history-point').count()==5
      for text in ('Completed','Failed','Cancelled','Earlier archived jobs'):assert text not in page.locator('#server-status').inner_text()
      assert 'private' not in page.locator('#server-status').inner_text().lower()
      assert page.evaluate('document.documentElement.scrollWidth <= innerWidth')
      page.screenshot(path=str(OUTPUT/f'{engine}-{width}-counts.png'))
      for path,count in [('/j/#'+keys['1'*32],1),('/b/#'+link_key(cfg.secret,'batch',batch),2),('/jobs/'+'3'*32+'#'+token_for(cfg.secret,'3'*32),1),('/#jobs='+','.join(jid+':'+token_for(cfg.secret,jid) for jid in ('1'*32,'2'*32)),2)]:
       page.goto(origin+path)
       deadline=time.monotonic()+10
       while time.monotonic()<deadline:
        if page.evaluate('tasks.length === '+str(count)+' && tasks.every(t => t.job && !t.fetchError)'):break
        page.wait_for_timeout(50)
       else:raise AssertionError('Private page did not load')
       page.wait_for_url('**/'+('j' if count==1 else 'b')+'/#*')
       assert page.evaluate('document.documentElement.scrollWidth <= innerWidth')
       if path.startswith('/j/'):
        store.update('1'*32,stage='Plotting sublibrary gene coverage (2 of 6 complete)')
        page.evaluate('refreshTask(tasks[0])')
        page.get_by_text('Plotting sublibrary gene coverage (2 of 6 complete)',exact=True).wait_for()
        assert page.locator('.cancel-analysis').evaluate("el => getComputedStyle(el).backgroundColor")=='rgb(255, 241, 240)'
        page.screenshot(path=str(OUTPUT/f'{engine}-{width}-progress.png'))
      # Slow private requests must show a result/loading view from the first paint;
      # an unavailable upload catalog must never block reports.
      page.route('**/api/catalog',lambda route: route.abort())
      def slow(route):
       page.wait_for_timeout(650)
       route.continue_()
      page.route('**/api/links/resolve',slow)
      page.route('**/api/jobs/*',slow)
      for path in ('/b/#'+link_key(cfg.secret,'batch',batch), '/#jobs='+','.join(jid+':'+token_for(cfg.secret,jid) for jid in ('1'*32,'2'*32))):
       page.goto(origin+path,wait_until='domcontentloaded')
       assert page.locator('#submission-view').is_hidden()
       assert page.locator('#page-loading').is_visible() or page.locator('#tasks-view').is_visible()
       wait_tasks(page,2)
       page.wait_for_url('**/b/#*')
       assert page.locator('#submission-view').is_hidden()
      page.unroute('**/api/links/resolve',slow);page.unroute('**/api/jobs/*',slow)
      # Conversion failure leaves the legacy link and reports usable.
      page.route('**/api/links/compact',lambda route: route.fulfill(status=503,json={'detail':'Try later'}))
      legacy='/jobs/'+'3'*32+'#'+token_for(cfg.secret,'3'*32)
      page.goto(origin+legacy)
      wait_tasks(page,1)
      assert page.url==origin+legacy
      page.goto(origin+'/b/#invalid')
      page.locator('#global-error').wait_for()
      assert page.locator('#submission-view').is_hidden()
      assert not context.cookies();assert not errors,errors
      for key in [*keys.values(),link_key(cfg.secret,'batch',batch)]:assert all(key not in url for url in urls)
      results.append({'browser':engine,'width':width,'daily_cumulative_counts':True,'no_homepage_flash':True,'short_and_legacy_links':True,'keys_absent_from_http_urls':True})
      context.close()
    finally:browser.close()
 finally:server.should_exit=True;thread.join(timeout=10);sock.close()
(OUTPUT/'validation.json').write_text(json.dumps(results,indent=2)+'\n');print(json.dumps(results,indent=2))
