import sys, urllib.request, urllib.parse, json, time
def s(search, filt='', n=25, field='search'):
    p={field:search,'per-page':n,'select':'id,doi,title,publication_year,authorships,primary_location,cited_by_count'}
    if filt: p['filter']=filt
    url='https://api.openalex.org/works?'+urllib.parse.urlencode(p)
    try:
        d=json.load(urllib.request.urlopen(url,timeout=60))
    except Exception as e:
        print('ERR',e); return
    print('## %s [%s] count=%s'%(search,filt,d['meta']['count']))
    for r in d['results']:
        au=', '.join(a['author']['display_name'] for a in r['authorships'][:3])
        src=((r.get('primary_location') or {}).get('source') or {}).get('display_name') or ''
        print(r['publication_year'],'|',(r['title'] or '')[:120],'|',au,'|',src[:35],'|',r.get('doi'),'| c=',r['cited_by_count'])
if __name__=='__main__':
    a=sys.argv[1:]
    i=0
    while i<len(a):
        s(a[i], a[i+1] if i+1<len(a) else ''); i+=2; time.sleep(1)
