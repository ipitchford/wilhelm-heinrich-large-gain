import sys, urllib.request, urllib.parse, json, time
def inv(a):
    if not a: return ''
    pos={}
    for w,ps in a.items():
        for p in ps: pos[p]=w
    return ' '.join(pos[i] for i in sorted(pos))
def work(doi_or_id, cites=True, n=60):
    key = doi_or_id if doi_or_id.startswith('W') else 'https://doi.org/'+doi_or_id
    try:
        d=json.load(urllib.request.urlopen('https://api.openalex.org/works/'+key,timeout=60))
    except Exception as e:
        print('ERR',doi_or_id,e); return
    print('##',d['id'],d['title'],d['publication_year'],'cited_by',d['cited_by_count'])
    print('ABSTRACT:',inv(d.get('abstract_inverted_index'))[:2500])
    oa=d.get('open_access',{}); print('OA:',oa.get('oa_url'))
    if cites:
        w=d['id'].split('/')[-1]
        c=json.load(urllib.request.urlopen('https://api.openalex.org/works?filter=cites:%s&per-page=%d&select=id,doi,title,publication_year,authorships,primary_location&sort=publication_year:desc'%(w,n),timeout=60))
        print('-- cited by (%d):'%c['meta']['count'])
        for r in c['results']:
            au=', '.join(a['author']['display_name'] for a in r['authorships'][:3])
            src=((r.get('primary_location') or {}).get('source') or {}).get('display_name') or ''
            print('  ',r['publication_year'],'|',(r['title'] or '')[:120],'|',au,'|',src[:30],'|',r.get('doi'))
if __name__=='__main__':
    for a in sys.argv[1:]:
        work(a); time.sleep(1)
