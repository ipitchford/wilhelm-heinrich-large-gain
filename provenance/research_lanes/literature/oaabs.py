import sys, urllib.request, json, time
def inv(a):
    if not a: return ''
    pos={}
    for w,ps in a.items():
        for p in ps: pos[p]=w
    return ' '.join(pos[i] for i in sorted(pos))
for doi in sys.argv[1:]:
    try:
        d=json.load(urllib.request.urlopen('https://api.openalex.org/works/https://doi.org/'+doi,timeout=60))
        print('##',d['title'],'|',d['publication_year'],'|',', '.join(a['author']['display_name'] for a in d['authorships'][:4]),'| cited',d['cited_by_count'],'| OA',d.get('open_access',{}).get('oa_url'))
        print('   ',inv(d.get('abstract_inverted_index'))[:1800])
    except Exception as e:
        print('ERR',doi,e)
    time.sleep(0.5)
