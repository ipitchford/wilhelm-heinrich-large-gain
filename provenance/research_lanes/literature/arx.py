import sys, urllib.request, urllib.parse, re, time
import xml.etree.ElementTree as ET
def q(query, n=25):
    url='http://export.arxiv.org/api/query?'+urllib.parse.urlencode({'search_query':query,'max_results':n,'sortBy':'submittedDate','sortOrder':'descending'})
    for attempt in range(3):
        try:
            data=urllib.request.urlopen(url,timeout=60).read()
            break
        except Exception as e:
            print('ERR',e); time.sleep(5); data=None
    if not data: return
    ns={'a':'http://www.w3.org/2005/Atom'}
    root=ET.fromstring(data)
    for e in root.findall('a:entry',ns):
        t=' '.join(e.find('a:title',ns).text.split())
        i=e.find('a:id',ns).text.split('/abs/')[-1]
        d=e.find('a:published',ns).text[:10]
        au=', '.join(a.find('a:name',ns).text for a in e.findall('a:author',ns)[:3])
        print(d,'|',i,'|',t[:130],'|',au)
if __name__=='__main__':
    for qq in sys.argv[1:]:
        print('##',qq); q(qq); time.sleep(3)
