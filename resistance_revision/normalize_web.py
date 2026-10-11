"""Draft-only mappings for the 76,786 website evidence rows. No release edits."""
import csv, json, re, hashlib, unicodedata, urllib.request, urllib.parse, urllib.error
import threading, time, concurrent.futures, collections
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
W=Path(__file__).resolve().parent
CACHE=W/'过程文件'; CACHE.mkdir(exist_ok=True)
csv.field_size_limit(10_000_000)
def read_tsv(p):
    with open(p,encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f,delimiter='\t'))
def write_tsv(p,rows,headers=None):
    if headers is None:headers=list(rows[0])
    with open(p,'w',encoding='utf-8-sig',newline='') as f:
        writer=csv.DictWriter(f,headers,delimiter='\t');writer.writeheader();writer.writerows(rows)
def save(p,x):p.write_text(json.dumps(x,ensure_ascii=False,indent=2),encoding='utf-8')
def key(s):
    s=unicodedata.normalize('NFKC',s).replace('–','-').replace('—','-').replace('’',"'")
    return re.sub(r'\s+',' ',s).strip().casefold()
def sha(p):
    h=hashlib.sha256()
    with open(p,'rb') as f:
        for b in iter(lambda:f.read(1048576),b''):h.update(b)
    return h.hexdigest()

def scope():
    src=ROOT/'input'
    rr=read_tsv(src/'evidence.tsv');pp=read_tsv(src/'provenance.tsv')
    assert len(rr)==len(pp)==119713
    rows=[]
    for r,p in zip(rr,pp):
        if r['Evidence type']=='Phosphoproteomics evidence' or r['Manual review']=='yes' and p['source_id']!='LIT_035077':
            rows.append({**p,**r})
    assert len(rows)==76786
    assert collections.Counter(r['Evidence type'] for r in rows)=={'Literature evidence':23757,'Phosphoproteomics evidence':53029}
    save(CACHE/'website_rows.json',rows)
    protected=list((ROOT/'protected').glob('*'))
    if not (CACHE/'protected_hashes.json').exists():save(CACHE/'protected_hashes.json',{str(p):sha(p) for p in protected})
    labels={}
    for field in ['Drug/Agent','Disease']:
        d={}
        for r in rows:
            label=r[field]; v=d.setdefault(label,{'label':label,'row_count':0,'literature_rows':0,'omics_rows':0,'sources':set(),'contexts':[]})
            v['row_count']+=1;v['literature_rows' if r['Evidence type']=='Literature evidence' else 'omics_rows']+=1
            v['sources'].add(r['Source'])
            if len(v['contexts'])<3:v['contexts'].append({'source':r['Source'],'detail':r['Details'],'drug':r['Drug/Agent'],'disease':r['Disease']})
        labels[field]=sorted([{**v,'sources':sorted(v['sources'])} for v in d.values()],key=lambda v:(-v['row_count'],v['label']))
    save(CACHE/'scoped_labels.json',labels)
    print('Scope complete:',len(rows),{k:len(v) for k,v in labels.items()},flush=True)
    return labels

RATE_LOCK=threading.Lock(); NEXT=0
def get_json(url):
    global NEXT
    for attempt in range(3):
        with RATE_LOCK:
            delay=max(0,NEXT-time.monotonic());NEXT=max(NEXT,time.monotonic())+.6
        if delay:time.sleep(delay)
        try:
            req=urllib.request.Request(url,headers={'User-Agent':'KinaseMD-normalization-draft/1.0','Accept':'application/json'})
            with urllib.request.urlopen(req,timeout=25) as f:return json.load(f)
        except urllib.error.HTTPError as e:
            if e.code==404:return {'not_found':True}
            if e.code not in [429,500,502,503,504]:return {'error':str(e)}
            retry=min(60,max(5,int(e.headers.get('Retry-After','10'))))
            with RATE_LOCK:NEXT=max(NEXT,time.monotonic()+retry)
        except Exception as e:
            if attempt==2:return {'error':str(e)}
            time.sleep(1+attempt)
    return {'error':'retry exhausted'}

def classify(label):
    s=key(label)
    if s in ['','-','n/a','na','not specified','none','unknown']:return 'missing'
    if re.search(r'\b(?:siRNA|shRNA|miRNA|microRNA|knockdown|knockout|overexpression|transfection|CRISPR|silencing)\b',label,re.I):return 'genetic intervention'
    if re.search(r'\b(?:irradiation|radiation|radiotherapy|hypoxia|heat shock|exercise|UVB|UVA|ultraviolet)\b',label,re.I):return 'physical intervention'
    if re.search(r'\b(?:extract|decoction|formula|conditioned medium|serum|lysate|essential oil)\b',label,re.I):return 'mixture or preparation'
    if re.search(r'^(?:compound|compund|compounds|complex)\s+[A-Za-z0-9-]+(?:\s|$)',s) or re.fullmatch(r'[A-Z]\d{0,2}|\d{1,2}[a-z]{0,2}|\d{2}[A-Z]|\d{1,2}[a-z]-[A-Z]\d',label):return 'study-specific compound'
    if re.search(r'\b(?:inhibitors|activators|agonists|antagonists|statins|NSAIDs|antioxidants)\b',label,re.I):return 'agent class'
    if re.fullmatch(r'(?:PI3K|AKT\d?|ERK\d?|MEK\d?|CDK\d+|EGFR|mTOR|PKA|PKC|JNK|p38|STAT\d|SGLT2|GPER1)\s+(?:inhibitor|activator|agonist|antagonist)',label,re.I):return 'agent class'
    if re.search(r'\b(?:inhibitor|activator|agonist|antagonist)\b',label,re.I) and not re.search(r'\d',label) and len(label.split())<=5 and not re.search(r'[/+;]|\band\b',label):return 'agent class'
    return 'agent'

def pieces(label):
    # Slash means "multiple or alias", never automatically "combination".
    # Only spaced plus signs; a stereochemical (+)- prefix is not a delimiter.
    parts=re.split(r'\s+(?:\+|/|and|plus|combined with|in combination with|with)\s+|;\s*',label,flags=re.I)
    return [p.strip() for p in parts if p.strip()]

def drug_queries(labels):
    cp=CACHE/'pubchem_names.json'
    cache=json.load(open(cp)) if cp.exists() else {}
    old=json.load(open(ROOT/'references/pubchem_mapping_cache.json'))
    for term,data in old.items():
        cache.setdefault(key(term),{'query':term,'response':data,'origin':'local PubChem snapshot 2026-09-17','url':'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/'+urllib.parse.quote(term,safe='')+'/synonyms/JSON'})
    queries=collections.OrderedDict()
    for v in labels['Drug/Agent']:
        label=v['label']
        if classify(label)!='agent':continue
        queries.setdefault(key(label),label)
        for p in pieces(label):
            if classify(p)=='agent':queries.setdefault(key(p),p)
        # Parenthetical generic/code names become candidates only; preserve full label.
        m=re.fullmatch(r'(.+?)\s+\(([^()]+)\)',label)
        if m:
            for p in m.groups():
                if classify(p)=='agent' and not re.fullmatch(r'\d+\s*(?:nM|uM|µM|mg|μM).*',p,re.I):queries.setdefault(key(p),p)
    # Exact synonyms already seen can avoid network requests; preserve candidate multiplicity.
    synonym_index=collections.defaultdict(set)
    for entry in cache.values():
        for info in entry.get('response',{}).get('InformationList',{}).get('Information',[]):
            for s in info.get('Synonym',[]):synonym_index[key(s)].add(info['CID'])
    pending=[(k,v) for k,v in queries.items() if k not in cache and k not in synonym_index]
    print('PubChem queries:',len(queries),'uncached:',len(pending),flush=True)
    def worker(item):
        k,label=item;url='https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/'+urllib.parse.quote(label,safe='')+'/synonyms/JSON'
        return k,{'query':label,'url':url,'origin':'PubChem live 2026-10-06','response':get_json(url)}
    done=0
    with concurrent.futures.ThreadPoolExecutor(max_workers=8) as pool:
        for k,entry in pool.map(worker,pending):
            cache[k]=entry;done+=1
            if done%100==0:
                save(cp,cache)
                found=sum(bool(e.get('response',{}).get('InformationList')) for e in cache.values())
                print('PubChem progress:',done,'/',len(pending),'successful names:',found,flush=True)
    save(cp,cache);print('PubChem name scan complete',flush=True)

if __name__=='__main__':
    labels=scope()
    drug_queries(labels)
