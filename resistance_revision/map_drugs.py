"""Combine database responses into conservative, reviewable drug mappings."""
from normalize_web import *

def ckey(s):
    s=key(s)
    # Equivalent spacing/hyphens ONLY in letter-prefix numeric development codes.
    if re.fullmatch(r'[a-z]{1,10}[ -]?\d{2,}',s):s=re.sub(r'[ -]','',s)
    return s

def main():
    bulk=CACHE/'pubchem_bulk';molecules={};alias=collections.defaultdict(set);direct={}
    for f in sorted(bulk.glob('synonyms_*.json')):
        d=json.load(open(f))
        for item in d.get('response',{}).get('InformationList',{}).get('Information',[]):
            cid=int(item['CID']);molecules.setdefault(cid,{'CID':cid})['synonyms']=item.get('Synonym',[])
    for f in sorted(bulk.glob('properties_*.json')):
        d=json.load(open(f))
        for item in d.get('response',{}).get('PropertyTable',{}).get('Properties',[]):molecules.setdefault(int(item['CID']),{}).update(item)
    cached=json.load(open(CACHE/'pubchem_names.json'))
    for term,d in cached.items():
        items=d.get('response',{}).get('InformationList',{}).get('Information',[])
        direct[ckey(term)]=set(int(x['CID']) for x in items)
        for item in items:
            cid=int(item['CID']);molecules.setdefault(cid,{'CID':cid}).setdefault('synonyms',item.get('Synonym',[]))
    for cid,info in molecules.items():
        for s in info.get('synonyms',[]):alias[ckey(s)].add(cid)
        if 'Title' in info:alias[ckey(info['Title'])].add(cid)
    def select(name):
        if classify(name)!='agent':return [],'',[]
        q=ckey(name)
        found=set(alias.get(q,[]))
        direct_ids=direct.get(q,set())
        if len(direct_ids)==1:return sorted(direct_ids),'PubChem name service: unique candidate',sorted(found-direct_ids)
        titles=[c for c in found if ckey(molecules[c].get('Title',''))==q]
        if len(titles)==1:return titles,'exact preferred compound title',sorted(found-set(titles))
        if len(found)==1:return sorted(found),'exact compound synonym',[]
        return [],'',sorted(found)
    def cname(cid):
        info=molecules[cid];ss=[s for s in info.get('synonyms',[]) if not re.fullmatch(r'\d{2,7}-\d{2}-\d|\d+',s)]
        # PubChem's leading synonym is often a readable generic name or code;
        # its Title can instead be a long systematic chemical name.
        title=info.get('Title','')
        return ss[0] if ss and (not title or len(ss[0])<len(title) or re.fullmatch(r'[A-Za-z]+-?\d+',ss[0])) else title or (ss[0] if ss else 'PubChem CID '+str(cid))
    rows=[];labels=json.load(open(CACHE/'scoped_labels.json'))['Drug/Agent']
    summary=json.load(open(CACHE/'pubchem_search_summary.json'));failed=set(summary['failed_batches'])
    failed_terms=set()
    for i in failed:
        d=json.load(open(bulk/f'search_{i:03d}.json'));failed_terms.update(ckey(x) for x in d['query_labels'])
    for v in labels:
        label=v['label'];typ=classify(label);chosen=[];candidates=[];normal=label;method='';note='';status='Unmapped'
        if typ!='agent' and len(pieces(label))==1:
            status='Retain original';method='non-single-drug / unresolved study code classification'
            note='Keep the intervention or class name; do not assign a single compound ID.'
        else:
            chosen,method,candidates=select(label)
            if chosen:status='Matched';normal=cname(chosen[0])
            else:
                ps=pieces(label);groups=[];methods=[];extras=[]
                if len(ps)>1:
                    for p in ps:
                        ids,m,cs=select(p);groups.append(ids);methods.append(m);extras+=cs
                    if all(len(x)==1 for x in groups):
                        chosen=list(dict.fromkeys(c for group in groups for c in group));candidates=list(dict.fromkeys(extras))
                        if len(chosen)==1:status='Alias label';normal=cname(chosen[0]);method='all components resolve to same CID'
                        else:
                            normal=' / '.join(cname(c) for c in chosen);method='component-wise compound mapping'
                            status='Combination' if re.search(r'\s\+\s|\b(?:plus|combined with|in combination with)\b',label,re.I) else 'Multiple agents'
                            note='Original agent grouping retained; slash/and does not alone demonstrate co-treatment.'
                    else:
                        known=list(dict.fromkeys(c for group in groups for c in group))
                        status='Partial mapping' if known else 'Needs review'
                        candidates=list(dict.fromkeys([*known,*extras,*candidates]));method='component-wise lookup, incomplete or ambiguous'
                        detail=[]
                        for p,ids in zip(ps,groups):detail.append(p+(' → CID:'+','.join(map(str,ids)) if ids else ' → unresolved'))
                        note='; '.join(detail)+'. No whole-label single CID assigned.'
                else:
                    # Parenthetical aliases are accepted only when full name fails and
                    # two independent parts converge on exactly the same CID.
                    m=re.fullmatch(r'(.+?)\s+\(([^()]+)\)',label)
                    if m:
                        a,ma,ca=select(m[1]);b,mb,cb=select(m[2])
                        if len(a)==len(b)==1 and a==b:
                            chosen=a;normal=cname(a[0]);status='Alias label';method='parenthetical aliases converge on same CID';candidates=list(set(ca+cb))
                        else:candidates=list(dict.fromkeys(candidates+a+b+ca+cb))
                    if not chosen:
                        status='Ambiguous' if candidates else 'Unmapped';method='multiple compound candidates' if candidates else 'no verified exact compound name/synonym'
                        note='No automatic salt, stereoisomer or formulation collapse.' if candidates else 'Keep original label; study-specific agents and biologics may lack a PubChem compound identifier.'
            if not chosen and ckey(label) in failed_terms:
                note+=' Exact-name batch retrieval failed; retry required.'
        if chosen and candidates and not note:note='Preferred exact/name-service mapping selected; other database synonym candidates retained for review.'
        if chosen and typ=='agent' and re.search(r'\b(?:hydrochloride|mesylate|maleate|malate|sulfate|sodium|chloride)\b',label,re.I):
            note+=' Salt/form explicitly preserved in original label; no parent conversion.'
        rows.append({'Original drug/agent':label,'Normalized drug/agent':normal,'Drug ID':' | '.join('CID:'+str(c) for c in chosen) or '-','Status':status,'Web rows':v['row_count'],'Literature rows':v['literature_rows'],'Omics rows':v['omics_rows'],'Agent type':typ,'Method':method,'Candidates':' | '.join('CID:'+str(c)+' '+cname(c) for c in candidates),'Notes':note.strip(),'Example source':v['sources'][0],'Database URL':' | '.join('https://pubchem.ncbi.nlm.nih.gov/compound/'+str(c) for c in chosen or candidates),'InChIKey':' | '.join(molecules[c].get('InChIKey','-') for c in chosen) or '-','SMILES':' | '.join(molecules[c].get('SMILES',molecules[c].get('IsomericSMILES','-')) for c in chosen) or '-'})
    # Product-label identity takes precedence over a coincidental compound title.
    # Parent drug grouping is explicit, with the malate form retained in Notes.
    sutent_source='https://dailymed.nlm.nih.gov/dailymed/drugInfo.cfm?setid=43a4d7f8-48ae-4a63-9108-2fa8e3ea9d9c'
    for r in rows:
        if key(r['Original drug/agent']) in ['sutent','sutent / sunitinib']:
            r.update({'Normalized drug/agent':'Sunitinib','Drug ID':'CID:5329102','Status':'Active ingredient','Method':'official product label + explicit active-moiety grouping','Candidates':'CID:6456015 Sunitinib Malate (product salt form)','Notes':'Sutent contains sunitinib malate (CID:6456015). Drug-name grouping uses the sunitinib active moiety (CID:5329102); the original formulation label is retained. CID:91826462 is not used.','Database URL':'https://pubchem.ncbi.nlm.nih.gov/compound/5329102 | '+sutent_source,'InChIKey':molecules[5329102].get('InChIKey','-'),'SMILES':molecules[5329102].get('SMILES','-')})
        if r['Original drug/agent']=='14H / MK2206':
            r['Notes']+=' 14H is a study-specific SGK inhibitor code; no structure is assigned from the unrelated same-name PubChem hit. See PMID:31461270.'
            r['Database URL']+=' | https://pubmed.ncbi.nlm.nih.gov/31461270/'
    save(CACHE/'pubchem_molecules.json',molecules);save(CACHE/'drug_mapping.json',rows);write_tsv(CACHE/'drug_mapping.tsv',rows)
    print('Drug statuses:',collections.Counter(r['Status'] for r in rows))
    cc=collections.Counter()
    for r in rows:cc[r['Status']]+=r['Web rows']
    print('Drug row coverage:',cc)
if __name__=='__main__':main()
