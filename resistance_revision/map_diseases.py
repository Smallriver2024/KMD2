"""Exact ontology matching, conservative context mapping, reviewable suggestions."""
from normalize_web import *
import difflib

REF=ROOT/'references/disease_ontology'
OLD=REF.parent/'output/literature_disease_bridge.tsv'
def dkey(s):
    s=key(s).replace('’',"'").replace("'",'').replace('leukaemia','leukemia').replace('tumour','tumor')
    s=re.sub(r'\s*\(disease\)$','',s)
    plural={'diseases':'disease','disorders':'disorder','cancers':'cancer','carcinomas':'carcinoma','neoplasms':'neoplasm','tumors':'tumor','adenocarcinomas':'adenocarcinoma','lymphomas':'lymphoma','sarcomas':'sarcoma','infections':'infection'}
    for a,b in plural.items():s=re.sub(r'\b'+a+r'\b',b,s)
    return re.sub(r'\s+',' ',re.sub(r'[-_,]+',' ',s)).strip(' .;')

def parse(path,onto):
    term={}; out=[]
    def emit():
        if term.get('id','').startswith(onto+':') and term.get('name') and not term.get('obsolete'):out.append(dict(term))
    version=''
    with open(path,encoding='utf-8') as f:
        for line in f:
            line=line.rstrip('\n')
            if line.startswith('data-version:'):version=line.split(': ',1)[1]
            if line=='[Term]':emit();term={'synonyms':[]}
            elif line.startswith('['):emit();term={}
            elif line.startswith('id: '):term['id']=line[4:]
            elif line.startswith('name: '):term['name']=line[6:]
            elif line=='is_obsolete: true':term['obsolete']=True
            elif line.startswith('synonym: '):
                m=re.match(r'synonym: "((?:[^"\\]|\\.)*)" EXACT\b',line)
                if m:term.setdefault('synonyms',[]).append(m.group(1).replace('\\"','"'))
        emit()
    return out,version

ABBR={'nsclc':'non-small cell lung cancer','sclc':'small cell lung cancer','hcc':'hepatocellular carcinoma','aml':'acute myeloid leukemia','cml':'chronic myeloid leukemia','cll':'chronic lymphocytic leukemia','all':'acute lymphoblastic leukemia','gbm':'glioblastoma','pdac':'pancreatic ductal adenocarcinoma','rcc':'renal cell carcinoma','tnbc':'triple-negative breast cancer','crc':'colorectal cancer','hnscc':'head and neck squamous cell carcinoma','oscc':'oral squamous cell carcinoma','dlbcl':'diffuse large B-cell lymphoma','gist':'gastrointestinal stromal tumor','t2dm':'type 2 diabetes mellitus','t1dm':'type 1 diabetes mellitus','mm':'multiple myeloma','copd':'chronic obstructive pulmonary disease','hcmv':'cytomegalovirus'}
def expand(s):
    for a,b in ABBR.items():s=re.sub(r'\b'+a+r'\b',b,s,flags=re.I)
    return s

def main():
    entries={}; indexes={};versions={}
    for onto,fn in [('MONDO','mondo.obo'),('DOID','doid.obo')]:
        terms,version=parse(REF/fn,onto);versions[onto]={'version':version,'file':str(REF/fn),'sha256':sha(REF/fn),'term_count':len(terms)}
        preferred=collections.defaultdict(set);synonyms=collections.defaultdict(set)
        for t in terms:
            entries[t['id']]=t
            preferred[dkey(t['name'])].add(t['id'])
            for s in t.get('synonyms',[]):synonyms[dkey(s)].add(t['id'])
        indexes[onto]=(preferred,synonyms)
    save(CACHE/'ontology_versions.json',versions)
    def exact(s):
        q=dkey(s)
        for onto in ['MONDO','DOID']:
            a,b=indexes[onto]
            if q in a:return sorted(a[q]),'preferred label'
            if q in b:
                ids=sorted(b[q]);human=[x for x in ids if 'non-human animal' not in entries[x]['name']]
                if human and len(human)<len(ids) and not re.search(r'\b(?:animal|mouse|murine|canine|feline|rat)\b',s,re.I):
                    return human,'exact synonym; human disease term'
                return ids,'exact synonym'
        return [],''
    prior=collections.defaultdict(list)
    for r in read_tsv(OLD):
        for variant in r['raw_variants'].split(' | ')+[r['display_label']]:
            if r['ontology_id'] in entries and r not in prior[dkey(variant)]:prior[dkey(variant)].append(r)
    labels=json.load(open(CACHE/'scoped_labels.json'))['Disease'];out=[]
    for v in labels:
        label=v['label']; ids,method=exact(label);status='Matched';note='';candidates=[]
        if not ids:
            expanded=expand(label);ids,method=exact(expanded)
            if ids:method='abbreviation expansion + '+method
        if len(ids)>1:status='Ambiguous';candidates=ids;ids=[];note='Exact synonym refers to multiple ontology terms; no automatic selection.'
        if not ids and status!='Ambiguous':
            components=re.split(r'\s*(?:/|;|\||\band\b)\s*',label,flags=re.I)
            # Resolve all components; retain the original context and membership.
            component_ids=[];complete=len(components)>1
            for p in components:
                found,_=exact(expand(p))
                if len(found)!=1:complete=False
                else:component_ids+=found
            if complete:
                ids=list(dict.fromkeys(component_ids));method='component-wise exact mapping';status='Multiple diseases' if len(ids)>1 else 'Alias label';note='Original multi-part label retained; no broad cancer-category collapse.'
        if not ids and status!='Ambiguous':
            stripped=re.sub(r'\b(?:human|patients? with)\s+','',label,flags=re.I)
            stripped=re.sub(r'\s+(?:cells?|cell lines?|cell lines? models?|xenografts?|models?)\b.*$','',stripped,flags=re.I)
            found,meth=exact(expand(stripped))
            if len(found)==1 and stripped!=label:
                ids=found;method='remove model descriptors + '+meth;status='Context mapped';note='Model descriptors remain in original label.'
        if not ids and status!='Ambiguous':
            suggestions=prior.get(dkey(label),[])
            candidates=list(dict.fromkeys(r['ontology_id'] for r in suggestions))
            # Link an explicit disease embedded in a model/resistance label to its
            # broader disease concept, without claiming subtype equivalence.
            safe_modifiers=r'\b(?:metastatic|advanced|resistant|resistance|refractory|mutant|positive|negative|driven|associated|induced|deficient|mice|mouse|model|xenograft|non|relapsed|recurrent)\b'
            embedded=[]
            if re.search(safe_modifiers,label,re.I):
                clean=dkey(expand(label))
                for candidate in candidates:
                    term=entries[candidate]
                    for alias in [term['name'],*term.get('synonyms',[])]:
                        ak=dkey(alias)
                        if len(ak)>6 and re.search(r'(?<!\w)'+re.escape(ak)+r'(?!\w)',clean):
                            embedded.append(candidate);break
                if len(set(embedded))==1:
                    ids=list(set(embedded));status='Broader disease';method='explicit embedded disease + context retained';candidates=[]
                    note='Links to the broader disease concept; original subtype, resistance or model wording is retained and is not an exact synonym.'
            # The previous LLM proposal is a candidate, not verified merely because ID exists.
            if ids:pass
            elif candidates:
                status='Needs review';method='previous proposed mapping, valid ontology ID';note='Prior proposal retained for review; subtype/broader-term equivalence not automatically accepted.'
            elif re.search(r'\b(?:healthy|normal tissue|normal cells|non.?diseas|HEK293|transplant|in vitro)\b',label,re.I):
                status='Model or context';method='context classification';note='Not treated as a disease diagnosis.'
            else:
                status='Unmapped';method='no exact label or exact-synonym match';note='Retain original term; review disease wording or model context.'
        if status=='Matched' and method=='':method='preferred label'
        urls=['https://www.ebi.ac.uk/ols4/ontologies/'+('mondo' if x.startswith('MONDO') else 'doid')+'/terms?obo_id='+x for x in ids or candidates]
        out.append({'Original disease':label,'Normalized disease':' | '.join(entries[x]['name'] for x in ids) or label,'Disease ID':' | '.join(ids) or '-', 'Status':status,'Web rows':v['row_count'],'Literature rows':v['literature_rows'],'Omics rows':v['omics_rows'],'Method':method,'Candidates':' | '.join(x+' '+entries[x]['name'] for x in candidates),'Notes':note,'Example source':v['sources'][0],'Database URL':' | '.join(urls)})
    write_tsv(CACHE/'disease_mapping.tsv',out);save(CACHE/'disease_mapping.json',out)
    stats=collections.Counter(r['Status'] for r in out);coverage=collections.Counter()
    for r in out:coverage[r['Status']]+=r['Web rows']
    print('Diseases complete:',stats,'row coverage:',coverage,flush=True)
    unresolved=[r for r in out if r['Status'] in ['Unmapped','Needs review','Ambiguous']]
    save(CACHE/'disease_review.json',unresolved)
if __name__=='__main__':main()
