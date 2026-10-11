"""Add simple grouped statistics to the existing download and eligible web figures.
Run with --prepare, then render updated figures with R, then --package.
"""
from pathlib import Path
from zipfile import ZipFile
from collections import defaultdict
import importlib.util,io,json,csv,shutil,hashlib,subprocess,sys
import pandas as pd
import numpy as np
from scipy import stats
import scipy

WORK=Path(__file__).resolve().parent
ROOT=WORK.parents[4]
PACKAGE=ROOT/'release_package.zip'
OLD_FIGURES=ROOT/'figure_metadata'
SOURCE=ROOT/'resistance_source'
STAGE=WORK/'staging';BACK=WORK/'backup'
TABLE='downloadfile/06_drug_resistance_1614_entries.tsv'
KEYS=['Kinase','Drug','Cancer name','Sub-tissue']
NEW_COLUMNS=['n_background','n_variant','n_reference','Grouped_FC','P_value_MWU','FDR_BH']

def text_table(rows,columns):
 out=io.StringIO(newline='');writer=csv.DictWriter(out,fieldnames=columns,delimiter='\t',lineterminator='\n');writer.writeheader();writer.writerows(rows);return out.getvalue().encode('utf-8')
def finite(value):return f'{float(value):.10g}' if pd.notna(value) else ''
def hashfile(path):
 h=hashlib.sha256()
 with path.open('rb') as f:
  for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
 return h.hexdigest()

def prepare():
 BACK.mkdir(exist_ok=True);STAGE.mkdir(exist_ok=True)
 backup=BACK/PACKAGE.name
 if not backup.exists():shutil.copy2(PACKAGE,backup)
 with ZipFile(backup) as z:
  original=z.read(TABLE);manifest=pd.read_csv(io.BytesIO(z.read('drug_resistance_figures/manifest.tsv')),sep='\t',dtype=str).fillna('')
 metadata=pd.read_csv(OLD_FIGURES/'figure_metadata.tsv',sep='\t');assert metadata.evidence_id.is_unique
 assert manifest['Evidence ID'].tolist()==metadata.evidence_id.tolist()
 reader=csv.DictReader(io.StringIO(original.decode('utf-8-sig')),delimiter='\t');columns=list(reader.fieldnames);download=list(reader)
 assert len(download)==len(metadata)==1614
 # Check stable identity against the original plotting payload, never rely on row order alone.
 for row,meta in zip(download,metadata.to_dict('records')):
  for col,mcol in [('Treatment','treatment'),('ID of cell line','cell_line_model_id'),('Cell name','cell_line_name'),('Hotspot','hotspot'),('Mutations','mutation')]:assert row[col]==str(meta[mcol]),(col,row,meta)
 spec=importlib.util.spec_from_file_location('source_resistance',SOURCE/'compute_drkinase_resistance.py');cdr=importlib.util.module_from_spec(spec);spec.loader.exec_module(cdr)
 model=cdr.load_model_metadata();hot,_=cdr.map_mutations_to_predicted_hotspots(model);response,_=cdr.load_response(model,cdr.load_cancer_code_map())
 cell=(response.groupby(KEYS+['Cell line ID'],dropna=False).agg(ic50=('IC50 (uM)','median'),n_measurements=('IC50 (uM)','size')).reset_index())
 hot_by=hot.groupby('kinase').response_cell_id.agg(set).to_dict()
 mutation_sets=hot.groupby(['kinase','AAchange']).response_cell_id.agg(set).to_dict()
 muts_by_kinase=defaultdict(dict)
 for (kin,mutation),ids in mutation_sets.items():muts_by_kinase[kin][mutation]=ids
 egfr=hot[hot.kinase.eq('EGFR')].groupby('response_cell_id').AAchange.agg(set)
 double=set(egfr[egfr.map(lambda s:{'p.L858R','p.T790M'}.issubset(s))].index)
 primary=set(egfr[egfr.map(lambda s:'p.L858R' in s and 'p.T790M' not in s)].index)
 background={tuple(k):frame for k,frame in cell.groupby(KEYS,dropna=False,sort=True)}
 test_rows=[];details=[];counts={};percell={}
 def compare(key,group_id,pos_ids,neg_ids,comparison):
  frame=background[key].set_index('Cell line ID');pos=frame.loc[sorted(pos_ids)];neg=frame.loc[sorted(neg_ids)]
  lookup=(comparison,*key,group_id);counts[lookup]=(len(pos),len(neg))
  if len(pos)<3 or len(neg)<3:return
  assert not pos_ids.intersection(neg_ids)
  x=pos.ic50.to_numpy(float);y=neg.ic50.to_numpy(float)
  result=stats.mannwhitneyu(x,y,alternative='two-sided',method='auto')
  test_rows.append(dict(zip(KEYS,key))|{'comparison':comparison,'Mutation':group_id,'n_variant':len(x),'n_reference':len(y),'variant_median_uM':np.median(x),'reference_median_uM':np.median(y),'Grouped_FC':np.median(x)/np.median(y),'P_value_MWU':result.pvalue})
  percell[lookup]=(pos,neg)
 # Include all eligible same-variant tests before selecting released FC-screen candidates.
 for key,frame in background.items():
  kin=key[0];ids=set(frame['Cell line ID']);neg_ids=ids-hot_by.get(kin,set())
  for mutation,carriers in muts_by_kinase.get(kin,{}).items():
   pos_ids=ids&carriers
   if pos_ids:compare(key,mutation,pos_ids,neg_ids,'Treatment-background')
  if kin=='EGFR' and ids&double:compare(key,'p.L858R+p.T790M',ids&double,ids&primary,'Primary-secondary')
 tests=pd.DataFrame(test_rows)
 tests['FDR_BH']=np.nan
 for comparison,idx in tests.groupby('comparison').groups.items():
  tests.loc[idx,'FDR_BH']=stats.false_discovery_control(tests.loc[idx,'P_value_MWU'].to_numpy(),method='bh')
 assert len(tests.query('comparison == "Treatment-background"'))==41
 lookup={(r['comparison'],*(r[k] for k in KEYS),r['Mutation']):r for r in tests.to_dict('records')}
 combined=[]
 for row,meta in zip(download,metadata.to_dict('records')):
  comparison='Primary-secondary' if meta['evidence_id'].startswith('KMD2_DRK_PS_') else 'Treatment-background'
  key=(meta['kinase'],meta['drug'],meta['cancer_name'],meta['sub_tissue'])
  mutation='p.L858R+p.T790M' if comparison=='Primary-secondary' else meta['mutation'];k=(comparison,*key,mutation)
  npos,nneg=counts[k];hit=lookup.get(k)
  ids=set(background[key]['Cell line ID'])
  nbg=len(ids&primary) if comparison=='Primary-secondary' else len(ids)
  assert nbg==int(meta['expected_g2_cells'])
  row['n_background']=nbg;row['n_variant']=npos;row['n_reference']=nneg
  for col in ['Grouped_FC','P_value_MWU','FDR_BH']:row[col]=finite(hit[col]) if hit else ''
  m=dict(meta);m['n_background']=nbg;m['n_variant']=npos;m['n_reference']=nneg
  for col in ['Grouped_FC','P_value_MWU','FDR_BH']:m[col]=hit[col] if hit else np.nan
  m['group_key']=json.dumps(k);combined.append(m)
  if hit:
   assert npos>=3 and nneg>=3
   for group,frame in zip(['Variant','Reference'],percell[k]):
    for cid,value in frame.iterrows():details.append({'evidence_id':meta['evidence_id'],'plot_group':group,'cell_id':cid,'ic50':value.ic50})
 updated=pd.DataFrame(combined);eligible=updated[updated.P_value_MWU.notna()]
 assert len(eligible)==91 and eligible.group_key.nunique()==26
 for old,new in zip(csv.DictReader(io.StringIO(original.decode('utf-8-sig')),delimiter='\t'),download):
  assert all(old[c]==new[c] for c in columns)
 destination=STAGE/TABLE;destination.parent.mkdir(parents=True,exist_ok=True);destination.write_bytes(text_table(download,columns+NEW_COLUMNS))
 tests.to_csv(WORK/'直接统计结果.tsv',sep='\t',index=False,float_format='%.12g')
 eligible.to_csv(WORK/'figure_metadata.tsv',sep='\t',index=False,float_format='%.12g')
 pd.DataFrame(details).to_csv(WORK/'group_cell_values.tsv',sep='\t',index=False,float_format='%.12g')
 # Original full plotting data retain their focal/background definitions.
 shutil.copy2(OLD_FIGURES/'figure_plot_detail.tsv.gz',WORK/'figure_plot_detail.tsv.gz')
 (STAGE/'drug_resistance_figures/png').mkdir(parents=True,exist_ok=True);(STAGE/'drug_resistance_figures/svg').mkdir(parents=True,exist_ok=True)
 summary={'scipy_version':scipy.__version__,'candidate_rows':1614,'candidate_combinations':1501,'eligible_candidate_combinations':26,'eligible_candidate_rows':91,'unchanged_candidate_figures':1523,'all_tests_for_BH':len(tests),'test':'Two-sided Mann-Whitney U, method=auto; one median raw-scale IC50 per independent cell line','reference':'No mapped predicted hotspot mutation in the target kinase; EGFR primary-secondary uses L858R-positive/T790M-negative reference','fdr_scope':'All eligible same-variant tests within matched kinase/drug/cancer/sub-tissue backgrounds, before FC-based candidate selection','download_columns_added':NEW_COLUMNS,'original_package_sha256':hashfile(backup)}
 (WORK/'summary.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2))
 print(json.dumps(summary,ensure_ascii=False,indent=2),flush=True)


def package():
 backup=BACK/PACKAGE.name;summary=json.loads((WORK/'summary.json').read_text());metadata=pd.read_csv(WORK/'figure_metadata.tsv',sep='\t');replacements={TABLE:(STAGE/TABLE).read_bytes()}
 with ZipFile(backup) as old:
  oldnames=old.namelist();checks=pd.read_csv(io.BytesIO(old.read('drug_resistance_figures/checksums.tsv')),sep='\t')
  for row in metadata.to_dict('records'):
   for col,kind in [('png_file','PNG'),('svg_file','SVG')]:
    member='drug_resistance_figures/'+row[col];file=STAGE/member;assert file.exists(),file
    data=file.read_bytes();assert len(data)>1000;replacements[member]=data
    match=checks['Evidence ID'].eq(row['evidence_id'])&checks.type.eq(kind);assert match.sum()==1
    checks.loc[match,'size_bytes']=len(data);checks.loc[match,'md5']=hashlib.md5(data).hexdigest()
  replacements['drug_resistance_figures/checksums.tsv']=checks.to_csv(sep='\t',index=False,lineterminator='\n').encode()
  # Write only replacement members; zip preserves other compressed members and all paths.
  for member,data in replacements.items():
   path=STAGE/member;path.parent.mkdir(parents=True,exist_ok=True);path.write_bytes(data)
  assert hashfile(PACKAGE)==summary['original_package_sha256'],'Current package changed since backup'
  temp=WORK/'待验证发布包.zip';shutil.copy2(PACKAGE,temp)
  subprocess.run(['zip','-q',str(temp.resolve()),*replacements],cwd=STAGE,check=True)
  with ZipFile(temp) as new:
   assert len(new.namelist())==len(set(new.namelist())) and set(new.namelist())==set(oldnames)
   assert new.testzip() is None
   for member in oldnames:
    if member in replacements:assert new.read(member)==replacements[member]
    else:assert new.getinfo(member).CRC==old.getinfo(member).CRC and new.getinfo(member).file_size==old.getinfo(member).file_size,member
   # Noneligible images must be byte-identical; eligible filenames remain unchanged.
   unchanged_pngs=sum(n.startswith('drug_resistance_figures/png/') and n not in replacements for n in oldnames)
   assert unchanged_pngs==1523
  temp.replace(PACKAGE)
 local=ROOT/'web_download/06_drug_resistance_1614_entries.tsv'
 if not (BACK/local.name).exists():shutil.copy2(local,BACK/local.name)
 local.write_bytes(replacements[TABLE])
 summary.update(package_updated=True,package_sha256=hashfile(PACKAGE),changed_members=len(replacements),eligible_pngs_replaced=91,eligible_svgs_replaced=91,noneligible_pngs_preserved=1523,other_members_preserved=True)
 (WORK/'summary.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2));print(json.dumps(summary,ensure_ascii=False,indent=2),flush=True)

if __name__=='__main__':
 if '--prepare' in sys.argv:prepare()
 elif '--package' in sys.argv and '--approved' in sys.argv:package()
 else:raise SystemExit('Use --prepare or --package')
