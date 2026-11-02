"""Read the deposited workbook and create disposable analysis input tables."""
from pathlib import Path
import csv, hashlib, json
from openpyxl import load_workbook

root=Path.cwd()
path=root/'data/Children_listening_effort_classroom_noise_meta-analysis.xlsx'
if not path.is_file():raise SystemExit('Download the Zenodo workbook into data/ before running the analysis.')
wb=load_workbook(path,read_only=True,data_only=True)
def records(name):
    rows=list(wb[name].values);return [dict(zip(rows[0],r)) for r in rows[1:]]
corpus=records('corpus')
study_number={}
for row in corpus:
    n=int(row['report_id'][1:]);study_number[row['study_id']]=min(n,study_number.get(row['study_id'],n))

target=root/'.cache/data';target.mkdir(parents=True,exist_ok=True)
def save(name,rows):
    with (target/name).open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    return rows
def mapped(row,rename):return {rename.get(k,k):v for k,v in row.items()}
effects=[]
for row in records('effect_sizes'):
    x=mapped(row,{'effect_id':'effect_ID','study_id':'study','noise_statistic_id':'experiment_group_id','control_statistic_id':'control_group_id','m_noise':'m1','sd_noise':'sd1','n_noise':'n1','m_control':'m2','sd_control':'sd2','n_control':'n2','SNR_dB':'SNR','age_years':'age'})
    x['journal_N']=study_number[x['study']]
    x['report_source']=str(int(x.pop('report_id')[1:]))+'.pdf'
    if x['outcome']=='ACC':x['outcome']='IS'
    effects.append(x)
save('analysis_table.csv',effects)
pairs=[]
for row in records('aligned_pairs'):
    rename={'pair_id':'pair_ID','study_id':'study','effect_id_ACC':'effect_ID_IS','effect_id_RT':'effect_ID_RT','SNR_dB':'SNR','age_years':'age','task_ACC':'task_IS','pair_type':'task_group'}
    for prefix in ['ACC','RT']:
        for a,b in [('m_noise','m1'),('sd_noise','sd1'),('n_noise','n1'),('m_control','m2'),('sd_control','sd2'),('n_control','n2')]:rename[prefix+'_'+a]=('IS' if prefix=='ACC' else 'RT')+'_'+b
    x=mapped(row,rename);x.pop('report_id');x['journal_N']=study_number[x['study']]
    if x['task_group']=='cross-task':x['task_group']='dual-task'
    pairs.append(x)
save('aligned_pairs.csv',pairs)
source=[]
for row in records('extracted_statistics'):
    x=mapped(row,{'statistic_id':'source_id','study_id':'study','SNR_dB':'SNR','age_years':'age'})
    x['report_source']=str(int(x.pop('report_id')[1:]))+'.pdf'
    if x['outcome']=='ACC':x['outcome']='IS'
    source.append(x)
save('source_rows.csv',source)
save('corpus_25_reports.csv',corpus)
ohat=[]
for row in records('risk_of_bias'):
    x=mapped(row,{'study_id':'study','reviewer_1':'R1_rating','reviewer_2':'R2_rating','consensus_rating':'final_rating','analysis_rating':'analysed_rating'})
    x['journal_N']=int(x.pop('report_id')[1:]);x['question']=x['question'].replace('ACC','IS');x['outcome']=x['outcome'].replace('ACC','IS')
    ohat.append(x)
save('ohat_items.csv',ohat)
save('search_flow.csv',records('search_flow'))
wb.close()
assert len(effects)==269 and len(pairs)==47 and len(source)==392 and len(corpus)==25 and len(ohat)==250
assert len(set(r['study'] for r in effects))==24
assert sum(r['in_condition_subset']=='TRUE' for r in effects)==238
assert sum(r['task_group']=='same-task' for r in pairs)==30
assert sum(r['task_group']=='dual-task' for r in pairs)==17
(root/'.cache/input-workbook.json').write_text(json.dumps({'file':path.name,'sha256':hashlib.sha256(path.read_bytes()).hexdigest()},indent=2)+'\n')
print('Loaded the deposited workbook: 269 effects, 47 pairs, 24 studies.')
