import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
import pandas as pd
import numpy as np
from datetime import datetime

# helper functions:
def get_ref_d(df, admin_typ):
    if admin_typ=='hosp':
        adm_id, intime, outtime = df['hadm_id'].iloc[0], df['admittime'].iloc[0], df['dischtime'].iloc[0]
        subj_id = df['subject_id'].iloc[0]
        if len(adm_id) > 1:
            order = [(t, out_t, aid) for t, out_t, aid in sorted(zip(intime, outtime, adm_id))]
            return subj_id, order #[(order[-1][0], order[-1][1], order[-1][2])]
        else:
            return subj_id, [(intime[0], outtime[0], adm_id[0])]
    elif admin_typ=='icu':
        adm_id, icu_id, intime, outtime = df['hadm_id'].iloc[0], df['stay_id'].iloc[0], df['intime'].iloc[0], df['outtime'].iloc[0]
        subj_id = df['subject_id'].iloc[0]
        if len(icu_id) > 1:
            order = [(t, out_t, iid, aid) for t, out_t, iid, aid in sorted(zip(intime, outtime, icu_id, adm_id))]
            return subj_id, [(order[-1][0], order[-1][1], order[-1][2], order[-1][3])]
        else:
            return subj_id, [(intime[0], outtime[0], icu_id[0], adm_id[0])]
    else:
        print('Admin type not found!')

def get_item_id(f, lab_f, d_item, d_labitem):
    item_id = []
    for name in f:
        key = d_item[d_item['label']==name].iloc[[0],:]
        item_id.append((name, key['itemid'].iloc[0], key['linksto'].iloc[0]))
    for name in lab_f:
        key = d_labitem[(d_labitem['label']==name) & (d_labitem['fluid']=='Blood')].iloc[[0],:]
        #key = d_labitem[d_labitem['label']==name].iloc[[0],:]
        item_id.append((name, key['itemid'].iloc[0], 'lab'))
    return item_id

def get_value(subj_id, icu_id, i_id, ref_d, typ):
    ref_d = datetime.strptime(ref_d, '%Y-%m-%d %H:%M:%S')
    FLAG, unit = False, 'NA'
    if typ == 'lab':
        tmp = lab.loc[(lab['subject_id']==subj_id) & (lab['itemid']==i_id)]
        day = [(datetime.strptime(d, '%Y-%m-%d %H:%M:%S')-ref_d).total_seconds() / 60 for d in list(tmp['charttime'])]
        num = list(tmp['valuenum'].astype('float'))
        # refine data after ref_d
        day = [hr for hr, _ in zip(day, num) if hr >= 0]
        num = [val for hr, val in zip(day, num) if hr >= 0]
        if len(num) > 0:
            unit = tmp['valueuom'].iloc[0]
        else:
            FLAG = True
    elif typ == 'chartevents':
        tmp = chart.loc[(chart['subject_id']==subj_id) & (chart['hadm_id']==icu_id) & (chart['itemid']==i_id)] # note the change of admin type!
        day = [(datetime.strptime(d, '%Y-%m-%d %H:%M:%S')-ref_d).total_seconds() / 60 for d in list(tmp['charttime'])]
        num = list(tmp['valuenum'].astype('float'))
        day = [hr for hr, _ in zip(day, num) if hr >= 0]
        num = [val for hr, val in zip(day, num) if hr >= 0]
        if len(num) > 0:
            unit = tmp['valueuom'].iloc[0]
        else:
            FLAG = True
    elif typ == 'inputevents':
        tmp = intake.loc[(intake['subject_id']==subj_id) & (intake['stay_id']==icu_id) & (intake['itemid']==i_id)]
        start_day = [(datetime.strptime(d, '%Y-%m-%d %H:%M:%S')-ref_d).total_seconds() / 60 for d in list(tmp['starttime'])]
        end_day = [(datetime.strptime(d, '%Y-%m-%d %H:%M:%S')-ref_d).total_seconds() / 60 for d in list(tmp['endtime'])]
        num = list(tmp['amount'].astype('float'))
        s_day = [s_hr for s_hr, _, _ in zip(start_day, end_day, num) if s_hr >= 0]
        e_day = [e_hr for s_hr, e_hr, _ in zip(start_day, end_day, num) if s_hr >= 0]
        num = [val for s_hr, _, val in zip(start_day, end_day, num) if s_hr >= 0]
        if len(num) > 0:
            unit = tmp['amountuom'].iloc[0]
        else:
            FLAG = True
    elif typ == 'outputevents':
        tmp = output.loc[(output['subject_id']==subj_id) & (output['stay_id']==icu_id) & (output['itemid']==i_id)]
        day = [(datetime.strptime(d, '%Y-%m-%d %H:%M:%S')-ref_d).total_seconds() / 60 for d in list(tmp['charttime'])]
        num = list(tmp['value'].astype('float'))
        day = [hr for hr, _ in zip(day, num) if hr >= 0]
        num = [val for hr, val in zip(day, num) if hr >= 0]
        if len(num) > 0:
            unit = tmp['valueuom'].iloc[0]
        else:
            FLAG = True
    else:
        print('Type not found!')
    day = [d for d, _ in sorted(zip(day, num))]
    num = [val for _, val in sorted(zip(day, num))]
    return day, num, FLAG, unit

def if_diagnose(adm_id, subj_id, diag_t):
    flag = False
    diag = diagnose[(diagnose['hadm_id']==adm_id) & (diagnose['subject_id']==subj_id)]
    for i in range(len(diag)):
        code = diag['icd_code'].iloc[i]
        ver = diag['icd_version'].iloc[i]
        name = d_diagnose['long_title'].loc[(d_procedure['icd_code']==code) & (d_procedure['icd_version']==ver)]
        for ele in list(name):
            if diag_t in ele:
                flag = True
    return flag

def if_procedure(subj_id, targets):
    pro_df = procedure[procedure['subject_id']==subj_id]
    flag = False
    for i in range(len(pro_df)):
        code, ver = pro_df.iloc[i,:]['icd_code'], pro_df.iloc[i,:]['icd_version']
        name = d_procedure['long_title'].loc[(d_procedure['icd_code']==code) & (d_procedure['icd_version']==ver)]
        for ele in list(name):
            for tgt in targets:
                if tgt in ele:
                    flag = True
    return flag

def get_lab_data(subj_id, admins, admin_typ, item_id, cut_hr, diag_t):
    subj_data = []
    if admin_typ=='hosp':
        for (ref_d, out_d, adm_id) in admins:
            flag = if_diagnose(adm_id, subj_id, diag_t)
            if flag:
                in_t = datetime.strptime(ref_d, '%Y-%m-%d %H:%M:%S')
                if isinstance(out_d, float):
                    continue
                out_hr = (datetime.strptime(out_d, '%Y-%m-%d %H:%M:%S') - in_t).total_seconds() / 3600
                if (out_hr < 0) | (out_hr > cut_hr):
                    continue
                for name, id, typ in item_id:
                    hour, num, _, _ = get_value(subj_id, adm_id, id, ref_d, typ)
                    for j in range(len(hour)):
                        subj_data.append([subj_id, name, hour[j], num[j], adm_id])
        # add demographic information
        s_df = patient[patient['subject_id']==subj_id]
        subj_data.append([subj_id, 'Age', np.nan, s_df['anchor_age'].values, adm_id])
    elif admin_typ=='icu':
        for (ref_d, out_d, icu_id, adm_id) in admins:
            flag = if_diagnose(adm_id, subj_id, diag_t)
            if flag:
                in_t = datetime.strptime(ref_d, '%Y-%m-%d %H:%M:%S')
                if isinstance(out_d, float):
                    continue
                out_hr = (datetime.strptime(out_d, '%Y-%m-%d %H:%M:%S') - in_t).total_seconds() / 3600
                if (out_hr < 0) | (out_hr > cut_hr):
                    continue
                for name, id, typ in item_id:
                    hour, num, _, _ = get_value(subj_id, icu_id, id, ref_d, typ)
                    for j in range(len(hour)):
                        subj_data.append([subj_id, name, hour[j], num[j], icu_id])
    else:
        print('Admin type not found!')
    return subj_data

def main(admin_typ, cut_hr, diag_t, proc_t):
    if admin_typ=='hosp':
        admin_time = admission.groupby('subject_id')['hadm_id','admittime','dischtime'].agg(list).reset_index()
    elif admin_typ=='icu':
        admin_time = icu_stay.groupby('subject_id')['hadm_id','stay_id','intime','outtime'].agg(list).reset_index()
    else:
        print('Admin type not found!')
    
    lab_features = ['Creatinine']
    # no chartevent features so []
    item_id = get_item_id([], lab_features, d_item, d_labitem)
    print('Start pulling data...')
    all_data = []
    for i in range(len(admin_time)):
        subj_id, admins = get_ref_d(admin_time.iloc[[i],:], admin_typ)
        proc_flag = if_procedure(subj_id, proc_t)
        if proc_flag:
            continue
        subj_data = get_lab_data(subj_id, admins, admin_typ, item_id, cut_hr, diag_t)
        if len(subj_data) > 0:
            all_data.append(pd.DataFrame(subj_data, columns=['SUBJECT_ID','FEATURE_NAME','RECORD_HR','VALUE','ICU_ID']))
    pd.concat(all_data).to_csv('AKI_DATA/AKI_data_raw.csv', index=False)
    return ;

patient = pd.read_csv('patients.csv')
icu_stay = pd.read_csv('icustays.csv')
admission = pd.read_csv('admissions.csv')
lab = pd.read_csv('labevents.csv')
d_labitem = pd.read_csv('d_labitems.csv')
#chart = pd.read_csv('chartevents.csv')
d_item = pd.read_csv('d_items.csv')
diagnose = pd.read_csv('diagnoses_icd.csv')
d_diagnose = pd.read_csv('d_icd_diagnoses.csv')
procedure = pd.read_csv('procedures_icd.csv')
d_procedure = pd.read_csv('d_icd_procedures.csv')

main('icu', 21*24, 'acute kidney injury', ['burn', 'renal dialysis', 'Hemodialysis', 'Peritoneal dialysis', 'ESRD', 'renal transplantation'])
