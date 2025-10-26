import pandas as pd
import numpy as np

def gen_base(df, df_demo, df_admin):
    # generate the baseline creatinine level for each subject
    base_dict = {'M':{'BLACK/AFRICAN AMERICAN':{(20,24):1.5, (25,29):1.5, (30,39):1.4, (40,54):1.3, (55,65):1.3, (65,np.inf):1.2},
                      'Others':{(20,24):1.3, (25,29):1.2, (30,39):1.2, (40,54):1.1, (55,65):1.1, (65,np.inf):1}},
                 'F':{'BLACK/AFRICAN AMERICAN':{(20,24):1.2, (25,29):1.1, (30,39):1.1, (40,54):1, (55,65):1, (65,np.inf):0.9},
                      'Others':{(20,24):1, (25,29):1, (30,39):0.9, (40,54):0.9, (55,65):0.8, (65,np.inf):0.8}}}
    df['BASEVALUE'] = np.nan
    subj_lst = np.unique(df['SUBJECT_ID'])
    for subj in subj_lst:
        age = df_demo.loc[df_demo['subject_id']==subj,'anchor_age'].values[0]
        gender = df_demo.loc[df_demo['subject_id']==subj,'gender'].values[0]
        ethnicity = df_admin.loc[df_admin['subject_id']==subj,'race'].values
        if 'BLACK/AFRICAN AMERICAN' in ethnicity:
            age_dict = base_dict[gender]['BLACK/AFRICAN AMERICAN']
        else:
            age_dict = base_dict[gender]['Others']
        if age < 20:
            age = 20
        for key in age_dict.keys():
            low_b, up_b = key
            if low_b <= age <= up_b:
                val = age_dict[key]
                df.loc[df['SUBJECT_ID']==subj,'BASEVALUE'] = val
    df.to_csv('AKI_DATA/AKI_raw_with_base.csv', index=False)
    return ;

def add_s3_lb(df):
    # generate severe AKI labels
    df['s3_lb'] = np.nan
    iid = np.unique(df['ID'])
    for i_id in iid:
        i_df = df[df['ID']==i_id]
        #i_base_df = df_base[df_base['ID']==i_id]
        #base = i_base_df['Base_SCr'].iloc[0]
        base = i_df['BASEVALUE'].iloc[0]
        if np.sum(i_df['Value'] >= 3 * base) > 0 or np.sum(i_df['Value'] >= 4.0) > 0:
            df.loc[df['ID']==i_id, 's3_lb'] = 1
        else:
            df.loc[df['ID']==i_id, 's3_lb'] = 0
    df.to_csv('AKI_DATA/AKI_raw_s3_labels.csv', index=False)
    return ;

def main():
    df_cr = pd.read_csv('AKI_DATA/AKI_data_raw.csv') # change path/name as needed
    df_demo = pd.read_csv('patients.csv') # from MIMIC-IV database
    df_admin = pd.read_csv('admissions.csv') # from MIMIC-IV database
    gen_base(df_cr, df_demo, df_admin)
    df_base = pd.read_csv('AKI_DATA/AKI_raw_with_base.csv')
    add_s3_lb(df_base)
    return ;

main()
