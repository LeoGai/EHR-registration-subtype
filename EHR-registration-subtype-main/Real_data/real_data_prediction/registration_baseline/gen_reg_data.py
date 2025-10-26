import time
import random
import math
import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split
from registration import *

def perform_registration(df, subj_id, scale_dict, te_sz):
    df_lb = pd.read_csv('AKI_DATA/AKI_raw_s3_labels.csv')
    df_lb = df_lb[['ID','s3_lb']]
    df_lb = df_lb.drop_duplicates(subset=['ID'])
    df_lb = df_lb.rename(columns={'ID':'SUBJECT_ID', 's3_lb':'cluster'})
    df_lb['cluster'] = df_lb['cluster'].astype(int)
    df_lb = df_lb[df_lb['SUBJECT_ID'].isin(subj_id)]

    # registration
    s_time = time.time()
    X_reg, template, features, mu_tmp, sigma_tmp = preprocess(df, scale_dict)
    e_time = time.time() - s_time
    X_reg['time'] = X_reg['time'].astype(int)
    X_reg['register_time'] = X_reg['register_time'].astype(int)
    X_reg.to_csv('AKI_DATA/AKI_reg_base.csv', index=False)
    print('Registration time in sec:', e_time)
    return ;

raw_data = pd.read_csv('AKI_DATA/AKI_data_raw.csv')
data = raw_data[['ID','Feature','Time','Value']].copy()
data = data.rename(columns={'ID':'SUBJECT_ID', 'Feature':'feature', 'Time':'time', 'Value':'response'})
data['time'] = round(data['time'])

subid = np.unique(data['SUBJECT_ID'])
scale_dict = dict()
for s_id in subid:
    temp = data[data['SUBJECT_ID']==s_id]
    scale_dict[s_id] = max(temp['time'])

perform_registration(data, subid, scale_dict, te_sz=0.3)
