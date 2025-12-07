import pandas as pd
import numpy as np
from predict_utils import *
import itertools
from multiprocessing import Pool
import torch

curve_range = 21
n_repeat = 500
te_sz = 0.3
raw_data = pd.read_csv(path) # change path as needed

def filter_data(raw_data, curve_range, name):
    subid = np.unique(raw_data['SUBJECT_ID'])
    data = []
    for i, id in enumerate(subid):
        temp = raw_data[raw_data['SUBJECT_ID']==id]
        if name == 'AKI_base_reg':
            data.append(temp[temp['register_time'] < curve_range])
        else:
            data.append(temp[temp['time'] < curve_range])
    data = pd.concat(data)
    return data, subid

def report_results(log_dict, pred_dict, length, name):
    # subgroup
    final_df = []
    groups = ['age_g1', 'age_g2', 'gender_g1', 'gender_g2']
    for gr in groups:
        gr_sc = []
        for metric in ['auc', 'pre', 'rec', 'acc', 'auprc']:
            val = log_dict[gr][metric]
            mu = np.mean(val)
            c25 = np.percentile(val, 2.5)
            c97 = np.percentile(val, 97.5)
            gr_sc.extend([mu, c25, c97])
        final_df.append(gr_sc)
    final_df = pd.DataFrame(final_df, columns=['auc_mu','auc_25','auc_97','pre_mu','pre_25','pre_97',
        'rec_mu','rec_25','rec_97','acc_mu','acc_25','acc_97','auprc_mu','auprc_25','auprc_97'], index=groups)
    #final_df.to_csv(path)
    # prediction
    pred_df = []
    for metric_p in ['auc', 'pre', 'rec', 'acc', 'auprc']:
        val_p = pred_dict[metric_p]
        mu_p = np.mean(val_p)
        c25_p = np.percentile(val_p, 2.5)
        c97_p = np.percentile(val_p, 97.5)
        pred_df.append([metric_p, mu_p, c25_p, c97_p])
    pred_df = pd.DataFrame(pred_df, columns=['metric','mean','c25','c97'])
    #pred_df.to_csv(path)
    return ;

def wrap_fcn(args, name='AKI_recover2'): # name options: AKI_origin, AKI_base_reg, AKI_recover, and AKI_recover2
    if name != 'AKI_base_reg':
        raw_data = raw_data[['ID', 'Feature', 'Time', 'Value']]
        raw_data = raw_data.rename(columns={'ID':'SUBJECT_ID', 'Feature':'feature', 'Time':'time', 'Value':'response'})
        raw_data['time'] = round(raw_data['time'])
    data, subid = filter_data(raw_data, curve_range, name)
    subid = np.unique(data['SUBJECT_ID'])
    sub_sc, pred_sc = test_prediction(data, True, None, subid, curve_range, te_sz, 'LR', n_repeat, 2, 10, 0.01, None, 200, args, name)
    #torch.save(sub_sc, path)
    report_results(sub_sc, pred_sc, args, name)
    return ;

def main():
    wrap_fcn(4)

if __name__ == '__main__':
    main()
