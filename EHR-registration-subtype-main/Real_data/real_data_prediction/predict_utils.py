import sys
import copy
from tqdm import tqdm
import math
import random
import pandas as pd
import numpy as np
from sklearn import metrics
from sklearn.model_selection import train_test_split
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.metrics import precision_recall_curve
from utils import *
from soft_clustering import *
import matplotlib.pyplot as plt

import torch
import torch.nn as nn
import torch.nn.functional as F
import torch.optim as optim

def compute_auc(y_te, y_proba, length, plot=False):
    n_class = np.unique(y_te)
    auroc, auprc = [], []
    plt.figure(figsize=(10,8))
    for i in n_class:
        y_true = y_te==i
        fpr, tpr, thres = metrics.roc_curve(y_true, y_proba[:,int(i)])
        pre, rec, thresholds = precision_recall_curve(y_true, y_proba[:,int(i)])
        auroc.append(metrics.auc(tpr, fpr))
        auprc.append(metrics.auc(rec, pre))
        plt.plot(rec, pre, label='Class'+str(int(i)))
    if plot:
        plt.legend(loc='best')
        #plt.savefig(path)
    return np.mean(auroc), np.mean(auprc)

def compute_feature(reg, df_lb, length, time, df_info):
    fea_name = np.unique(reg['feature'])
    df = pd.DataFrame()
    for name in fea_name:
        val = reg[reg['feature']==name]['response'].values
        val = np.nan_to_num(val)
        reg.loc[reg['feature']==name, 'response'] = val
        fea_data = gen_subj_data(reg, name, df_lb['SUBJECT_ID'].values, length, time)
        df = pd.concat([df, pd.DataFrame(fea_data)], axis=1)
    return df.values, df_lb['cluster'].values

def metric_compute(y_proba, y_hat, y_te, length):
    auc = roc_auc_score(y_te, y_proba[:,1])
    pre = metrics.precision_score(y_te, y_hat, average='macro')
    rec = metrics.recall_score(y_te, y_hat, average='macro')
    acc = metrics.accuracy_score(y_te, y_hat)
    _, auprc = compute_auc(y_te, y_proba, length, plot=False)
    return auc, pre, rec, acc, auprc

def subgroup_analysis(stay_list, y_proba, y_hat, y_te, length, name, run_id):
    df_demo = pd.read_csv('patients.csv')
    df_icu = pd.read_csv('icustays.csv')
    age_list, gender_list = [], []
    for stay in stay_list:
        subj = df_icu.loc[df_icu['stay_id']==stay,'subject_id'].iloc[0]
        df_age = df_demo.loc[df_demo['subject_id']==subj,'anchor_age']
        age = df_age.values[0]
        gender = df_demo.loc[df_demo['subject_id']==subj,'gender'].values[0] # M or F
        age_list.append(age)
        gender_list.append(gender)
    age_list = np.array(age_list)
    gender_list = np.array([0 if x=='F' else 1 for x in gender_list])
    # age
    age_g1 = (y_proba[age_list < 60], y_hat[age_list < 60], y_te[age_list < 60])
    age_g2 = (y_proba[age_list >= 60], y_hat[age_list >= 60], y_te[age_list >= 60])
    #torch.save((age_g1, age_g2), path)
    age_sc_g1 = metric_compute(age_g1[0], age_g1[1], age_g1[2], length)
    age_sc_g2 = metric_compute(age_g2[0], age_g2[1], age_g2[2], length)
    # gender
    gender_g1 = (y_proba[gender_list==0], y_hat[gender_list==0], y_te[gender_list==0])
    gender_g2 = (y_proba[gender_list==1], y_hat[gender_list==1], y_te[gender_list==1])
    #torch.save((gender_g1, gender_g2), path)
    gender_sc_g1 = metric_compute(gender_g1[0], gender_g1[1], gender_g1[2], length)
    gender_sc_g2 = metric_compute(gender_g2[0], gender_g2[1], gender_g2[2], length)
    return age_sc_g1, age_sc_g2, gender_sc_g1, gender_sc_g2

def test_prediction(df, lb, icd_df, subj_id, curve_range, te_sz, method, n_repeat, lmbd, beta, gamma, K, t_max, length, name):
    pred_dict = {'auc':[], 'pre':[], 'rec':[], 'acc':[], 'auprc':[]}
    log_dict = {'auc':[], 'pre':[], 'rec':[], 'acc':[], 'auprc':[]}
    groups = ['age_g1', 'age_g2', 'gender_g1', 'gender_g2']
    subgr_dict = {g: copy.deepcopy(log_dict) for g in groups}
    if name != 'AKI_base_reg':
        X_reg, _, _, _, _ = preprocess(df, curve_range, icd_df, beta, gamma, lmbd, t_max, lb, plot=False)
        X_reg['time'] = X_reg['time'].astype(int)
    # create labels
    df_lb = pd.read_csv('AKI_DATA/AKI_raw_s3_labels.csv')
    df_lb = df_lb[['ID','s3_lb']]
    df_lb = df_lb.drop_duplicates(subset=['ID'])
    df_lb = df_lb.rename(columns={'ID':'SUBJECT_ID', 's3_lb':'cluster'})
    df_lb['cluster'] = df_lb['cluster'].astype(int)
    df_lb = df_lb[df_lb['SUBJECT_ID'].isin(subj_id)]
    gr0_subj = df_lb.loc[df_lb['cluster']==0, 'SUBJECT_ID']
    gr1_subj = df_lb.loc[df_lb['cluster']==1, 'SUBJECT_ID']
    for i in tqdm(range(n_repeat)):
        # start prediction
        tr0, te0 = train_test_split(gr0_subj, test_size=te_sz, random_state=i) #random_state=42
        tr1, te1 = train_test_split(gr1_subj, test_size=te_sz, random_state=i)
        tr_subj, te_subj = pd.concat([tr0, tr1]), pd.concat([te0, te1])
        if name == 'AKI_base_reg':
            X_train, X_test = df[df['SUBJECT_ID'].isin(tr_subj)], df[df['SUBJECT_ID'].isin(te_subj)]
        else:
            X_train, X_test = X_reg[X_reg['SUBJECT_ID'].isin(tr_subj)], X_reg[X_reg['SUBJECT_ID'].isin(te_subj)]
        y_train, y_test = df_lb[df_lb['SUBJECT_ID'].isin(tr_subj)], df_lb[df_lb['SUBJECT_ID'].isin(te_subj)]
        if name == 'AKI_base_reg':
            X_tr, y_tr = compute_feature(X_train, y_train, length, 'register_time', None)
            X_te, y_te = compute_feature(X_test, y_test, length, 'register_time', None)
        else:
            X_tr, y_tr = compute_feature(X_train, y_train, length, 'time', None)
            X_te, y_te = compute_feature(X_test, y_test, length, 'time', None)

        if method == 'LR':
            clf = LogisticRegression().fit(X_tr, y_tr)
            y_proba = clf.predict_proba(X_te)
            y_hat = np.argmax(y_proba, axis=1)
            #'''
            # subgroup analysis
            age_sc_g1, age_sc_g2, gender_sc_g1, gender_sc_g2 = subgroup_analysis(y_test['SUBJECT_ID'], y_proba, y_hat, y_te,
                    length, name, i)
            for gr in groups:
                for j, metric in enumerate(['auc', 'pre', 'rec', 'acc', 'auprc']):
                    if gr == 'age_g1':
                        subgr_dict[gr][metric].append(age_sc_g1[j])
                    elif gr == 'age_g2':
                        subgr_dict[gr][metric].append(age_sc_g2[j])
                    elif gr == 'gender_g1':
                        subgr_dict[gr][metric].append(gender_sc_g1[j])
                    elif gr == 'gender_g2':
                        subgr_dict[gr][metric].append(gender_sc_g2[j])
                    else:
                        print('Group name not in the list!')
            #'''
            #'''
            auc = roc_auc_score(y_te, y_proba[:,1]) #multi_class='ovo'
            pre = metrics.precision_score(y_te, y_hat, average='macro')
            rec = metrics.recall_score(y_te, y_hat, average='macro')
            acc = metrics.accuracy_score(y_te, y_hat)
            _, auprc = compute_auc(y_te, y_proba, length, plot=False)
            pred_dict['auc'].append(auc)
            pred_dict['pre'].append(pre)
            pred_dict['rec'].append(rec)
            pred_dict['acc'].append(acc)
            pred_dict['auprc'].append(auprc)
            #'''
        else:
            print('Method not defined yet!')
    return subgr_dict, pred_dict
