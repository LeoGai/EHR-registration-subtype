import sys
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

def test_prediction(df, lb, icd_df, subj_id, curve_range, te_sz, method, n_repeat, lmbd, beta, gamma, K, t_max, length, name):
    log_dict = {'auc':[], 'pre':[], 'rec':[], 'acc':[], 'auprc':[]}
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
    for i in range(n_repeat):
        # start prediction
        tr0, te0 = train_test_split(gr0_subj, test_size=te_sz) #random_state=42
        tr1, te1 = train_test_split(gr1_subj, test_size=te_sz)
        tr_subj, te_subj = pd.concat([tr0, tr1]), pd.concat([te0, te1])
        X_train, X_test = X_reg[X_reg['SUBJECT_ID'].isin(tr_subj)], X_reg[X_reg['SUBJECT_ID'].isin(te_subj)]
        y_train, y_test = df_lb[df_lb['SUBJECT_ID'].isin(tr_subj)], df_lb[df_lb['SUBJECT_ID'].isin(te_subj)]
        X_tr, y_tr = compute_feature(X_train, y_train, length, 'time', None)
        X_te, y_te = compute_feature(X_test, y_test, length, 'time', None)
        if method == 'LR':
            clf = LogisticRegression().fit(X_tr, y_tr)
            y_proba = clf.predict_proba(X_te)
            auc = roc_auc_score(y_te, y_proba[:,1])
            y_hat = np.argmax(y_proba, axis=1)
            pre = metrics.precision_score(y_te, y_hat, average='macro')
            rec = metrics.recall_score(y_te, y_hat, average='macro')
            acc = metrics.accuracy_score(y_te, y_hat)
            auroc, auprc = compute_auc(y_te, y_proba, length, plot=False)
            #torch.save(clf, 'AKI_DATA/'+name+'_'+str(length)+'.pt')
            cm = metrics.plot_confusion_matrix(clf, X_te, y_te, normalize=None)
            #cm.figure_.savefig('AKI_DATA/'+name+'_'+str(length)+'_cmat.pdf')
            log_dict['auc'].append(auc)
            log_dict['pre'].append(pre)
            log_dict['rec'].append(rec)
            log_dict['acc'].append(acc)
            log_dict['auprc'].append(auprc)
        else:
            print('Method not defined yet!')
    return log_dict
