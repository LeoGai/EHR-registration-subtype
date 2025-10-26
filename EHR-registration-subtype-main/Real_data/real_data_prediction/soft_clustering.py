import sys
import math
import pandas as pd
import numpy as np
from sklearn.cluster import KMeans
from scipy import interpolate
from scipy.spatial import distance
from scipy.stats import f_oneway
from scipy.io import savemat
import statsmodels.api as sm
import torch
import time
import string
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from utils import *
from soft_cluster_utils import *
from hard_cluster_utils import init_temp

def normalize(df, train=True, mu_tmp=None, sigma_tmp=None, full_features=None, if_filter=True):
    features = np.unique(df['feature'])
    if train:
        if if_filter==False:
            featureGroup = df.groupby('feature')
            mu_lst, std_lst = featureGroup.transform('mean').response.values, featureGroup.transform('std').response.values
            mu_tmp = featureGroup['response'].agg(['mean'])['mean'].to_numpy().reshape((1,len(features)))
            sigma_tmp = featureGroup['response'].agg(['std'])['std'].to_numpy().reshape((1,len(features)))
            df.loc[:,'response'] = (df.loc[:,'response'].values - mu_lst) / std_lst
            return df, mu_tmp, sigma_tmp, features
        else:
            featureGroup = df.groupby('feature')
            mu_lst, std_lst = featureGroup.transform('mean').response.values, featureGroup.transform('std').response.values
            mu_tmp = featureGroup['response'].agg(['mean'])['mean'].to_numpy().reshape((1,len(features)))
            sigma_tmp = featureGroup['response'].agg(['std'])['std'].to_numpy().reshape((1,len(features)))
            upperlim = mu_lst + 3 * std_lst
            lowerlim = mu_lst - 3 * std_lst
            lowidx = df[df.loc[:,'response'].values > upperlim].index
            upidx = df[df.loc[:,'response'].values < lowerlim].index
            df = df.drop(upidx|lowidx)
            #df.to_csv('lab_drg_870_872_filtered_Feb_iv.csv', index=False)
            featureGroup = df.groupby('feature')
            mu_lst, std_lst = featureGroup.transform('mean').response.values, featureGroup.transform('std').response.values
            mu_tmp = featureGroup['response'].agg(['mean'])['mean'].to_numpy().reshape((1,len(features)))
            sigma_tmp = featureGroup['response'].agg(['std'])['std'].to_numpy().reshape((1,len(features)))
            df.loc[:,'response'] = (df.loc[:,'response'].values - mu_lst) / std_lst
            return df, mu_tmp, sigma_tmp, features
    else:
        for name in features:
            assert name in full_features
            f_mu = mu_tmp[0, full_features==name][0]
            f_std = sigma_tmp[0, full_features==name][0]
            df.loc[df['feature']==name,'response'] = (df.loc[df['feature']==name,'response'].values - f_mu) / f_std
        return df, features

def scale(df, subid, sc_dict):
    datlist, new_id = [], []
    for i, id in enumerate(subid):
        temp = df[df['SUBJECT_ID']==id]
        #temp['response'] = temp['response'] / np.maximum(np.ones(temp['time'].shape), (sc_dict[id] - temp['time']))
        if len(temp) > 0:
            datlist.append(temp)
            new_id.append(id)
    return datlist, new_id

def makeTensor(datlist, features, curve_range):
    num_sub = len(datlist)
    dat_tensor, raw_indices = [], []
    for j in features:
        record_j = np.nan * np.zeros((num_sub, curve_range)) # feature slice
        for i in range(num_sub):
            temp = datlist[i][datlist[i]['feature']==j]
            t_step = [int(t) for t in temp['time'].values] # record the t_step with records
            record_j[i, t_step] = temp['response'].values
        dat_tensor.append(record_j)
    dat_tensor = np.concatenate(dat_tensor).reshape(len(features), num_sub, -1) # feature, subject, time
    for i in range(num_sub):
        mat = dat_tensor[:,i,:]
        index = ~np.isnan(mat) # boolean mask
        raw_indices.append(index)
    return dat_tensor, np.array(raw_indices)

def intextplt(seq):
    x_bool = ~np.isnan(seq)
    if sum(x_bool) == 0:
        return seq
    else:
        if np.isnan(seq[0]):
            seq[0] = seq[x_bool][0]
        if np.isnan(seq[-1]):
            seq[-1] = seq[x_bool][-1]
        # update after filling ends
        x_bool = ~np.isnan(seq)
        x = np.arange(len(seq))[x_bool]
        y = seq[x_bool]
        f = interpolate.interp1d(x, y, kind='linear')
        seq = f(np.arange(len(seq)))
        return seq

def makeDF(tensor_est, subid, features):
    num_t = tensor_est.shape[-1]
    df_filled = []
    for i, id in enumerate(subid):
        for j, feat in enumerate(features):
            df_feat = pd.DataFrame(columns=['SUBJECT_ID','feature','time','response'])
            df_feat['SUBJECT_ID'] = [id] * num_t
            df_feat['feature'] = [feat] * num_t
            df_feat['time'] = np.arange(num_t)
            df_feat['response'] = tensor_est[j,i,:]
            df_filled.append(df_feat)
    df_filled = pd.concat(df_filled)
    return df_filled

def preprocess(df, curve_range, icd_df, beta, gamma, eta, t_max, label=False, plot=True):
    df = df.dropna(subset=['response'])
    df, mu, sigma, features = normalize(df, train=True, mu_tmp=None, sigma_tmp=None, full_features=None)
    subid = np.unique(df['SUBJECT_ID'])
    datlist, _ = scale(df, subid, sc_dict=None)
    # missing value imputation
    dat_tensor, _ = makeTensor(datlist, features, curve_range)
    
    omega = ~np.isnan(dat_tensor)
    _, tensor_est = init_temp(dat_tensor, omega, None, len(features), curve_range, K=3)
    
    df_filled = makeDF(tensor_est, subid, features)
    # return data
    register_data = []
    for i, id in enumerate(subid):
        temp = df_filled[df_filled['SUBJECT_ID']==id]
        register_data.append(temp)
    register_data = pd.concat(register_data)
    if label == False:
        print('Not applicable for this task!')
        return ;
    else:
        return register_data, None, None, None, None
