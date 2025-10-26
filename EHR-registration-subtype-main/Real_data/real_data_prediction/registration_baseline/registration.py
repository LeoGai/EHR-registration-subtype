import pandas as pd
import numpy as np
from itertools import compress
from scipy import interpolate
from cp_als_method import cp_als
import time

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

def MVFR(datlist, features, curve_range=140, n_iter=10, method='mean'):
    num_sub = len(datlist)
    num_feature = len(features)
    t_max, t_min = np.zeros(num_sub), np.zeros(num_sub)
    for i in range(num_sub):
        ts = datlist[i]['time'].values
        t_max[i] = max(ts)
        t_min[i] = min(ts)

    t0_hat = np.zeros(num_sub) # shifts
    if method == 'decomposition':
        for it in range(n_iter):
            # update template curve
            dat_tensor = []
            for j in features:
                record_j = np.nan * np.zeros((num_sub, curve_range)) # feature slice
                for i in range(num_sub):
                    temp = datlist[i][datlist[i]['feature']==j].copy()
                    temp['time'] = temp['time'].values + t0_hat[i]
                    temp[temp['time'] > curve_range-1]['time'] = curve_range-1
                    temp[temp['time'] < 0]['time'] = 0
                    t_step = [int(t) for t in temp['time'].values]
                    record_j[i, t_step] = temp['response'].values
                dat_tensor.append(record_j)
            dat_tensor = np.concatenate(dat_tensor).reshape(num_feature, num_sub, -1) # feature, subject, time

            # compute template curve
            omega = np.isnan(dat_tensor)
            _, dat_tensor_recover = cp_als(dat_tensor, 10, omega)

            # update t0_hat
            total_loss = 0
            for i in range(num_sub):
                loss = np.inf
                for t0_candidate in range(0, int(curve_range-t_length[i])): # modified
                    temp = datlist[i].copy()
                    temp['time'] = temp['time'].values + t0_candidate
                    this_loss = 0
                    for idx, j in enumerate(features):
                        record_j = temp[temp['feature'] == j]
                        this_loss += np.sum((record_j['response'].values - dat_tensor_recover[idx,i,record_j['time'].values.astype(int)])**2)
                    if this_loss <= loss:
                        loss = this_loss
                        t0_hat[i] = t0_candidate
                total_loss += loss
            print(total_loss/num_sub)
        return t0_hat

    else:
        mu_hat = np.zeros((num_feature, curve_range)) # template
        record_save = [] # lst of features->lst of subject df
        for it in range(n_iter):
        # update template curve
            #s_time = time.time()
            for idx, j in enumerate(features):
                mu_hat[idx,:] = np.nan
                record_j, r_j = [], []
                for i in range(num_sub):
                    s_df = datlist[i][datlist[i]['feature']==j]
                    r_j.append(s_df)
                    temp = s_df.copy()
                    temp['time'] = temp['time'].values + t0_hat[i]
                    record_j.append(temp)
                # compute template curve
                record_j = pd.concat(record_j)
                record_save.append(r_j)
                t_pts = np.unique(record_j['time'].values)
                for t in t_pts:
                    mu_hat[idx,int(t)] = np.mean(record_j[record_j['time']==t]['response'].values)
                mu_hat[idx] = intextplt(mu_hat[idx])

            # update t0_hat
            total_loss = 0
            for i in range(num_sub):
                loss = np.inf
                for t0_candidate in range(int(-t_min[i]), int(curve_range-t_max[i])): # modified
                    this_loss = 0
                    for idx, j in enumerate(features):
                        record_j = record_save[idx][i]
                        record_j['time'] = record_j['time'].values + t0_candidate
                        this_loss += np.sum((record_j['response'].values - mu_hat[idx,record_j['time'].values.astype(int)])**2)
                        record_j['time'] = record_j['time'].values - t0_candidate
                    if this_loss <= loss:
                        loss = this_loss
                        t0_hat[i] = t0_candidate
                total_loss += loss
            #print(total_loss/num_sub)
            #e_time = time.time() - s_time
            #print('Time in sec: ', e_time)
        return mu_hat, t0_hat

def preprocess(df, sc_dict):
    df = df.dropna(subset=['response'])
    t_range = 21 #120
    curve_range = 28 #140
    subid = np.unique(df['SUBJECT_ID'].values)
    num_sub = len(subid)
    features = np.unique(df['feature'].values)
    datlist = []

    # normalization
    featureGroup = df.groupby('feature')
    mu_lst, std_lst = featureGroup.transform('mean').response.values, featureGroup.transform('std').response.values
    mu_tmp = featureGroup['response'].agg(['mean'])['mean'].to_numpy().reshape((1,len(features)))
    sigma_tmp = featureGroup['response'].agg(['std'])['std'].to_numpy().reshape((1,len(features)))
    df.loc[:,'response'] = (df.loc[:,'response'].values - mu_lst) / std_lst
    # re-scale
    for i, id in enumerate(subid):
        temp = df[df['SUBJECT_ID']==id]
        temp['response'] = temp['response'] / np.maximum(np.ones(temp['time'].shape), (sc_dict[id] - temp['time']))
        datlist.append(temp.drop(columns=['SUBJECT_ID']))
    #datlist = [x.drop(columns=['SUBJECT_ID']) for _, x in df.groupby('SUBJECT_ID')]

    # update features
    template, shifts = MVFR(datlist, features, curve_range, method='mean')
    register_data = []
    for i in range(num_sub):
        temp = datlist[i]
        temp['register_time'] = temp['time'].values + shifts[i]
        temp['SUBJECT_ID'] = [subid[i]] * len(temp)
        register_data.append(temp)
    register_data = pd.concat(register_data)
    return register_data, template, features, mu_tmp, sigma_tmp

def fit_template(df, template, full_features, mu_tmp, sigma_tmp, sc_dict):
    df = df.dropna(subset=['response'])
    t_range = 21 #120
    curve_range = 28 #140
    subid = np.unique(df['SUBJECT_ID'].values)
    num_sub = len(subid)
    datlist = []

    # normalization
    features = np.unique(df['feature'].values)
    for name in features:
        assert name in full_features
        f_mu = mu_tmp[0, full_features==name][0]
        f_std = sigma_tmp[0, full_features==name][0]
        df.loc[df['feature']==name,'response'] = (df.loc[df['feature']==name,'response'].values - f_mu) / f_std
    # re-scale
    for i, id in enumerate(subid):
        temp = df[df['SUBJECT_ID']==id]
        temp['response'] = temp['response'] / np.maximum(np.ones(temp['time'].shape), (sc_dict[id] - temp['time']))
        datlist.append(temp.drop(columns=['SUBJECT_ID']))
    #datlist = [x.drop(columns=['SUBJECT_ID']) for _, x in df.groupby('SUBJECT_ID')]

    # update features
    t_length = [(max(temp['time'].values), min(temp['time'].values)) for temp in datlist]
    # find best t0
    t0_hat = np.zeros(num_sub) # shifts
    total_loss = 0
    for i in range(num_sub):
        loss = np.inf
        temp = datlist[i]
        for t0_candidate in range(int(-t_length[i][1]), int(curve_range - t_length[i][0])): # modified
            temp['time'] = temp['time'].values + t0_candidate
            this_loss = 0
            for j in features:
                record_j = temp[temp['feature']==j]
                this_loss += np.sum((record_j['response'].values - template[full_features==j,record_j['time'].values.astype(int)])**2)
            temp['time'] = temp['time'].values - t0_candidate
            if this_loss <= loss:
                loss = this_loss
                t0_hat[i] = t0_candidate
        total_loss += loss
    print('Test: ', total_loss/num_sub)

    # return registered data
    register_data = []
    for i in range(num_sub):
        temp = datlist[i]
        temp['register_time'] = temp['time'].values + t0_hat[i]
        temp['SUBJECT_ID'] = [subid[i]] * len(temp)
        register_data.append(temp)
    register_data = pd.concat(register_data)
    return register_data
