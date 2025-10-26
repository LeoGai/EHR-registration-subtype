import numpy as np
import pandas as pd
from sklearn.cluster import KMeans
from scipy.spatial import distance
import torch
import torch.nn as nn
import torch.optim as optim
from utils import *

def init_temp(dat_tensor, omega, dat_med, num_fea, curve_range, K, method='ALS'):
    k_template = []
    tensor_est = np.zeros(dat_tensor.shape)
    dat_tensor_temp = dat_tensor.copy()
    dat_tensor_temp[~omega] = 0
    for t in range(num_fea):
        matrix = dat_tensor_temp[t,:,:].copy()
        if method == 'ALS':
            alpha = np.nanmean(matrix, axis=1)
            beta = np.nanmean(matrix, axis=0)
            alpha[np.isnan(alpha)] = 0
            beta[np.isnan(beta)] = 0
            alpha = np.expand_dims(alpha, axis=1)
            alpha = np.repeat(alpha, curve_range, axis=1)
            matrix = matrix - alpha
            matrix_est = ALS_reg_completion(matrix,omega[t,:,:],rank=min(curve_range,11),time_range=curve_range,iter_num=5,penalty=1)
            matrix_est = matrix_est + alpha
        else:
            print('Method does not exist!')
        tensor_est[t,:,:] = matrix_est
    # count obs
    count = np.sum(omega, axis=2)
    count_allfea = np.sum(count, axis=0)
    threshold = np.mean(count_allfea) #100
    idx = np.where(count_allfea > threshold)[0]
    template_i = idx[np.random.randint(np.shape(idx)[0])]
    k_template.append(tensor_est[:,template_i,:])
    # select subjects
    subnum = np.sum(count_allfea > threshold)
    for i in range(1,K):
        distance = np.zeros(tensor_est.shape[1])
        for j in range(tensor_est.shape[1]):
            if count_allfea[j] < threshold:
                continue
            loss = tensor_est[:,j,:] - k_template[-1]
            distance[j] = np.trace(loss.dot(loss.T))
        prob = distance / np.sum(distance)
        np.nan_to_num(prob, nan=0)
        template_i = np.random.choice(tensor_est.shape[1], p=prob)
        k_template.append(tensor_est[:,template_i,:])
    return k_template, tensor_est
