import pandas as pd
import numpy as np
from predict_utils import *
import itertools
from multiprocessing import Pool
import torch

curve_range = 21
n_repeat = 500
te_sz = 0.3
raw_data = pd.read_csv('AKI_DATA/AKI_data_recovered.csv') # change path as needed
raw_data = raw_data[['ID', 'Feature', 'Time', 'Value']]
raw_data = raw_data.rename(columns={'ID':'SUBJECT_ID', 'Feature':'feature', 'Time':'time', 'Value':'response'})
raw_data['time'] = round(raw_data['time'])

def filter_data(raw_data, curve_range):
    subid = np.unique(raw_data['SUBJECT_ID'])
    data = []
    for i, id in enumerate(subid):
        temp = raw_data[raw_data['SUBJECT_ID']==id]
        data.append(temp[temp['time'] < curve_range])
    data = pd.concat(data)
    return data, subid

def report_results(log_dict, length):
    for metric in ['auc','pre','rec','acc','auprc']:
        val = log_dict[metric]
        c25 = np.percentile(val, 2.5)
        c97 = np.percentile(val, 97.5)
        print('n_days:', length, 'metric:', metric, 'mean:', np.mean(val), '95CI:', c25, c97)
    return ;

def wrap_fcn(args, name='recover2'):
    data, subid = filter_data(raw_data, curve_range)
    scores = test_prediction(data, True, None, subid, curve_range, te_sz, 'LR', n_repeat, 2, 10, 0.01, None, 200, args, name)
    report_results(scores, args)
    torch.save(scores, 'AKI_DATA/'+name+'_500rep_'+str(args)+'.pt')
    return ;

def main(num_job=4):
    combinations = [3, 4]
    with Pool(num_job) as pool:
        results = pool.map(wrap_fcn, combinations)

if __name__ == '__main__':
    main()
