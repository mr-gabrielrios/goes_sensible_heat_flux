'''
Error Analysis Script
Objective: Contain functions to perform error analysis on sensible heat flux model against ground observations and WRF data
'''

from datetime import datetime
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import rcParams as rc
from matplotlib import colors
from scipy import stats
from scipy.stats import gaussian_kde
import os

######################################################

def data_read(url):
    df = pd.DataFrame()
    file_list = sorted(os.listdir(url), key = lambda x: x.split()[0])
    n, c = [len(file_list), 0]
    
    for file in file_list:
        if c < n:
            # print(file)
            temp = pd.read_csv(os.path.join(url, file))
            df = df.append(temp)
        c += 1
    return df

######################################################

### Nash-Sutcliffe Coefficient (NSC) function
def nsc(date, meas, pred):
    terms = [0] * 2
    
    for i in range(0, len(date)):
        if ~np.isnan(meas[i]):
            terms[0] += (meas[i] - np.nanmean(meas))**2
        else:
            continue
        
    for j in range(0, len(date)):
        if ~np.isnan(meas[j]) and ~np.isnan(pred[j]):
            terms[1] += (meas[j] - pred[j])**2
        else:
            continue
    
    e = (terms[0] - terms[1])/terms[0]
    
    return e

######################################################

### Root Mean Squared Error (RMSE) function
def rmse(date, meas, pred):
    num = 0
    
    for i in range(0, len(date)):
        if ~np.isnan(meas[i]) and ~np.isnan(pred[i]):
            num += (meas[i] - pred[i])**2
        else:
            continue
        
    r = np.sqrt(num/len(date))
    
    return r

######################################################

### Mean Bias Error (MBE) function
def mbe(date, meas, pred):
    num = 0
    for i in range(0, len(date)):
        if ~np.isnan(meas[i]) and ~np.isnan(pred[i]):
            num += meas[i] - pred[i]
        else:
            continue
        
    m = num/len(date)

    return m

######################################################

### R^2 Regression Coefficient function
def r_sq(data, var):
    var_meas = var + '_obs'
    
    meas = np.array(data[var_meas])
    pred = np.array(data[var])
    mask = ~np.isnan(data[var_meas]) & ~np.isnan(data[var])
    
    res = stats.linregress(meas[mask], pred[mask])
    
    return res, res.rvalue

def error_plot(var, meas, pred, errs, rgsn, date_range):
    
    [start_date, end_date, start_hr, end_hr] = date_range
    
    err_names = ['RMSE', 'NSC', 'MBE', '$R^2$', 'N']
    a = np.linspace(np.nanmin(meas), np.nanmax(meas), 100)
    b = a
    
    # Filter out bad data (non-finite, such as nans and infs)
    meas = np.asarray(meas)
    pred = np.asarray(pred)
    df = pd.DataFrame({'meas': meas, 
                       'pred': pred})
    df = df.dropna()
    
    x = np.asarray(df['meas'].tolist())
    y = np.asarray(df['pred'].tolist())
    
    ### Heatmap generation
    # Calculate the point density
    xy = np.vstack([x,y])
    z = gaussian_kde(xy)(xy)
    # Sort the points by density, so that the densest points are plotted last
    idx = z.argsort()
    x, y, z = x[idx], y[idx], z[idx]
    
    ### Set up figure    
    # Plot formatting
    rc['font.family'] = 'Segoe UI'    
    fig, ax = plt.subplots(figsize=[8, 6], dpi=144)
    ax.set_aspect('equal')
    ax.set_xlabel('Observed $' + var + '$')
    ax.set_ylabel('Model $' + var + '$')
    # Plot data
    ax.scatter(x, y, c=z, cmap='jet', s=50, edgecolor='')
    ax.plot(a, b, linestyle='--', color='k')
        
    fig.suptitle('Error Data for New York City, $' + var + '$', x=0.52, y=0.95, fontweight='semibold')
    note = '%s to %s, during the hours %s to %s UTC' % (start_date.strftime('%Y-%m-%d'), end_date.strftime('%Y-%m-%d'), (str(start_hr).zfill(2) + '00'), (str(end_hr).zfill(2) + '00'))
    ax.set_title(note, fontsize='medium')
    
    
    text_str = [err_names[i] + ': ' + "{:4.2f}".format(errs[i]) + ' \n' for i in range(0, len(errs)-1)]
    text_str.append(err_names[-1] + ': ' + "{:4d}".format(df.shape[0]) + '\n') 
    
    for i in range(0, len(text_str)):
        ax.annotate(text_str[i], xy=(0.02, 0.92-i/24), xycoords='axes fraction', fontsize='small')


def main(var, data=None, start_hr=0, end_hr=24, plot_type='heatmap'):
    
    if data is None:
        df = data_read('storage')
    else:
        df = data
    
    safe_data = df
    
    df['date'] = pd.to_datetime(df['date'], format='%Y-%m-%d %H:%M:%S')
    
    start_date, end_date = [datetime(year=2020, month=3, day=1), datetime(year=2020, month=5, day=31)]
    date_range = [start_date, end_date, start_hr, end_hr]
    
    df = df[(df['date'] > start_date) & (df['date'] < end_date)]    
    df = df[((df['date'].dt.hour >= start_hr) & (df['date'].dt.hour < end_hr))]
        
    var_meas = var + '_obs'
    
    date = df['date'].tolist()
    meas = df[var_meas].tolist()
    pred = df[var].tolist()
    
    errs = [None] * 5
    errs[0] = rmse(date, meas, pred)
    errs[1] = nsc(date, meas, pred)
    errs[2] = mbe(date, meas, pred)
    rgsn, errs[3] = r_sq(df, var)
    errs[4] = len(date)
    error_plot(var, meas, pred, errs, rgsn, date_range)
    
    