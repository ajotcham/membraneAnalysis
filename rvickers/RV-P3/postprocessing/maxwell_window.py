#!/usr/bin/env python
import re
import numpy as np 
import pandas as pd
import os
import gc
import gzip
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import sys
idx = pd.IndexSlice

def sorted_alphanumeric(data):
    convert = lambda text: int(text) if text.isdigit() else text.lower()
    alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ] 
    return sorted(data, key=alphanum_key)

def col_count_consec(x):
    return x * (x.groupby((x.shift() != x).cumsum()).cumcount())

def get_dis(x,dim='x'):
    length = len(x.loc[idx[:,'tube']].values)
    ind = np.array([i for i in range(length)] - x.loc[idx[:,'tube']].values).astype(int)
    return x.loc[idx[:,dim]].values - x.loc[idx[:,dim]].values[ind]

def get_residence_time(x,rt_dist):
    cutDF = pd.cut(x.loc[idx[:,'z']],np.arange(0,rt_dist*100,rt_dist),labels=False).fillna(-1)
    deltaDF = cutDF[cutDF.diff() > 0]
    boundsDF = deltaDF[deltaDF.diff() != 0].index
    times = [((int(boundsDF[i])-5000), int(boundsDF[i+1]), int(boundsDF[i+1])-(int(boundsDF[i])-5000)) for i in range(len(boundsDF)-1) if (cutDF[boundsDF[i]] < cutDF[boundsDF[i+1]])]
    if not times:
        times=None
    return times

def get_residence_time_and_y_diff(x,rt_dist):
    cutDF = pd.cut(x.loc[idx[:,'z']],np.arange(0,rt_dist*100,rt_dist),labels=False).fillna(-1)
    deltaDF = cutDF[cutDF.diff() > 0]
    boundsDF = deltaDF[deltaDF.diff() != 0].index
    times = []
    for i in range(len(boundsDF)-1): 
        if (cutDF[boundsDF[i]] < cutDF[boundsDF[i+1]]):
            start = boundsDF[i]-5000
            stop = boundsDF[i+1]
            y_avg = x.loc[idx[start:stop,'ycorr']].mean()
            x_dis = x.loc[idx[start:stop,'x']].diff().fillna(0).droplevel(1) ** 2
            y_dis = x.loc[idx[start:stop,'y']].diff().fillna(0).droplevel(1) ** 2
            tot_dist = np.sqrt((x_dis + y_dis)).sum()
            times.append([start, stop, stop-start, y_avg, tot_dist])
    if not times:
        times=None
    return times

def get_MSD(x,dim):
    return x.iloc[:,x.columns.get_level_values(1).isin(('tube',dim))].droplevel(0,axis=1).groupby(['tube']).agg(['count','sum'])

def get_maxwell_dist_ng(x):
    y = x.droplevel(0,axis=1)
    return y[y['tube']].groupby('v').count()

def get_maxwell_dist_out(x):
    y = x.droplevel(0,axis=1)
    return y[y['tube']==False].groupby('v').count()

x_min = 0
x_max = 97.7228
y_min = 0
y_max = 166.129
# x_min = 0
# x_max = 119.544
# y_min = 0
# y_max = 138.038
x_len = x_max - x_min
y_len = y_max - y_min


ang3_to_cm3 = (1e8)**3
avo_num = 6.0221409e+23

gap_min = 60
gap_max = 106

gap_len = gap_max - gap_min
nGapBins = int((gap_len * 0.5) + 1)
gapBinSize = gap_len/nGapBins

z_min = -128.921
z_max = 571.635
z_len = z_max - z_min

nBins = int((z_len * 0.5) + 1)
bin_size = z_len/nBins

min_timestep = 0
max_timestep = 800000000

inlet_z = 105
outlet_z = 280
rt_dist = '35'

gap_size = sys.argv[1]
pressure_diffs = ['0','100','200','500','1000']
if gap_size == '05':
    averaging_window = (5000000,100000000)
    offsety = 149.97
elif gap_size == '06':
    averaging_window = (5000000,100000000)
    offsety = 151.47
elif gap_size == '10':
    averaging_window = (5000000,100000000)
    offsety = 156.357
    cutoffs = [3.65,3.63]
elif gap_size == '20':
    averaging_window = (5000000,100000000)
    offsety = 166.129
elif gap_size == '30':
    averaging_window = (5000000,100000000)
    offsety = 175.901
elif gap_size == '50':
    averaging_window = (3000000,100000000)
    offsety = 195.446
    cutoffs = [14.57,14.37]
elif gap_size == '100':
    averaging_window = (3000000,100000000)
    offsety = 244.307
cl = 70 + int(gap_size)/2

cl_lo = cl-1
cl_hi = cl+1

ang3_to_cm3 = (1e8)**3
avo_num = 6.0221409e+23
windows = [10,100,500,1000]
for pressure_diff in pressure_diffs:
    molDF = pd.read_csv('./csv/pagap'+gap_size+'_nve'+pressure_diff+'_mols_v_theta.csv',index_col=[0],header=[0,1])
    timesteps = [str(i) for i in range(averaging_window[0],averaging_window[1]+5000,5000) if str(i) in molDF.columns.get_level_values(0).to_list()]
    molDF = molDF[timesteps]
    molDF.rename(columns=lambda x: int(x), level=0, inplace=True)
    vDF = (molDF.iloc[:,molDF.columns.get_level_values(1)=='vx'].rename({'vx':'v'}, axis=1) ** 2 + molDF.iloc[:,molDF.columns.get_level_values(1)=='vy'].rename({'vy':'v'}, axis=1) ** 2 + molDF.iloc[:,molDF.columns.get_level_values(1)=='vz'].rename({'vz':'v'}, axis=1) ** 2) ** 0.5
    tubeDF = molDF.iloc[:,molDF.columns.get_level_values(1)=='tube']
    ppDF = pd.concat([vDF,tubeDF], axis=1).sort_index(axis=1)
    trim_vDF=ppDF.iloc[:,ppDF.columns.get_level_values(1)=='v'].droplevel(1,axis=1)[ppDF.iloc[:,ppDF.columns.get_level_values(1)=='tube'].droplevel(1,axis=1)]
    arrays = []
    for window in windows:
        windowedtrim_vDF=trim_vDF.rolling(window,axis=1).mean()
        molecular_speeds=np.array(windowedtrim_vDF.values).flatten()
        arrays.append(molecular_speeds[~np.isnan(molecular_speeds)])
    windowedDF = pd.DataFrame(arrays).T
    windowedDF.columns=windows
    windowedDF.to_csv('./csv_maxwell_windowed/pagap'+gap_size+'_nve'+pressure_diff+'.csv')