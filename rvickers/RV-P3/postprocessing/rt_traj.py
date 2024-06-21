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
            thetax_avg = x.loc[idx[start:stop,'thetax']].mean()
            thetay_avg = x.loc[idx[start:stop,'thetay']].mean()
            thetaz_avg = x.loc[idx[start:stop,'thetaz']].mean()
            times.append([start, stop, stop-start,y_avg,thetax_avg,thetay_avg,thetaz_avg])
    if not times:
        times=None
    return times

def get_MSD(x,dim):
    return x.iloc[:,x.columns.get_level_values(1).isin(('tube',dim))].droplevel(0,axis=1).groupby(['tube']).agg(['count','sum'])

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
rt_dist = sys.argv[2]

gap_size = sys.argv[1]
pressure_diffs = ['100','1000']
if gap_size == '05':
    averaging_window = (5000000,100000000)
    offsety = 149.97
elif gap_size == '06':
    averaging_window = (5000000,100000000)
    offsety = 151.47
elif gap_size == '10':
    averaging_window = (5000000,100000000)
    offsety = 156.357
elif gap_size == '20':
    averaging_window = (5000000,100000000)
    offsety = 166.129
elif gap_size == '30':
    averaging_window = (5000000,100000000)
    offsety = 175.901
elif gap_size == '50':
    averaging_window = (5000000,100000000)
    offsety = 195.446

cl = 70 + int(gap_size)/2
cl_lo = cl-1
cl_hi = cl+1

ang3_to_cm3 = (1e8)**3
avo_num = 6.0221409e+23

for pressure_diff in pressure_diffs:
    molDF = pd.read_csv('./csv/pagap'+gap_size+'_nve'+pressure_diff+'_mols_v_theta.csv',index_col=[0],header=[0,1])
    timesteps = [str(i) for i in range(averaging_window[0],averaging_window[1]+5000,5000) if str(i) in molDF.columns.get_level_values(0).to_list()]
    molDF = molDF[timesteps]
    molDF.rename(columns=lambda x: int(x), level=0, inplace=True)
    r_sample_df = pd.read_csv('./selected_samples/pa_gap{}_nve{}_r.csv'.format(gap_size,pressure_diff),index_col=[0])
    rt_sample_df = pd.read_csv('./selected_samples/pa_gap{}_nve{}_rt.csv'.format(gap_size,pressure_diff),index_col=[0])
    mols = []
    dfs = []
    for i,row in r_sample_df.iterrows():
      mol_id = int(row['mol'])
      mols.append(str(mol_id)+'_'+str(i))
      start_ts = int(row['start'])
      end_ts = int(row['end'])
      mol_entry = molDF.loc[mol_id]
      timestep_values = mol_entry.index.get_level_values(0)
      timestep_span = mol_entry[(timestep_values >= start_ts) & (timestep_values <= end_ts)]
      x_locs = timestep_span.index.get_level_values(1)=='x'
      y_locs = timestep_span.index.get_level_values(1)=='y'
      z_locs = timestep_span.index.get_level_values(1)=='z'
      x_thes = timestep_span.index.get_level_values(1)=='thetax'
      y_thes = timestep_span.index.get_level_values(1)=='thetay'
      z_thes = timestep_span.index.get_level_values(1)=='thetaz'
      sampleDF = pd.DataFrame()
      sampleDF['x'] = timestep_span[x_locs].values
      sampleDF['y'] = timestep_span[y_locs].values
      sampleDF['z'] = timestep_span[z_locs].values
      sampleDF['thetax'] = timestep_span[x_thes].values
      sampleDF['thetay'] = timestep_span[y_thes].values
      sampleDF['thetaz'] = timestep_span[z_thes].values
      dfs.append(sampleDF)
    trajectoryDF = pd.concat(dfs, axis=1, keys=mols)
    trajectoryDF.to_csv('./selected_trajs/pa_gap{}_nve{}_r.csv'.format(gap_size,pressure_diff))
    
    mols = []
    dfs = []    
    for i,row in rt_sample_df.iterrows():
      mol_id = int(row['mol'])
      mols.append(str(mol_id)+'_'+str(i))
      start_ts = int(row['start'])
      end_ts = int(row['end'])
      mol_entry = molDF.loc[mol_id]
      timestep_values = mol_entry.index.get_level_values(0)
      timestep_span = mol_entry[(timestep_values >= start_ts) & (timestep_values <= end_ts)]
      x_locs = timestep_span.index.get_level_values(1)=='x'
      y_locs = timestep_span.index.get_level_values(1)=='y'
      z_locs = timestep_span.index.get_level_values(1)=='z'
      x_thes = timestep_span.index.get_level_values(1)=='thetax'
      y_thes = timestep_span.index.get_level_values(1)=='thetay'
      z_thes = timestep_span.index.get_level_values(1)=='thetaz'
      sampleDF = pd.DataFrame()
      sampleDF['x'] = timestep_span[x_locs].values
      sampleDF['y'] = timestep_span[y_locs].values
      sampleDF['z'] = timestep_span[z_locs].values
      sampleDF['thetax'] = timestep_span[x_thes].values
      sampleDF['thetay'] = timestep_span[y_thes].values
      sampleDF['thetaz'] = timestep_span[z_thes].values
      dfs.append(sampleDF)
    trajectoryDF1 = pd.concat(dfs, axis=1, keys=mols)
    trajectoryDF1.to_csv('./selected_trajs/pa_gap{}_nve{}_rt.csv'.format(gap_size,pressure_diff))