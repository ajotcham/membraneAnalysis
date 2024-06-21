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

# Set parent directory containing allsimulations to be analyzed 

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
    cutoffs = [3.65,3.63]
elif gap_size == '20':
    averaging_window = (5000000,100000000)
    offsety = 166.129
elif gap_size == '30':
    averaging_window = (5000000,100000000)
    offsety = 175.901
elif gap_size == '50':
    averaging_window = (5000000,100000000)
    offsety = 195.446
    cutoffs = [14.57,14.37]
elif gap_size == '100':
    averaging_window = (3000000,100000000)
    offsety = 244.307
cl = 70 + int(gap_size)/2

for pressure_diff,cutoff in zip(pressure_diffs,cutoffs):
    molDF = pd.read_csv('./csv_intersect/pagap'+gap_size+'_nve'+pressure_diff+'_mols.csv',index_col=[0],header=[0,1])
    timesteps = [str(i) for i in range(averaging_window[0],averaging_window[1]+5000,5000) if str(i) in molDF.columns.get_level_values(0).to_list()]
    molDF = molDF[timesteps]
    molDF.rename(columns=lambda x: int(x), level=0, inplace=True)
    r_sample_df = pd.read_csv('./random_samples/pa_gap{}_nve{}_dy.csv'.format(gap_size,pressure_diff),index_col=[0])
    rt_sample_df = pd.read_csv('./random_samples/pa_gap{}_nve{}_rt.csv'.format(gap_size,pressure_diff),index_col=[0])
    ycorrDF = (molDF.T.loc[idx[:,('y')],idx[:]]-offsety*np.floor(molDF.T.loc[idx[:,('y')],idx[:]]/offsety)).T.rename({'y':'ycorr'},axis=1)
    print('pre-concat')
    ycorrDF = abs(ycorrDF-cl)
    molDF = pd.concat([molDF,ycorrDF],axis=1).sort_index(axis=1)
    thetaDFarr=[]
    mol_ids=[]
    for dffy in [r_sample_df, rt_sample_df]:
      for i, row in dffy.iterrows():
        mol_id = int(row['mol'])
        mol_ids.append(mol_id)
        mol_entry = molDF.loc[mol_id]
        alltube=np.array(mol_entry.T.droplevel(0).loc['tube'].values).flatten()
        allycorr=np.array(mol_entry.T.droplevel(0).loc['ycorr'].values).flatten()
        allxys=np.array(mol_entry.T.droplevel(0).loc['thetaxy'].values).flatten()
        allxzs=np.array(mol_entry.T.droplevel(0).loc['thetaxz'].values).flatten()
        allyzs=np.array(mol_entry.T.droplevel(0).loc['thetayz'].values).flatten()
        allycorr=allycorr[alltube==True]
        allxys=allxys[alltube==True]
        allxzs=allxzs[alltube==True]
        allyzs=allyzs[alltube==True]
        thetaDF = pd.DataFrame([allycorr,allxys,allxzs,allyzs]).T
        thetaDF.columns = ['ycorr','thetaxy','thetaxz','thetayz']
        thetaDFarr.append(thetaDF)
    thetaDFall= pd.concat(thetaDFarr,axis=1,keys=mol_ids)
    thetaDFall.to_csv('./csv_intersect_clean/pagap'+gap_size+'_nve'+pressure_diff+'_sample2.csv')