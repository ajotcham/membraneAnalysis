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

def get_MSD(x,dim):
    return x.iloc[:,x.columns.get_level_values(1).isin(('tube',dim))].droplevel(0,axis=1).groupby(['tube']).agg(['count','sum'])
x_min = 0.173378
x_max = 175.422
y_min = 0.173378
y_max = 175.422
# x_min = 0
# x_max = 119.544
# y_min = 0
# y_max = 138.038
x_len = x_max - x_min
y_len = y_max - y_min


ang3_to_cm3 = (1e8)**3
avo_num = 6.0221409e+23

gap_min = y_min
gap_max = y_max

gap_len = gap_max - gap_min
nGapBins = int((gap_len * 0.5) + 1)
gapBinSize = gap_len/nGapBins

gap_size = sys.argv[1]
pressure_diffs = ['100']
stage = sys.argv[2]

if gap_size == '03_90':
    inlet_z = -75
    outlet_z = 75
    z_min = -287
    z_max = 237
    averaging_window = (71505000,1000000000)
    indexes = ['0','1','2','3','4','5']
    molRanges = [[0,40000],[40000,80000],[80000,120000],[120000,160000],[160000,200000],[200000,240000]]
elif gap_size == 'nacl03_90':
    inlet_z = -75
    outlet_z = 75
    z_min = -287
    z_max = 237
    averaging_window = (505000,1000000000)
    indexes = ['0','1','2','3','4']
    molRanges = [[0,40000],[40000,80000],[80000,120000],[120000,160000],[160000,200000],[200000,240000]]
elif gap_size == '02_98_v3':
    inlet_z = 2
    outlet_z = 172
    z_min = -175
    z_max = 354.573
    averaging_window = (5000000,1000000000)

z_len = z_max - z_min

nBins = int((z_len * 0.5) + 1)
bin_size = z_len/nBins

for j,molRange in enumerate(molRanges):
  beg = molRange[0]
  truncMolDFs = []
  for indy in indexes:
    molDF = pd.read_csv('./'+stage+'_csv_mols_v2/64x_'+gap_size+'_mols_v_theta_'+indy+'.csv',index_col=[0],header=[0,1])
    end = min(molRange[1],len(molDF))
    truncMolDF = molDF[beg:end]
    truncMolDFs.append(truncMolDF[:])
    del molDF
    del truncMolDF
    gc.collect()
  truncMolDF_allT = pd.concat(truncMolDFs,axis=1)
  truncMolDF_allT.to_csv('./'+stage+'_csv_mols_allt_v2/64x_'+gap_size+'_'+str(j)+'.csv')
  del truncMolDFs
  del truncMolDF_allT
  gc.collect()