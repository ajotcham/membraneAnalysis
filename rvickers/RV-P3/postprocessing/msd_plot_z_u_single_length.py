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
filenum = sys.argv[3]
length = sys.argv[4]

if gap_size == '03_90':
    inlet_z = -75
    outlet_z = 75
    z_min = -287
    z_max = 237
    averaging_window = (71505000,1000000000)
elif gap_size == 'nacl03_90':
    inlet_z = -75
    outlet_z = 75
    z_min = -287
    z_max = 237
    averaging_window = (505000,1000000000)
elif gap_size == '02_98_v3':
    inlet_z = 2
    outlet_z = 172
    z_min = -175
    z_max = 354.573
    averaging_window = (5000000,1000000000)

z_len = z_max - z_min

nBins = int((z_len * 0.5) + 1)
bin_size = z_len/nBins


# Set parent directory containing allsimulations to be analyzed
parentDir = './'+stage+'_csv_mols_allt_v2/'
file = parentDir+'64x_'+gap_size+'_'+str(filenum)+'.csv'
print('starting: '+file)
molDF = pd.read_csv(file,index_col=[0],header=[0,1])
timesteps = [str(i) for i in range(averaging_window[0],averaging_window[1]+5000,5000) if str(i) in molDF.columns.get_level_values(0).to_list()]
molDF = molDF[timesteps]
molDF.rename(columns=lambda x: int(x), level=0, inplace=True)
# put MSD into bins based on how long mol has been in the tube
molDF.iloc[:,molDF.columns.get_level_values(1)=='tube'] = molDF.iloc[:,molDF.columns.get_level_values(1)=='tube'].T.apply(col_count_consec).T
molDF.iloc[:,molDF.columns.get_level_values(1)=='x'] = molDF.T.apply(get_dis,dim='x').T
molDF.iloc[:,molDF.columns.get_level_values(1)=='y'] = molDF.T.apply(get_dis,dim='y').T
molDF.iloc[:,molDF.columns.get_level_values(1)=='z'] = molDF.T.apply(get_dis,dim='z').T
mask=(molDF.iloc[:,molDF.columns.get_level_values(1)=='tube'].T==molDF.iloc[:,molDF.columns.get_level_values(1)=='tube'].max(axis=1).values).T
mask2=molDF.iloc[:,molDF.columns.get_level_values(1)=='z'].droplevel(1,axis=1)[mask.droplevel(1,axis=1)].max(axis=1)
molDF = molDF[mask2.between(float(length),float(length)+25)]
np.savetxt('./'+stage+'_msd_samples/125_'+str(filenum)+'.csv',molDF.index.values)
SDDFz = molDF.iloc[:,molDF.columns.get_level_values(1)=='z'].rename({'z':'sdz'}, axis=1) ** 2
SDDF = molDF.iloc[:,molDF.columns.get_level_values(1)=='x'].rename({'x':'sd'}, axis=1) ** 2 + molDF.iloc[:,molDF.columns.get_level_values(1)=='y'].rename({'y':'sd'}, axis=1) ** 2 + SDDFz.rename({'sdz':'sd'}, axis=1)
molDF = pd.concat([molDF, SDDF, SDDFz], axis=1).sort_index(axis=1)

MSDDF = molDF.groupby(level=0,axis=1).apply(get_MSD,'sd').droplevel([0,1], axis=1).groupby(level=0,axis=1).sum()
MSDDFz = molDF.groupby(level=0,axis=1).apply(get_MSD,'sdz').droplevel([0,1], axis=1).groupby(level=0,axis=1).sum()
MDDFz = molDF.groupby(level=0,axis=1).apply(get_MSD,'z').droplevel([0,1], axis=1).groupby(level=0,axis=1).sum()

MSDDF['avg'] = MSDDF['sum']/MSDDF['count']
MSDDFz['avg'] = MSDDFz['sum']/MSDDFz['count']
MDDFz['avg'] = MDDFz['sum']/MDDFz['count']

MSDDF.to_csv('./'+stage+'_msd/64x_'+gap_size+'_molsMSDDF_'+str(length)+'_'+str(filenum)+'.csv')
MSDDFz.to_csv('./'+stage+'_msd/64x_'+gap_size+'_molsMSDDF_z_'+str(length)+'_'+str(filenum)+'.csv')
MDDFz.to_csv('./'+stage+'_msd/64x_'+gap_size+'_molsMDDF_z_'+str(length)+'_'+str(filenum)+'.csv')
