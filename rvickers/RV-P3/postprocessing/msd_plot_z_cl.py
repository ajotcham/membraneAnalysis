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
pressure_diffs = ['0','100','200','500','1000']
if gap_size == '05':
    averaging_window = (5000000,100000000)
elif gap_size == '06':
    averaging_window = (5000000,100000000)
elif gap_size == '10':
    averaging_window = (5000000,100000000)
elif gap_size == '20':
    averaging_window = (5000000,100000000)
elif gap_size == '30':
    averaging_window = (5000000,100000000)

ylo = 0
yhi = 149.97
cl = 70 + int(gap_size)/2
cl_lo = cl-2.5
cl_hi = cl+2.5
offsety = 149.97


for pressure_diff in pressure_diffs:
    molDF = pd.read_csv('./csv/pagap'+gap_size+'_nve'+pressure_diff+'_mols_u.csv',index_col=[0],header=[0,1])
    timesteps = [str(i) for i in range(averaging_window[0],averaging_window[1]+5000,5000) if str(i) in molDF.columns.get_level_values(0).to_list()]
    molDF = molDF[timesteps]
    # put MSD into bins based on how long mol has been in the tube
    molDF.iloc[:,molDF.columns.get_level_values(1)=='tube'] = molDF.iloc[:,molDF.columns.get_level_values(1)=='tube'].T.apply(col_count_consec).T
    ycorrDF = (molDF.T.loc[idx[:,('y')],idx[:]]-offsety*np.floor(molDF.T.loc[idx[:,('y')],idx[:]]/offsety)).T.rename({'y':'ycorr'},axis=1)
    molDF = pd.concat([molDF,ycorrDF],axis=1).sort_index(axis=1)
    clDF = ((molDF.iloc[:,molDF.columns.get_level_values(1)=='ycorr'] > cl_lo) & (molDF.iloc[:,molDF.columns.get_level_values(1)=='ycorr'] < cl_hi)).rename({'ycorr':'cl'},axis=1)
    molDF = pd.concat([molDF,clDF],axis=1).sort_index(axis=1)
    molDF.iloc[:,molDF.columns.get_level_values(1)=='cl'] = molDF.iloc[:,molDF.columns.get_level_values(1)=='cl'].T.apply(col_count_consec).T
    molDF.iloc[:,molDF.columns.get_level_values(1)=='tube'] = molDF.T.loc[idx[:,('tube','cl')],idx[:]].groupby(level=0).min().T
    molDF.iloc[:,molDF.columns.get_level_values(1)=='x'] = molDF.T.apply(get_dis,dim='x').T
    molDF.iloc[:,molDF.columns.get_level_values(1)=='y'] = molDF.T.apply(get_dis,dim='y').T
    molDF.iloc[:,molDF.columns.get_level_values(1)=='z'] = molDF.T.apply(get_dis,dim='z').T
    
    SDDFz = molDF.iloc[:,molDF.columns.get_level_values(1)=='z'].rename({'z':'sdz'}, axis=1) ** 2
    SDDF = molDF.iloc[:,molDF.columns.get_level_values(1)=='x'].rename({'x':'sd'}, axis=1) ** 2 + molDF.iloc[:,molDF.columns.get_level_values(1)=='y'].rename({'y':'sd'}, axis=1) ** 2 + SDDFz.rename({'sdz':'sd'}, axis=1)
    molDF = pd.concat([molDF, SDDF, SDDFz], axis=1).sort_index(axis=1)
    
    MSDDF = molDF.groupby(level=0,axis=1).apply(get_MSD,'sd').droplevel([0,1], axis=1).groupby(level=0,axis=1).sum()
    MSDDFz = molDF.groupby(level=0,axis=1).apply(get_MSD,'sdz').droplevel([0,1], axis=1).groupby(level=0,axis=1).sum()
    MDDFz = molDF.groupby(level=0,axis=1).apply(get_MSD,'z').droplevel([0,1], axis=1).groupby(level=0,axis=1).sum()
    MSDDF['avg'] = MSDDF['sum']/MSDDF['count']
    MSDDFz['avg'] = MSDDFz['sum']/MSDDFz['count']
    MDDFz['avg'] = MDDFz['sum']/MDDFz['count']
    
    MSDDF.to_csv('./csv/pagap'+gap_size+'_nve'+pressure_diff+'_molsMSDDF_cl.csv')
    MSDDFz.to_csv('./csv/pagap'+gap_size+'_nve'+pressure_diff+'_molsMSDDF_z_cl.csv')
    MDDFz.to_csv('./csv/pagap'+gap_size+'_nve'+pressure_diff+'_molsMDDF_z_cl.csv')
    
    del molDF
    del SDDFz
    del SDDF
    del MSDDF
    del MSDDFz
    del MDDFz
    gc.collect()