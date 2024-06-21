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
    averaging_window = (5000000,20000000)
elif gap_size == '06':
    averaging_window = (5000000,20000000)
elif gap_size == '10':
    averaging_window = (5000000,20000000)
elif gap_size == '20':
    averaging_window = (5000000,20000000)
elif gap_size == '30':
    averaging_window = (5000000,20000000)


for pressure_diff in pressure_diffs:
    molDF = pd.read_csv('./csv/pagap'+gap_size+'_nve'+pressure_diff+'_mols.csv',index_col=[0],header=[0,1])
    timesteps = [str(i) for i in range(averaging_window[0],averaging_window[1]+5000,5000) if str(i) in molDF.columns.get_level_values(0).to_list()]
    steps_list = np.array_split(timesteps,int(sys.argv[2]))
    molDFS = [molDF[steps.tolist()] for steps in steps_list]
    i=0
    for molDF in molDFS:
        # put MSD into bins based on how long mol has been in the tube
        zDF = molDF.iloc[:,molDF.columns.get_level_values(1)=='z'].T
        tubeDF = molDF.iloc[:,molDF.columns.get_level_values(1)=='tube'].T

        consecDF = tubeDF.apply(col_count_consec)
        print('post consec create')
        DF = pd.concat([zDF.T,consecDF.T],axis=1)
        DF = DF.T
        # xDisDF = DF.apply(get_dis,dim='x')
        # yDisDF = DF.apply(get_dis,dim='y')
        zDisDF = DF.apply(get_dis,dim='z')
        print('post displace get')
        #to replicate the other version, in which only the mols that stay in the gap the entire time are utilized:

        # SDDF = (xDisDF**2 + yDisDF**2 + zDisDF**2)[consecDF.loc[:,(consecDF==197).any()].columns.tolist()]
        # consecDF = pd.DataFrame(consecDF.values,columns=xDisDF.columns)[consecDF.loc[:,(consecDF==197).any()].columns.tolist()]
        SDDF = zDisDF ** 2
        # SDDF.to_csv('./csv/pagap'+gap_size+'_nve'+pressure_diff+'_molsSDDF_z.csv')
        consecDF = pd.DataFrame(consecDF.values,columns=zDisDF.columns)
        # consecDF.to_csv('./csv/pagap'+gap_size+'_nve'+pressure_diff+'_molsconsec.csv')
        MSD = np.zeros((len(SDDF.values),3))
        for bins, values in zip(consecDF.values, SDDF.values):
            for binny, value in zip(bins,values):
                MSD[binny][0] += 1
                MSD[binny][1] += value
        for binny in MSD:
            binny[2] = binny[1]/binny[0]

        MSDDF = pd.DataFrame(MSD,columns=['count','sum','avg'])
        MSDDF.to_csv('./csv/pagap'+gap_size+'_nve'+pressure_diff+'_molsMSDDF_z_'+sys.argv[2]+'_'+str(i)+'.csv')
        
        # keep only those mols that stayed in the tube for the entire time, keeps N constant, but will be a problem 
        # for larger openings and higher delta P
        del tubeDF
        del consecDF
        del DF
        del zDF
        del zDisDF
        del SDDF
        del MSD
        del MSDDF
        gc.collect()
        i+=1
