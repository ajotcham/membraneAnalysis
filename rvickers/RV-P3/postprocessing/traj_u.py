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
    molDF = pd.read_csv('./csv/pagap'+gap_size+'_nve'+pressure_diff+'_mols_u.csv',index_col=[0],header=[0,1])
    timesteps = [str(i) for i in range(averaging_window[0],averaging_window[1]+5000,5000) if str(i) in molDF.columns.get_level_values(0).to_list()]
    molDF = molDF[timesteps]
    for timestep in timesteps:
        molDF = molDF.loc[(molDF[timestep]['z'] < outlet_z) & (molDF[timestep]['z'] > inlet_z)]
    molDF.to_csv('./csv/pagap'+gap_size+'_nve'+pressure_diff+'_mols_u_trimmed.csv')

# xdisDF = molDF.iloc[:,molDF.columns.get_level_values(1)=='x'].T
# ydisDF = molDF.iloc[:,molDF.columns.get_level_values(1)=='y'].T
# zdisDF = molDF.iloc[:,molDF.columns.get_level_values(1)=='z'].T
# xdisDF = xdisDF.sub(xdisDF.iloc[0,:])
# ydisDF = ydisDF.sub(ydisDF.iloc[0,:])
# zdisDF = zdisDF.sub(zdisDF.iloc[0,:])
# disDF = pd.concat([xdisDF.T,ydisDF.T,zdisDF.T],axis=1)
# sqDisSeries = []
# for timestep in timesteps[1:]:
#     sqDisSeries.append((disDF[timestep]['x']**2 + disDF[timestep]['y']**2 + disDF[timestep]['z']**2).rename(timestep))
# ASDDF = pd.concat(sqDisSeries, axis=1)
# compare_fig, compare_ax=plt.subplots()
# ASDDF.mean().reset_index().plot(y=0,logy=True,logx=True,ax=compare_ax,label='Never Exited Gap')
# MSDDF.plot(y='avg',logy=True,logx=True,ax=compare_ax,label='Exited Gap')
# compare_ax.set_xlabel('Lag time within gap (ts)')
# compare_ax.set_ylabel('Mean Squared Displacement (Å$^2$)')
# compare_ax.set_title('MSD vs. lag time of water within 20 Å gap of polyamide')
# plt.savefig('./plots/pagap_msd_test.jpg',bbox_inches='tight')
