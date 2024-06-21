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

def col_group_consec(x):
    return  (x * (x != x.shift())).cumsum()

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
            start = int(boundsDF[i])-5000
            stop = int(boundsDF[i+1])
            y_avg = x.loc[idx[str(start):str(stop),'ycorr']].mean()
            times.append([start, stop, stop-start, y_avg])
    if not times:
        times=None
    return times

def get_MSD(x,dim):
    return x.iloc[:,x.columns.get_level_values(1).isin(('tube',dim))].droplevel(0,axis=1).groupby(['tube']).agg(['count','sum'])

def get_MSD_yavg(x,cl):
    return abs(x.xs('ycorr',level=1, drop_level=True).groupby(x.xs('tubeocc',level=1, drop_level=True)).transform('mean')-cl)

def get_MSD_yavg_out(x,cl):
    return abs(x.xs('ycorr',level=1, drop_level=True).groupby(x.xs('tubeocc',level=1, drop_level=True)).mean()-cl)


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

gap_size = sys.argv[1]
pressure_diffs = ['0','100','200','500','1000']
if gap_size == '05':
    averaging_window = (5000000,100000000)
    offsety = 149.97
    cuts = [[0,1.3],[1.3,100]]
elif gap_size == '06':
    averaging_window = (5000000,100000000)
    offsety = 151.47
    cuts = [[0,1.1],[1.1,100]]
elif gap_size == '10':
    averaging_window = (5000000,100000000)
    offsety = 156.357
    cuts = [[0,1.15],[1.15,100]]
elif gap_size == '20':
    averaging_window = (5000000,100000000)
    offsety = 166.129
    cuts = [[0,2.75],[2.75,100]]
elif gap_size == '30':
    averaging_window = (5000000,100000000)
    offsety = 175.901
    cuts = [[0,5],[5,100]]
cut_side = int(sys.argv[2])
cl = 70 + int(gap_size)/2
cl_lo = cl-1
cl_hi = cl+1

ang3_to_cm3 = (1e8)**3
avo_num = 6.0221409e+23

for pressure_diff in pressure_diffs:
    print(pressure_diff)
    molDF = pd.read_csv('./csv/pagap'+gap_size+'_nve'+pressure_diff+'_mols_u.csv',index_col=[0],header=[0,1])
    timesteps = [str(i) for i in range(averaging_window[0],averaging_window[1]+5000,5000) if str(i) in molDF.columns.get_level_values(0).to_list()]
    molDF = molDF[timesteps]
    ycorrDF = (molDF.T.loc[idx[:,('y')],idx[:]]-offsety*np.floor(molDF.T.loc[idx[:,('y')],idx[:]]/offsety)).T.rename({'y':'ycorr'},axis=1)
    molDF = pd.concat([molDF,ycorrDF],axis=1).sort_index(axis=1)
    tube_occDF = molDF.iloc[:,molDF.columns.get_level_values(1)=='tube'].T.apply(col_group_consec).T.rename({'tube':'tubeocc'},axis=1)
    molDF = pd.concat([molDF,tube_occDF],axis=1).sort_index(axis=1)
    molDF.iloc[:,molDF.columns.get_level_values(1)=='tube'] = molDF.iloc[:,molDF.columns.get_level_values(1)=='tube'].T.apply(col_count_consec).T
    yDF = molDF.T.apply(get_MSD_yavg,args=[cl]).T
    yDFout = molDF.T.apply(get_MSD_yavg_out,args=[cl]).T
    yDFout.to_csv('./csv_y0/pagap'+gap_size+'_nve'+pressure_diff+'_ydists.csv')
    yDF.columns = pd.MultiIndex.from_product([yDF.columns, ['y avg']])
    molDF = pd.concat([molDF,yDF],axis=1).sort_index(axis=1)
    
    tubeyDF = molDF.loc[:,idx[:,'tube']].droplevel(1,axis=1)[(molDF.loc[:,idx[:,'y avg']].droplevel(1,axis=1) >= cuts[cut_side][0]) & (molDF.loc[:,idx[:,'y avg']].droplevel(1,axis=1) < cuts[cut_side][1])].fillna(0)
    tubeyDF.columns = pd.MultiIndex.from_product([tubeyDF.columns, ['tube']])
    molDF.loc[:,idx[:,'tube']] = tubeyDF
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

    MSDDF.to_csv('./csv_y'+str(cut_side)+'/pagap'+gap_size+'_nve'+pressure_diff+'_molsMSDDF.csv')
    MSDDFz.to_csv('./csv_y'+str(cut_side)+'/pagap'+gap_size+'_nve'+pressure_diff+'_molsMSDDF_z.csv')
    MDDFz.to_csv('./csv_y'+str(cut_side)+'/pagap'+gap_size+'_nve'+pressure_diff+'_molsMDDF_z.csv')

    del molDF
    del SDDFz
    del SDDF
    del MSDDF
    del MSDDFz
    del MDDFz
    gc.collect()