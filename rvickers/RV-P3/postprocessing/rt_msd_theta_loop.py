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
            x_avg = x.loc[idx[start:stop,'xcorr']].mean()
            y_avg = x.loc[idx[start:stop,'ycorr']].mean()
            thetax_avg = x.loc[idx[start:stop,'thetax']].mean()
            thetay_avg = x.loc[idx[start:stop,'thetay']].mean()
            thetaz_avg = x.loc[idx[start:stop,'thetaz']].mean()
            times.append([start, stop, stop-start,x_avg,y_avg,thetax_avg,thetay_avg,thetaz_avg])
    if not times:
        times=None
    return times

def get_theta_hist(x,rt_dist):
    cutDF = pd.cut(x.loc[idx[:,'z']],np.arange(0,rt_dist*100,rt_dist),labels=False).fillna(-1)
    deltaDF = cutDF[cutDF.diff() > 0]
    boundsDF = deltaDF[deltaDF.diff() != 0].index
    hist_x= []
    hist_y= []
    hist_z= []
    for i in range(len(boundsDF)-1): 
        if (cutDF[boundsDF[i]] < cutDF[boundsDF[i+1]]):
            start = boundsDF[i]-5000
            stop = boundsDF[i+1]
            thetax = x.loc[idx[start:stop,'thetax']]
            thetay = x.loc[idx[start:stop,'thetay']]
            thetaz = x.loc[idx[start:stop,'thetaz']]
            hist_x.append(thetax)
            hist_y.append(thetay)
            hist_z.append(thetaz)
    allxhist=np.concatenate(hist_x)
    allyhist=np.concatenate(hist_y)
    allzhist=np.concatenate(hist_z)
    return allxhist,allyhist,allzhist

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
stage = sys.argv[2]
rt_dist = int(sys.argv[3])
if gap_size == '03_90':
    inlet_z = -110
    outlet_z = 110
    z_min = -287
    z_max = 237
    offsetx=175.422
    offsety=175.422
    averaging_window = (71505000,1000000000)
elif gap_size == 'nacl03_90':
    inlet_z = -110
    outlet_z = 110
    z_min = -287
    z_max = 237
    offsetx=175.422
    offsety=175.422
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

parentDir = './'+stage+'_csv_mols_allt/'
for filenum in range(0,12):
    file = parentDir+'64x_'+gap_size+'_'+str(filenum)+'.csv'
    i=0
    rtDF_arr = []
    print('started: '+file)
    molDF = pd.read_csv(file,index_col=[0],header=[0,1])
    timesteps = [str(i) for i in range(averaging_window[0],averaging_window[1]+5000,5000) if str(i) in molDF.columns.get_level_values(0).to_list()]
    molDF = molDF[timesteps]
    molDF.rename(columns=lambda x: int(x), level=0, inplace=True)
    xcorrDF = (molDF.T.loc[idx[:,('x')],idx[:]]-offsetx*np.floor(molDF.T.loc[idx[:,('x')],idx[:]]/offsetx)).T.rename({'x':'xcorr'},axis=1)
    ycorrDF = (molDF.T.loc[idx[:,('y')],idx[:]]-offsety*np.floor(molDF.T.loc[idx[:,('y')],idx[:]]/offsety)).T.rename({'y':'ycorr'},axis=1)
    print('pre-concat')
    molDF = pd.concat([molDF,xcorrDF,ycorrDF],axis=1).sort_index(axis=1)
    print('post-concat')
    del ycorrDF
    gc.collect()
    print('post-del')
    molDF.iloc[:,molDF.columns.get_level_values(1)=='tube'] = molDF.iloc[:,molDF.columns.get_level_values(1)=='tube'].T.apply(col_count_consec).T
    molDF.iloc[:,molDF.columns.get_level_values(1)=='x'] = molDF.T.apply(get_dis,dim='x').T
    molDF.iloc[:,molDF.columns.get_level_values(1)=='y'] = molDF.T.apply(get_dis,dim='y').T
    molDF.iloc[:,molDF.columns.get_level_values(1)=='z'] = molDF.T.apply(get_dis,dim='z').T

    rtDF = molDF.T.apply(get_residence_time_and_y_diff,args=[int(rt_dist)]).dropna()
    rearr_rtDF = []
    for array1,index in zip(rtDF,rtDF.index):
        for array2 in array1:
            rearr_rtDF.append([index,array2[0],array2[1],array2[2],array2[3],array2[4],array2[5],array2[6],array2[7]])
    rearr_rtDF=pd.DataFrame(rearr_rtDF,columns=['mol','start','end','Residence time','avg x','avg y','avg thetax', 'avg thetay', 'avg thetaz']).set_index('mol')
    rearr_rtDF.to_csv('./'+stage+'_rt/64x_'+gap_size+'_rt_theta_'+str(filenum)+'.csv')
    print('finished: '+file)
