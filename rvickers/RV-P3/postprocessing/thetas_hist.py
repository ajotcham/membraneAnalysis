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
    inlet_z = -110
    outlet_z = 110
    z_min = -287
    z_max = 237
    averaging_window = (71505000,1000000000)
    offsetx = 175.422
    offsety = 175.422
elif gap_size == 'nacl03_90':
    inlet_z = -110
    outlet_z = 110
    z_min = -287
    z_max = 237
    averaging_window = (505000,1000000000)
    offsetx = 175.422
    offsety = 175.422
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
parentDir = './'+stage+'_thetas/'
i=0
thetahist_arr = []
for file in os.listdir(parentDir):
    print('starting: '+file)
    # if file in ['64x_03_90_0.csv','64x_03_90_7.csv','64x_03_90_8.csv','64x_03_90_5.csv','64x_03_90_6.csv','64x_03_90_9.csv','64x_03_90_1.csv']:
    #     continue
    molDF = pd.read_csv(parentDir+file,index_col=[0],header=[0,1])
    timesteps = [str(i) for i in range(averaging_window[0],averaging_window[1]+5000,5000) if str(i) in molDF.columns.get_level_values(0).to_list()]
    molDF = molDF[timesteps]
    xcorrDF = (molDF.T.loc[idx[:,('x')],idx[:]]-offsetx*np.floor(molDF.T.loc[idx[:,('x')],idx[:]]/offsetx)).T.rename({'x':'xcorr'},axis=1)
    ycorrDF = (molDF.T.loc[idx[:,('y')],idx[:]]-offsety*np.floor(molDF.T.loc[idx[:,('y')],idx[:]]/offsety)).T.rename({'y':'ycorr'},axis=1)
    print('pre-concat')
    molDF = pd.concat([molDF,xcorrDF,ycorrDF],axis=1).sort_index(axis=1)
    del xcorrDF
    del ycorrDF
    gc.collect()
    alltube=np.array(molDF.T.droplevel(0).loc['tube'].values).flatten()
    allxcorr=np.array(molDF.T.droplevel(0).loc['xcorr'].values).flatten()
    allycorr=np.array(molDF.T.droplevel(0).loc['ycorr'].values).flatten()
    allz=np.array(molDF.T.droplevel(0).loc['z'].values).flatten()
    # allxys=np.array(molDF.T.droplevel(0).loc['thetaxy'].values).flatten()
    # allxzs=np.array(molDF.T.droplevel(0).loc['thetaxz'].values).flatten()
    # allyzs=np.array(molDF.T.droplevel(0).loc['thetayz'].values).flatten()
    del molDF
    gc.collect()
    allxcorr=allxcorr[alltube==True]
    allycorr=allycorr[alltube==True]
    allz=allz[alltube==True]
    # allxys=allxys[alltube==True]
    # allxzs=allxzs[alltube==True]
    # allyzs=allyzs[alltube==True]
    thetaDF = pd.DataFrame([allxcorr,allycorr,allz]).T
    thetaDF.columns = ['xcorr','ycorr','z']
    thetahist_arr.append(thetaDF[:])
    del allxcorr
    del allycorr
    del allz
    del thetaDF
    gc.collect()
    thetahist_tot = pd.concat(thetahist_arr)
    thetahist_tot.to_csv('./'+stage+'_intersection/64x_'+gap_size+'_intersection.csv')
    print('finished: '+file)