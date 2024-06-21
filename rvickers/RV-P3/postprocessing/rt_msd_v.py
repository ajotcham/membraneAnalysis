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

def get_residence_time_and_y_diff(x,rt_dist,gap_size):
    cutDF = pd.cut(x.loc[idx[:,'z']],np.arange(0,rt_dist*100,rt_dist),labels=False).fillna(-1)
    deltaDF = cutDF[cutDF.diff() > 0]
    boundsDF = deltaDF[deltaDF.diff() != 0].index
    times = []
    for i in range(len(boundsDF)-1): 
        if (cutDF[boundsDF[i]] < cutDF[boundsDF[i+1]]):
            start = boundsDF[i]-5000
            stop = boundsDF[i+1]
            total = x.loc[idx[start:stop,'ycorrabs']].count()
            interim = x.loc[idx[start:stop,'ycorrabs']]
            outsi = interim[interim > float(gap_size)/2].count()
            y_avg = x.loc[idx[start:stop,'ycorrabs']].mean()
            x_dis = x.loc[idx[start:stop,'x']].diff().fillna(0).droplevel(1) ** 2
            y_dis = x.loc[idx[start:stop,'y']].diff().fillna(0).droplevel(1) ** 2
            tot_dist = np.sqrt((x_dis + y_dis)).sum()
            v_avg = x.loc[idx[start:stop,'v']].mean()
            times.append([start, stop, stop-start, y_avg, tot_dist, v_avg, outsi/total])
    if not times:
        times=None
    return times

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
rt_dist = int(sys.argv[3])
if gap_size == '03_90':
    inlet_z = -110
    outlet_z = 110
    z_min = -287
    z_max = 237
    offsety=175.422
    averaging_window = (71505000,1000000000)
elif gap_size == 'nacl03_90':
    inlet_z = -110
    outlet_z = 110
    z_min = -287
    z_max = 237
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

for pressure_diff in pressure_diffs:
    molDF = pd.read_csv('./csv/pagap'+gap_size+'_nve'+pressure_diff+'_mols_v_theta.csv',index_col=[0],header=[0,1])
    timesteps = [str(i) for i in range(averaging_window[0],averaging_window[1]+5000,5000) if str(i) in molDF.columns.get_level_values(0).to_list()]
    molDF = molDF[timesteps]
    molDF.rename(columns=lambda x: int(x), level=0, inplace=True)
    ycorrDF = (molDF.T.loc[idx[:,('y')],idx[:]]-offsety*np.floor(molDF.T.loc[idx[:,('y')],idx[:]]/offsety)).T.rename({'y':'ycorrabs'},axis=1)
    vDF = (molDF.iloc[:,molDF.columns.get_level_values(1)=='vx'].rename({'vx':'v'}, axis=1) ** 2 + molDF.iloc[:,molDF.columns.get_level_values(1)=='vy'].rename({'vy':'v'}, axis=1) ** 2 + molDF.iloc[:,molDF.columns.get_level_values(1)=='vz'].rename({'vz':'v'}, axis=1) ** 2) ** 0.5
    print('pre-concat')
    molDF = pd.concat([molDF,ycorrDF,vDF],axis=1).sort_index(axis=1)
    print('post-concat')
    del ycorrDF
    del vDF
    gc.collect()
    print('post-del')
    molDF.iloc[:,molDF.columns.get_level_values(1)=='tube'] = molDF.iloc[:,molDF.columns.get_level_values(1)=='tube'].T.apply(col_count_consec).T
    molDF.iloc[:,molDF.columns.get_level_values(1)=='x'] = molDF.T.apply(get_dis,dim='x').T
    molDF.iloc[:,molDF.columns.get_level_values(1)=='y'] = molDF.T.apply(get_dis,dim='y').T
    molDF.iloc[:,molDF.columns.get_level_values(1)=='z'] = molDF.T.apply(get_dis,dim='z').T

    rtDF = molDF.T.apply(get_residence_time_and_y_diff,args=[int(rt_dist),float(gap_size)]).dropna()
    rearr_rtDF = []
    for array1,index in zip(rtDF,rtDF.index):
        for array2 in array1:
            rearr_rtDF.append([index,array2[0],array2[1],array2[2],array2[3],array2[4],array2[5],array2[6]])
    del rtDF
    gc.collect()
    rearr_rtDF=pd.DataFrame(rearr_rtDF,columns=['mol','start','end','Residence time','avg absydiff', 'xy dist','avg v','ratio']).set_index('mol')
    # msds were looking weird when trimmed, see what happens without trimming of mols that don't go the designated distance at least one time.
    # molDF = molDF.loc[list(set(rearr_rtDF.index.values))] 
    if molDF.empty:
        print('Pressure {0} had no molecules move from designated start to finish, skipping.'.format(pressure_diff))
        continue
    else:
        print('Pressure {0} had some molecules crossed the finish, continuing'.format(pressure_diff))
    rearr_rtDF.to_csv('./'+stage+'_csv_rt/64x_'+gap_size+'_rt_v.csv')

    # SDDFz = molDF.iloc[:,molDF.columns.get_level_values(1)=='z'].rename({'z':'sdz'}, axis=1) ** 2
    # SDDF = molDF.iloc[:,molDF.columns.get_level_values(1)=='x'].rename({'x':'sd'}, axis=1) ** 2 + molDF.iloc[:,molDF.columns.get_level_values(1)=='y'].rename({'y':'sd'}, axis=1) ** 2 + SDDFz.rename({'sdz':'sd'}, axis=1)
    # molDF = pd.concat([molDF, SDDF, SDDFz], axis=1).sort_index(axis=1)
    
    # MSDDF = molDF.groupby(level=0,axis=1).apply(get_MSD,'sd').droplevel([0,1], axis=1).groupby(level=0,axis=1).sum()
    # MSDDFz = molDF.groupby(level=0,axis=1).apply(get_MSD,'sdz').droplevel([0,1], axis=1).groupby(level=0,axis=1).sum()
    # MDDFz = molDF.groupby(level=0,axis=1).apply(get_MSD,'z').droplevel([0,1], axis=1).groupby(level=0,axis=1).sum()
    # MSDDF['avg'] = MSDDF['sum']/MSDDF['count']
    # MSDDFz['avg'] = MSDDFz['sum']/MSDDFz['count']
    # MDDFz['avg'] = MDDFz['sum']/MDDFz['count']
    
    # MSDDF.to_csv('./csv'+rt_dist+'/pagap'+gap_size+'_nve'+pressure_diff+'_molsMSDDF.csv')
    # MSDDFz.to_csv('./csv'+rt_dist+'/pagap'+gap_size+'_nve'+pressure_diff+'_molsMSDDF_z.csv')
    # MDDFz.to_csv('./csv'+rt_dist+'/pagap'+gap_size+'_nve'+pressure_diff+'_molsMDDF_z.csv')
    
    del molDF
    del rearr_rtDF
    # del SDDFz
    # del SDDF
    # del MSDDF
    # del MSDDFz
    # del MDDFz
    gc.collect()