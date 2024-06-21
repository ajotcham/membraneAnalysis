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
pressure_diffs = ['500','1000']
# pressure_diffs = ['100','200','500','1000']
for pressure_diff in pressure_diffs:
    parentDir = "../pa_gap/pa_gap"+gap_size+"_nve"+pressure_diff+"/"
    allWaterArr = []
    # loop over all simulation folders
    for simDir in sorted_alphanumeric(os.listdir(parentDir)):
        waterArr = []
        # data is expected to be in a folder titled pressure
        if simDir != 'pressure':
            continue
        dataDir = parentDir+simDir
        for step in sorted_alphanumeric(os.listdir(dataDir)):
            waterArr.append(dataDir+"/"+step)
        if waterArr:
            allWaterArr.append(waterArr)
        else:
            print("No water array created for: " + dataDir)
#     print(allWaterArr)
    i=0
    # loop through arrays of file names containing data to be analyzed
    for waterArr in allWaterArr:
        timesteps = []
        waterdfs = []
        for waterFile in waterArr:
            print(waterFile)
            with gzip.open(waterFile) as f:
                f.readline().rstrip()
                timestep = f.readline().rstrip().decode("utf-8")
            if int(timestep) > 10000000:
                continue
            timesteps.append(timestep)
            idf = pd.read_csv(waterFile).iloc[7:,:]
            dfCols = idf.iloc[0,].str.split(' ')[0]
            del dfCols[0:2]
            df = idf.iloc[1:,:]['ITEM: TIMESTEP'].str.split(' ', expand=True)
            df.set_axis(dfCols,axis=1,inplace=True)
            df.reset_index(drop=True, inplace=True)
            df = df.apply(pd.to_numeric)
            df = df[df.type.isin([15,16])]
            df = df[df.vz != 0]
            waterdfs.append(df)
    print('pre dfMol create')
    i=0
    dfMolSeries = []
    for df in waterdfs:
        dfMol = pd.DataFrame()
        df['xmass'] = df['mass'] * df['x']
        df['ymass'] = df['mass'] * df['y']
        df['zmass'] = df['mass'] * df['z']

        # create per molecule dataframe to analyze at finer scale
        # atoms of a molecule are expected to be consecutive
        dfMol['mass'] = df['mass'].groupby(df['mol']).sum()
        dfMol['xmass'] = df['xmass'].groupby(df['mol']).sum()
        dfMol['ymass'] = df['ymass'].groupby(df['mol']).sum()
        dfMol['zmass'] = df['zmass'].groupby(df['mol']).sum()
        dfMol['x'] = dfMol['xmass'] / dfMol['mass']
        dfMol['y'] = dfMol['ymass'] / dfMol['mass']
        dfMol['z'] = dfMol['zmass'] / dfMol['mass']
        dfMol['tube'] = (dfMol['z'] > inlet_z) & (dfMol['z'] < outlet_z)
        dfMolSeries.append(dfMol[['x','y','z','tube']])
        i+=1
    print('post dfMol create')
    molDF = pd.concat(dfMolSeries,axis=1,keys=timesteps)
    molDF.to_csv('./csv/pagap'+gap_size+'_nve'+pressure_diff+'_mols.csv')
    print('post moldf create and save')
    # # put MSD into bins based on how long mol has been in the tube
    # xDF = molDF.iloc[:,molDF.columns.get_level_values(1)=='x'].T
    # yDF = molDF.iloc[:,molDF.columns.get_level_values(1)=='y'].T
    # zDF = molDF.iloc[:,molDF.columns.get_level_values(1)=='z'].T
    # tubeDF = molDF.iloc[:,molDF.columns.get_level_values(1)=='tube'].T

    # consecDF = tubeDF.apply(col_count_consec)
    # print('post consec create')
    # DF = pd.concat([xDF.T,yDF.T,zDF.T,consecDF.T],axis=1)
    # DF = DF.T
    # xDisDF = DF.apply(get_dis,dim='x')
    # yDisDF = DF.apply(get_dis,dim='y')
    # zDisDF = DF.apply(get_dis,dim='z')
    # print('post displace get')
    # #to replicate the other version, in which only the mols that stay in the gap the entire time are utilized:

    # # SDDF = (xDisDF**2 + yDisDF**2 + zDisDF**2)[consecDF.loc[:,(consecDF==197).any()].columns.tolist()]
    # # consecDF = pd.DataFrame(consecDF.values,columns=xDisDF.columns)[consecDF.loc[:,(consecDF==197).any()].columns.tolist()]
    # SDDF = (xDisDF**2 + yDisDF**2 + zDisDF**2)
    # consecDF = pd.DataFrame(consecDF.values,columns=xDisDF.columns)
    # MSD = np.zeros((len(SDDF.values),3))
    # for bins, values in zip(consecDF.values, SDDF.values):
    #     for binny, value in zip(bins,values):
    #         MSD[binny][0] += 1
    #         MSD[binny][1] += value
    # for binny in MSD:
    #     binny[2] = binny[1]/binny[0]

    # MSDDF = pd.DataFrame(MSD,columns=['count','sum','avg'])
    
    # # keep only those mols that stayed in the tube for the entire time, keeps N constant, but will be a problem 
    # # for larger openings and higher delta P

    # for timestep in timesteps:
    #     molDF = molDF.loc[(molDF[timestep]['z'] < outlet_z) & (molDF[timestep]['z'] > inlet_z)]

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


    del molDF
    # del MSDDF
    # del SDDF
    # del ASDDF
    # del xdisDF
    # del ydisDF
    # del zdisDF
    # del xDisDF
    # del yDisDF
    # del zDisDF
    # del xDF
    # del yDF
    # del zDF
    # del tubeDF
    # del consecDF
    # del DF
    # del disDF
    gc.collect()