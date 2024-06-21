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
for pressure_diff in pressure_diffs:
    existingMolDF = pd.read_csv('./csv/pagap'+gap_size+'_nve'+pressure_diff+'_mols.csv',index_col=[0],header=[0,1])
    prevMaxtimestep = int(existingMolDF.columns[-1][0])
    print('Using ' + str(prevMaxtimestep)+' as max previous timestep')
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
            if (int(step.split('.')[1]) > prevMaxtimestep):
                waterArr.append(dataDir+"/"+step)
        if waterArr:
            allWaterArr.append(waterArr)
        else:
            print("No water array created for: " + dataDir)
#     print(allWaterArr)

    if (allWaterArr):
        pass
    else:
        continue
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
            idf = pd.read_csv(waterFile).iloc[7:,:]
            dfCols = idf.iloc[0,].str.split(' ')[0]
            timesteps.append(timestep)
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
    combinedMolDF = pd.merge(existingMolDF, molDF, how='outer', left_index=True, right_index=True)
    combinedMolDF.to_csv('./csv/pagap'+gap_size+'_nve'+pressure_diff+'_mols.csv')

    del combinedMolDF
    del existingMolDF
    del molDF
    gc.collect()