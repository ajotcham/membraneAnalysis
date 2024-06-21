#!/usr/bin/env python
print('Starting analysis')
import re
def sorted_alphanumeric(data):
    convert = lambda text: int(text) if text.isdigit() else text.lower()
    alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ] 
    return sorted(data, key=alphanum_key)
import numpy as np 
import pandas as pd
import os
import gc
import gzip
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import sys
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

# Set parent directory containing allsimulations to be analyzed 

gap_sizes = [sys.argv[1]]
pressure_diffs = ['100']
stage = sys.argv[2]

for gap_size in gap_sizes:
    print('Starting gap ' + gap_size)
    for pressure_diff in pressure_diffs:
        print('Starting pressure ' + pressure_diff)
        parentDir = "../64x/"+gap_size+"/"+stage+"/"
        allWaterArr = []
        allMembrArr = []
        # loop over all simulation folders
        for simDir in sorted_alphanumeric(os.listdir(parentDir)):
            # data is expected to be in a folder titled pressure
            if simDir != 'pressure':
                continue
            waterArr = []
            dataDir = parentDir+simDir
            for step in sorted_alphanumeric(os.listdir(dataDir)):
                waterArr.append(dataDir+"/"+step)
            if waterArr:
                allWaterArr.append(waterArr)
            else:
                print("No water array created for: " + dataDir)
        
        for simDir in sorted_alphanumeric(os.listdir(parentDir)):
            
            # data is expected to be in a folder titled pressure
            if simDir != 'membrane':
                continue
            membrArr = []
            dataDir = parentDir+simDir
            for step in sorted_alphanumeric(os.listdir(dataDir)):
                membrArr.append(dataDir+"/"+step)
            if membrArr:
                allMembrArr.append(membrArr)
            else:
                print("No water array created for: " + dataDir)
        
        min_atom_positions = []
        max_atom_positions = []
        # loop through arrays of file names containing data to be analyzed
        for waterArr in allWaterArr:
            timesteps = []
            waterdfs = []
            membrdfs = []
#             print(waterArr[0])  
        i=0
        print(waterArr)
        print(membrArr)
        for waterFile,membrFile in zip(waterArr,membrArr):
            if i > 1000:
                continue
            #print(waterFile)
            try:
                idf = pd.read_csv(waterFile)
                jdf = pd.read_csv(membrFile)
            except:
                print('Skip because gz issue')
                continue
            timestep = idf.iloc[0].values[0]
            print(timestep)
            timesteps.append(timestep)
            idf = idf.iloc[7:,:]
            dfCols = idf.iloc[0,].str.split(' ')[0]
            del dfCols[0:2]
            df = idf.iloc[1:,:]['ITEM: TIMESTEP'].str.split(' ', expand=True)
            df.set_axis(dfCols,axis=1,inplace=True)
            df.reset_index(drop=True, inplace=True)
            df = df.apply(pd.to_numeric)
            min_atom_positions.append(df['z'].min())
            max_atom_positions.append(df['z'].max())
            waterdfs.append(df)

            jdf = jdf.iloc[7:,:]
            dfCols = jdf.iloc[0,].str.split(' ')[0]
            del dfCols[0:2]
            df = jdf.iloc[1:,:]['ITEM: TIMESTEP'].str.split(' ', expand=True)
            df.set_axis(dfCols,axis=1,inplace=True)
            df.reset_index(drop=True, inplace=True)
            df = df.apply(pd.to_numeric)
            min_atom_positions.append(df['z'].min())
            max_atom_positions.append(df['z'].max())
            membrdfs.append(df)
            i+=1
        wpressureSeries = []
        wmassSeries = []
        wvolumeSeries = []

        mpressureSeries = []
        mmassSeries = []
        mvolumeSeries = []

        ## gap10
        # gap_min = 21
        # gap_max = 32

        ## gap46
        gap_min = y_min
        gap_max = y_max

        gap_len = gap_max - gap_min
        nGapBins = int((gap_len * 0.5) + 1)
        gapBinSize = gap_len/nGapBins

        if gap_size == '03_90':
            inlet_z = -110
            outlet_z = 110
            z_min = -287
            z_max = 237
        elif gap_size == 'nacl03_90':
            inlet_z = -110
            outlet_z = 110
            z_min = -287
            z_max = 237
        elif gap_size == '02_98_v3':
            inlet_z = 2
            outlet_z = 172
            z_min = -175
            z_max = 354.573

        z_len = z_max - z_min

        nBins = int((z_len * 0.5) + 1)
        bin_size = z_len/nBins

        min_timestep = 0
        max_timestep = 800000000

        i=0
        print(len(waterdfs))
        print(len(membrdfs))
        for wdf,mdf in zip(waterdfs,membrdfs):
            if ((int(timesteps[i]) < min_timestep) or (int(timesteps[i]) > max_timestep)):
                i+=1
                continue
            wdf['bin'] = pd.cut(wdf['z'], np.linspace(z_min,z_max,nBins),labels=False)*bin_size + (z_min)
            mdf['bin'] = pd.cut(mdf['z'], np.linspace(z_min,z_max,nBins),labels=False)*bin_size + (z_min)
            stepwDF_sum = wdf.groupby(['bin'], dropna=True).sum()
            stepmDF_sum = mdf.groupby(['bin'], dropna=True).sum()

            wpressureSeries.append(stepwDF_sum['v_peratompress'].rename(timesteps[i]))
            wmassSeries.append(stepwDF_sum['mass'].rename(timesteps[i]))
            wvolumeSeries.append(stepwDF_sum['c_peratomvol[1]'].rename(timesteps[i]))

            mpressureSeries.append(stepmDF_sum['v_peratompress'].rename(timesteps[i]))
            mmassSeries.append(stepmDF_sum['mass'].rename(timesteps[i]))
            mvolumeSeries.append(stepmDF_sum['c_peratomvol[1]'].rename(timesteps[i]))

            i += 1

        wpressureSeriesDF = pd.concat(wpressureSeries, axis=1, ignore_index=False)
        wmassSeriesDF = pd.concat(wmassSeries,axis=1, ignore_index=False)
        wvolumeSeriesDF = pd.concat(wvolumeSeries, axis=1, ignore_index=False)

        mpressureSeriesDF = pd.concat(mpressureSeries, axis=1, ignore_index=False)
        mmassSeriesDF = pd.concat(mmassSeries,axis=1, ignore_index=False)
        mvolumeSeriesDF = pd.concat(mvolumeSeries, axis=1, ignore_index=False)


        wpressureSeriesDF.index = np.round(wpressureSeriesDF.index, 3)
        wmassSeriesDF.index = np.round(wmassSeriesDF.index, 3)
        wvolumeSeriesDF.index = np.round(wvolumeSeriesDF.index, 3)

        mpressureSeriesDF.index = np.round(mpressureSeriesDF.index, 3)
        mmassSeriesDF.index = np.round(mmassSeriesDF.index, 3)
        mvolumeSeriesDF.index = np.round(mvolumeSeriesDF.index, 3)

        wpressureSeriesDF.index.names = ['Distance in z (Å)']
        wmassSeriesDF.index.names = ['Distance in z (Å)']
        wvolumeSeriesDF.index.names = ['Distance in z (Å)']

        mpressureSeriesDF.index.names = ['Distance in z (Å)']
        mmassSeriesDF.index.names = ['Distance in z (Å)']
        mvolumeSeriesDF.index.names = ['Distance in z (Å)']

        AVGwpressureSeriesDF = pd.DataFrame()
        AVGwmassSeriesDF = pd.DataFrame()
        AVGwvolumeSeriesDF = pd.DataFrame()

        AVGmpressureSeriesDF = pd.DataFrame()
        AVGmmassSeriesDF = pd.DataFrame()
        AVGmvolumeSeriesDF = pd.DataFrame()

        AVGwpressureSeriesDF['mean'] = wpressureSeriesDF.mean(axis=1)
        AVGwpressureSeriesDF['stdev'] = wpressureSeriesDF.std(axis=1)

        AVGwmassSeriesDF['mean'] = wmassSeriesDF.mean(axis=1)
        AVGwmassSeriesDF['stdev'] = wmassSeriesDF.std(axis=1)

        AVGwvolumeSeriesDF['mean'] = wvolumeSeriesDF.mean(axis=1)
        AVGwvolumeSeriesDF['stdev'] = wvolumeSeriesDF.std(axis=1)

        AVGmpressureSeriesDF['mean'] = mpressureSeriesDF.mean(axis=1)
        AVGmpressureSeriesDF['stdev'] = mpressureSeriesDF.std(axis=1)

        AVGmmassSeriesDF['mean'] = mmassSeriesDF.mean(axis=1)
        AVGmmassSeriesDF['stdev'] = mmassSeriesDF.std(axis=1)

        AVGmvolumeSeriesDF['mean'] = mvolumeSeriesDF.mean(axis=1)
        AVGmvolumeSeriesDF['stdev'] = mvolumeSeriesDF.std(axis=1)

        AVGwpressureSeriesDF.index.names = ['Distance in z (Å)']
        AVGwmassSeriesDF.index.names = ['Distance in z (Å)']
        AVGwvolumeSeriesDF.index.names = ['Distance in z (Å)']

        AVGmpressureSeriesDF.index.names = ['Distance in z (Å)']
        AVGmmassSeriesDF.index.names = ['Distance in z (Å)']
        AVGmvolumeSeriesDF.index.names = ['Distance in z (Å)']  

        wpressureSeriesDF.to_csv('./'+stage+'_csv_all/64x_'+gap_size+'_pressure_w.csv')
        AVGwpressureSeriesDF.to_csv('./'+stage+'_csv_all/64x_'+gap_size+'_pressure_w_avg.csv')

        wmassSeriesDF.to_csv('./'+stage+'_csv_all/64x_'+gap_size+'_mass_w.csv')
        AVGwmassSeriesDF.to_csv('./'+stage+'_csv_all/64x_'+gap_size+'_mass_w_avg.csv')

        wvolumeSeriesDF.to_csv('./'+stage+'_csv_all/64x_'+gap_size+'_volume_w.csv')
        AVGwvolumeSeriesDF.to_csv('./'+stage+'_csv_all/64x_'+gap_size+'_volume_w_avg.csv')

        mpressureSeriesDF.to_csv('./'+stage+'_csv_all/64x_'+gap_size+'_pressure_m.csv')
        AVGmpressureSeriesDF.to_csv('./'+stage+'_csv_all/64x_'+gap_size+'_pressure_m_avg.csv')

        mmassSeriesDF.to_csv('./'+stage+'_csv_all/64x_'+gap_size+'_mass_m.csv')
        AVGmmassSeriesDF.to_csv('./'+stage+'_csv_all/64x_'+gap_size+'_mass_m_avg.csv')

        mvolumeSeriesDF.to_csv('./'+stage+'_csv_all/64x_'+gap_size+'_volume_m.csv')
        AVGmvolumeSeriesDF.to_csv('./'+stage+'_csv_all/64x_'+gap_size+'_volume_m_avg.csv')

