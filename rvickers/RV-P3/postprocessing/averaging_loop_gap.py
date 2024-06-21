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
x_min = 0
x_max = 97.7228
y_min = 0
y_max = 244.307
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
pressure_diffs = ['0','100']

for gap_size in gap_sizes:
    print('Starting gap ' + gap_size)
    for pressure_diff in pressure_diffs:
        print('Starting pressure ' + pressure_diff)

        existing_density_gap_csv = "./csv/pagap"+gap_size+"_nve"+pressure_diff+"_density_gap.csv"
        existing_vx_gap_csv = "./csv/pagap"+gap_size+"_nve"+pressure_diff+"_vx_gap.csv"
        existing_vy_gap_csv = "./csv/pagap"+gap_size+"_nve"+pressure_diff+"_vy_gap.csv"
        existing_vz_gap_csv = "./csv/pagap"+gap_size+"_nve"+pressure_diff+"_vz_gap.csv"

        existingDensityGapDF = pd.read_csv(existing_density_gap_csv,index_col=0,dtype=np.float64)
        existingDensityGapDF.index = np.round(existingDensityGapDF.index,3)

        existingvxGapDF = pd.read_csv(existing_vx_gap_csv,index_col=0,dtype=np.float64)
        existingvxGapDF.index = np.round(existingvxGapDF.index,3)

        existingvyGapDF = pd.read_csv(existing_vy_gap_csv,index_col=0,dtype=np.float64)
        existingvyGapDF.index = np.round(existingvyGapDF.index,3)

        existingvzGapDF = pd.read_csv(existing_vz_gap_csv,index_col=0,dtype=np.float64)
        existingvzGapDF.index = np.round(existingvzGapDF.index,3)

        prevMaxtimestep = int(existingvxGapDF.columns[-1])
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
            
        if (allWaterArr):
            pass
        else:
            continue
        i=0
        min_atom_positions = []
        max_atom_positions = []
        # loop through arrays of file names containing data to be analyzed
        for waterArr in allWaterArr:
            timesteps = []
            waterdfs = []
#             print(waterArr[0])  
            i=0
            for waterFile in waterArr:
                #print(waterFile)
                if i > 1500:
                    continue
                try:
                    idf = pd.read_csv(waterFile)
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
                i+=1

        gapDensitySeries = []
        gapvxSeries = []
        gapvySeries = []
        gapvzSeries = []

        ## gap10
        # gap_min = 21
        # gap_max = 32

        ## gap46
        gap_min = 60
        gap_max = 200

        gap_len = gap_max - gap_min
        nGapBins = int((gap_len * 10) + 1)
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

        i=0
        for df in waterdfs:
            if ((int(timesteps[i]) < min_timestep) or (int(timesteps[i]) > max_timestep)):
                i+=1
                continue
            dfMol = pd.DataFrame()

            df['xmass'] = df['mass'] * df['x']
            df['ymass'] = df['mass'] * df['y']
            df['zmass'] = df['mass'] * df['z']

            df['vxmass'] = df['mass'] * df['vx']
            df['vymass'] = df['mass'] * df['vy']
            df['vzmass'] = df['mass'] * df['vz']

            # create per molecule dataframe to analyze at finer scale
            # atoms of a molecule are expected to be consecutive

            dfMol['mass'] = df['mass'].groupby(df.index // 3).sum()
            dfMol['xmass'] = df['xmass'].groupby(df.index // 3).sum()
            dfMol['ymass'] = df['ymass'].groupby(df.index // 3).sum()
            dfMol['zmass'] = df['zmass'].groupby(df.index // 3).sum()
            dfMol['vxmass'] = df['vxmass'].groupby(df.index // 3).sum()
            dfMol['vymass'] = df['vymass'].groupby(df.index // 3).sum()
            dfMol['vzmass'] = df['vzmass'].groupby(df.index // 3).sum()
            dfMol['x'] = dfMol['xmass'] / dfMol['mass']
            dfMol['y'] = dfMol['ymass'] / dfMol['mass']
            dfMol['z'] = dfMol['zmass'] / dfMol['mass']
            dfMol['vx'] = dfMol['vxmass'] / dfMol['mass']
            dfMol['vy'] = dfMol['vymass'] / dfMol['mass']
            dfMol['vz'] = dfMol['vzmass'] / dfMol['mass']
            dfMol['vol'] = df['c_peratomvol[1]'].groupby(df.index // 3).sum()
            dfMol['density'] = (dfMol['mass'] / (x_len * y_len * z_len/nBins)) * ang3_to_cm3 / avo_num
            dfMol['ybin'] = pd.cut(dfMol['y'].loc[(dfMol['z'] > inlet_z) & (dfMol['z'] < outlet_z)], np.linspace(gap_min,gap_max,nGapBins), labels=False) * gapBinSize + gap_min
            gapDF_mean = dfMol.groupby(['ybin'], dropna=True).mean()
            gapDF_sum = dfMol.groupby(['ybin'], dropna=True).sum()
            gapDF_sum['ybin_count'] = dfMol.groupby(['ybin'], dropna=True)['ybin'].count()

            gapDensitySeries.append((gapDF_sum['mass'].rename(timesteps[i])/(gapDF_sum['vol'].rename(timesteps[i])))* ang3_to_cm3 / avo_num)
            gapvxSeries.append(gapDF_mean['vx'].rename(timesteps[i]))
            gapvySeries.append(gapDF_mean['vy'].rename(timesteps[i]))
            gapvzSeries.append(gapDF_mean['vz'].rename(timesteps[i]))
            i += 1

        gapDensitySeriesDF = pd.concat(gapDensitySeries,axis=1, ignore_index=False)
        gapvxSeriesDF = pd.concat(gapvxSeries, axis=1, ignore_index=False)
        gapvySeriesDF = pd.concat(gapvySeries, axis=1, ignore_index=False)                    
        gapvzSeriesDF = pd.concat(gapvzSeries, axis=1, ignore_index=False)

        gapDensitySeriesDF.index = np.round(gapDensitySeriesDF.index, 3)
        gapvxSeriesDF.index = np.round(gapvxSeriesDF.index, 3)
        gapvySeriesDF.index = np.round(gapvySeriesDF.index, 3)
        gapvzSeriesDF.index = np.round(gapvzSeriesDF.index, 3)

        gapvxSeriesDF.index.names = ['Distance in y (Å)']
        gapvySeriesDF.index.names = ['Distance in y (Å)']
        gapvzSeriesDF.index.names = ['Distance in y (Å)']
        densityGapDF = pd.merge(existingDensityGapDF,gapDensitySeriesDF,how='outer',left_index=True, right_index=True)
        vxGapDF = pd.merge(existingvxGapDF,gapvxSeriesDF,how='outer',left_index=True, right_index=True)
        vyGapDF = pd.merge(existingvyGapDF,gapvySeriesDF,how='outer',left_index=True, right_index=True)
        vzGapDF = pd.merge(existingvzGapDF,gapvzSeriesDF,how='outer',left_index=True, right_index=True)
        AVGgapDensitySeriesDF = pd.DataFrame()

        AVGgapvxSeriesDF = pd.DataFrame()
        AVGgapvySeriesDF = pd.DataFrame()
        AVGgapvzSeriesDF = pd.DataFrame()            
        

        AVGgapDensitySeriesDF['mean'] = densityGapDF.mean(axis=1)
        AVGgapDensitySeriesDF['stdev'] = densityGapDF.std(axis=1)

        AVGgapvxSeriesDF['mean'] = vxGapDF.mean(axis=1)
        AVGgapvxSeriesDF['stdev'] = vxGapDF.std(axis=1)
        AVGgapvySeriesDF['mean'] = vyGapDF.mean(axis=1)
        AVGgapvySeriesDF['stdev'] = vyGapDF.std(axis=1)
        AVGgapvzSeriesDF['mean'] = vzGapDF.mean(axis=1)
        AVGgapvzSeriesDF['stdev'] = vzGapDF.std(axis=1)
        

        AVGgapDensitySeriesDF.index.names = ['Distance in y (Å)']  
        
        AVGgapvxSeriesDF.index.names = ['Distance in y (Å)']  
        AVGgapvySeriesDF.index.names = ['Distance in y (Å)']  
        AVGgapvzSeriesDF.index.names = ['Distance in y (Å)'] 
        densityGapDF.to_csv('./csv/pagap'+gap_size+'_nve'+pressure_diff+'_density_gap.csv')
        AVGgapDensitySeriesDF.to_csv('./csv/pagap'+gap_size+'_nve'+pressure_diff+'_density_gap_avg.csv')
        vxGapDF.to_csv('./csv/pagap'+gap_size+'_nve'+pressure_diff+'_vx_gap.csv')
        AVGgapvxSeriesDF.to_csv('./csv/pagap'+gap_size+'_nve'+pressure_diff+'_vx_gap_avg.csv')
        vyGapDF.to_csv('./csv/pagap'+gap_size+'_nve'+pressure_diff+'_vy_gap.csv')
        AVGgapvySeriesDF.to_csv('./csv/pagap'+gap_size+'_nve'+pressure_diff+'_vy_gap_avg.csv')
        vzGapDF.to_csv('./csv/pagap'+gap_size+'_nve'+pressure_diff+'_vz_gap.csv')
        AVGgapvzSeriesDF.to_csv('./csv/pagap'+gap_size+'_nve'+pressure_diff+'_vz_gap_avg.csv')
       
        