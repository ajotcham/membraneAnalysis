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

def get_xyz_rot_from_mat(rot_mat):
    thetaX = np.arctan2(-rot_mat[2][1],rot_mat[2][2])
    thetaY = np.arctan2(-rot_mat[2][0],(rot_mat[2][1] ** 2 + rot_mat[2][2] ** 2) ** 0.5)
    thetaZ = np.arctan2(-rot_mat[1][0],rot_mat[0][0])
    return(np.array([thetaX,thetaY,thetaZ]))

def get_mol_thetas(mol):
    ref_mol = np.array([[ 0., 0, 0.],
                    [ 0., -0.593649, -0.819341],
                    [ 0., -0.593649, 0.819341]])
    A=ref_mol.T
    B=mol.T
    centroid_A = np.mean(A, axis=1)
    centroid_B = np.mean(B, axis=1)
    centroid_A = centroid_A.reshape(-1, 1)
    centroid_B = centroid_B.reshape(-1, 1)
    Am = A - centroid_A
    Bm = B - centroid_B
    H = Am @ np.transpose(Bm)
    U, S, Vt = np.linalg.svd(H)
    R = Vt.T @ U.T
    if np.linalg.det(R) < 0:
        Vt[2,:] *= -1
        R = Vt.T @ U.T
    return get_xyz_rot_from_mat(R)*360/(2*np.pi)

def find_thetas(x):
    mol = np.array(x.values)
    return pd.Series(get_mol_thetas(mol),index=['thetax','thetay','thetaz'])


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
for pressure_diff in pressure_diffs:
    parentDir = "../64x/"+gap_size+"/"+stage+"/"
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
    # loop through arrays of file names containing data to be analyzed
    for waterArr in allWaterArr:
        timesteps = []
        waterdfs = []
        i=0
        for waterFile in waterArr:
            if i > 700:
                continue
            try:
                idf = pd.read_csv(waterFile)
            except:
                print('Skip because gz issue')
                continue
            timestep = idf.iloc[0].values[0]
            idf = idf.iloc[7:,:]
            dfCols = idf.iloc[0,].str.split(' ')[0]
            if 'xu' not in dfCols:
                continue
            print(waterFile)
            timesteps.append(timestep)
            del dfCols[0:2]
            df = idf.iloc[1:,:]['ITEM: TIMESTEP'].str.split(' ', expand=True)
            df.set_axis(dfCols,axis=1,inplace=True)
            df.reset_index(drop=True, inplace=True)
            df = df.apply(pd.to_numeric)
            df = df[df.type.isin([15,16])]
            df = df[df.vz != 0]
            waterdfs.append(df)
            i+=1
    print('pre dfMol create')
    i=0
    dfMolSeries = []
    for df in waterdfs:
        dfMol = pd.DataFrame()
        df['xmass'] = df['mass'] * df['xu']
        df['ymass'] = df['mass'] * df['yu']
        df['zmass'] = df['mass'] * df['zu']

        df['vxmass'] = df['mass'] * df['vx']
        df['vymass'] = df['mass'] * df['vy']
        df['vzmass'] = df['mass'] * df['vz']
        df = df.sort_values(by=['mol','id'])
        # create per molecule dataframe to analyze at finer scale
        # atoms of a molecule are expected to be consecutive
        grouped_df = df.groupby(df['mol'])
        dfMol[['thetax','thetay','thetaz']] = grouped_df[['x','y','z']].apply(find_thetas)
        dfMol['mass'] = grouped_df['mass'].sum()
        dfMol['xmass'] = grouped_df['xmass'].sum()
        dfMol['ymass'] = grouped_df['ymass'].sum()
        dfMol['zmass'] = grouped_df['zmass'].sum()

        dfMol['vxmass'] = grouped_df['vxmass'].sum()
        dfMol['vymass'] = grouped_df['vymass'].sum()
        dfMol['vzmass'] = grouped_df['vzmass'].sum()

        dfMol['x'] = dfMol['xmass'] / dfMol['mass']
        dfMol['y'] = dfMol['ymass'] / dfMol['mass']
        dfMol['z'] = dfMol['zmass'] / dfMol['mass']

        dfMol['vx'] = dfMol['vxmass'] / dfMol['mass']
        dfMol['vy'] = dfMol['vymass'] / dfMol['mass']
        dfMol['vz'] = dfMol['vzmass'] / dfMol['mass']

        dfMol['tube'] = (dfMol['z'] > inlet_z) & (dfMol['z'] < outlet_z)
        dfMolSeries.append(dfMol[['x','y','z','vx','vy','vz','thetax','thetay','thetaz','tube']])
        i+=1
    print('post dfMol create')
    molDF = pd.concat(dfMolSeries,axis=1,keys=timesteps)
    molDF.to_csv('./'+stage+'_csv_mols/64x_'+gap_size+'_mols_v_theta.csv')
    print('post moldf create and save')

    del molDF
    gc.collect()
