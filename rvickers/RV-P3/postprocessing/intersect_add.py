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
    thetaX = np.arctan2(rot_mat[2][1],rot_mat[2][2])
    thetaY = np.arctan2(-rot_mat[2][0],(rot_mat[2][1] ** 2 + rot_mat[2][2] ** 2) ** 0.5)
    thetaZ = np.arctan2(rot_mat[1][0],rot_mat[0][0])
    return(np.array([thetaX,thetaY,thetaZ]))
# def get_mol_thetas(mol):
# #     print(mol)
#     # ref_mol = np.array([[ 0., 0, 0.],
#     #                 [ 0., -0.593649, -0.819341],
#     #                 [ 0., -0.593649, 0.819341]])
#     ref_mol = np.array([[ 0., 0, 0],
#                         [  -0.819341 , -0.593649,0.],
#                         [  0.819341, -0.593649,0.],
#                         ])
#     A=ref_mol.T
#     B=mol.T
    
#     centroid_A = np.mean(A, axis=1)
#     centroid_B = np.mean(B, axis=1)
#     centroid_A = centroid_A.reshape(-1, 1)
#     centroid_B = centroid_B.reshape(-1, 1)
#     Am = A - centroid_A
#     Bm = B - centroid_B
# #     print(Am)
# #     print(Bm)
#     H = Am @ np.transpose(Bm)
#     U, S, Vt = np.linalg.svd(H)
#     R = Vt.T @ U.T
# #     print(R)
#     if np.linalg.det(R) < 0:
#         Vt[2,:] *= -1
#         R = Vt.T @ U.T
#     r=RTRANS.from_matrix(R)
#     return r.as_euler('zyx', degrees=True)

# def get_mol_thetas1(mol):
# #     print(mol)
#     ref_mol = np.array([[ 0., 0, 0.],
#                     [ 0., -0.593649, -0.819341],
#                     [ 0., -0.593649, 0.819341]])
# #     ref_mol = np.array([[ 0., 0, 0],
# #                         [  -0.819341 , -0.593649,0.],
# #                         [  0.819341, -0.593649,0.],
# #                         ])
#     A=ref_mol.T
#     B=mol.T
    
#     centroid_A = np.mean(A, axis=1)
#     centroid_B = np.mean(B, axis=1)
#     centroid_A = centroid_A.reshape(-1, 1)
#     centroid_B = centroid_B.reshape(-1, 1)
#     Am = A - centroid_A
#     Bm = B - centroid_B
# #     print(Am)
# #     print(Bm)
#     H = Am @ np.transpose(Bm)
#     U, S, Vt = np.linalg.svd(H)
#     R = Vt.T @ U.T
# #     print(R)
#     if np.linalg.det(R) < 0:
#         Vt[2,:] *= -1
#         R = Vt.T @ U.T
#     r=RTRANS.from_matrix(R)
#     return r.as_euler('zyx', degrees=True)

def find_angle(vec1,vec2):
    numer = np.dot(vec1,vec2)
    denumer = np.linalg.norm(vec1) * np.linalg.norm(vec2)
    rad_ang = np.arccos(numer/denumer)
    deg_ang = rad_ang*360/(2*np.pi)
#     return min(deg_ang,180-deg_ang)
#     return deg_ang
    return deg_ang
def get_mol_thetas_planes(mol):
    xy_norm = np.array([0,0,1])
    xz_norm = np.array([0,1,0])
    yz_norm = np.array([1,0,0])
    B=mol.T
    centroid_B = np.mean(B, axis=1)
    centroid_B = centroid_B.reshape(-1, 1)
    Bm = B - centroid_B
#     print(Bm)
#     print(Bm[:,1]-Bm[:,0])
#     print(Bm[:,2]-Bm[:,0])
    Bm_norm=np.cross(Bm[:,2]-Bm[:,0],Bm[:,1]-Bm[:,0])
    theta_xy = find_angle(Bm_norm,xy_norm)
    theta_xz = find_angle(Bm_norm,xz_norm)
    theta_yz = find_angle(Bm_norm,yz_norm)
    return np.array([theta_xy,theta_xz,theta_yz])
    
# def find_thetas(x):
#     mol = np.array(x.values)
#     return pd.Series(get_mol_thetas(mol),index=['thetax','thetay','thetaz'])

# def find_thetas1(x):
#     mol = np.array(x.values)
#     return pd.Series(get_mol_thetas1(mol),index=['thetax1','thetay1','thetaz1'])

def find_thetas_planes(x):
    mol = np.array(x.values)
    return pd.Series(get_mol_thetas_planes(mol),index=['thetaxy','thetaxz','thetayz'])


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
pressure_diffs = ['100','1000']
# pressure_diffs = ['500','1000']
for pressure_diff in pressure_diffs:
    existingMolDF = pd.read_csv('./csv_intersect/pagap'+gap_size+'_nve'+pressure_diff+'_mols.csv',index_col=[0],header=[0,1])
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
        i=0
        for waterFile in waterArr:
            if i > 500:
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

        df = df.sort_values(by=['mol','id'])
        # create per molecule dataframe to analyze at finer scale
        # atoms of a molecule are expected to be consecutive
        grouped_df = df.groupby(df['mol'])
        dfMol[['thetaxy','thetaxz','thetayz']] = grouped_df[['x','y','z']].apply(find_thetas_planes)
        dfMol['mass'] = grouped_df['mass'].sum()
        dfMol['xmass'] = grouped_df['xmass'].sum()
        dfMol['ymass'] = grouped_df['ymass'].sum()
        dfMol['zmass'] = grouped_df['zmass'].sum()

        dfMol['x'] = dfMol['xmass'] / dfMol['mass']
        dfMol['y'] = dfMol['ymass'] / dfMol['mass']
        dfMol['z'] = dfMol['zmass'] / dfMol['mass']

        dfMol['tube'] = (dfMol['z'] > inlet_z) & (dfMol['z'] < outlet_z)
        dfMolSeries.append(dfMol[['x','y','z','thetaxy','thetaxz','thetayz','tube']])
        i+=1
    print('post dfMol create')
    molDF = pd.concat(dfMolSeries,axis=1,keys=timesteps)
    combinedMolDF = pd.merge(existingMolDF, molDF, how='outer', left_index=True, right_index=True)
    combinedMolDF.to_csv('./csv_intersect/pagap'+gap_size+'_nve'+pressure_diff+'_mols.csv')
    
    del combinedMolDF
    del existingMolDF
    del molDF
    gc.collect()