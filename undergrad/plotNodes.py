import numpy as np
from numpy import diff
import time
import pmmoto
import warnings
from matplotlib import pyplot as plt
import glob
import csv
import re
from numpy import prod
plt.rcParams.update({'font.size': 16})
###Want to read in CSV data in a 2x2 shaped array

filename = 'final_results/grid_independence_data_final.csv'

# Define a regular expression pattern to extract numeric values including the identifier
pattern = re.compile(r'(\d+)\s+\(([\d.,\s]+)\)')

vol = 4.65454e6
nodes = []
mink0 = []
mink1 = []
mink2 = []
mink3 = []

with open(filename, 'r') as file:
    for line in file:
        # Split the line into identifier and tuple components
        parts = line.split()

        # Extract identifier and numeric values
        node = int(parts[0])
        MF0 = [float(x) for x in parts[1][1:-1].split(',')]
        MF1 = [float(x) for x in parts[2][0:-1].split(',')]
        MF2 = [float(x) for x in parts[3][0:-1].split(',')]
        MF3 = [float(x) for x in parts[4][0:-1].split(',')]

        # Append values to the lists
        nodes.append(node)
        mink0.append(MF0)
        mink1.append(MF1)
        mink2.append(MF2)
        mink3.append(MF3)

mink0 = np.array(mink0)
mink1 = np.array(mink1)
mink2 = np.array(mink2)
mink3 = np.array(mink3)
#Normalize minks by volume so divide each MF list by volume to give a porosity
mink0norm = mink0/vol
mink1norm = mink1/vol
mink2norm = mink2/vol
mink3norm = mink3/vol
###Now we have a nodes list from n = n to n = nMax with steps of 100 and lists of each MF corresponding to node size
        
flag_plotMF = False
if flag_plotMF:
    ###Plot nodes vs each MF, we need every MF to reach a rate of 0 in change per node count
    plt.style.use('seaborn-v0_8-whitegrid')
    plt.figure(figsize=(10,6))
    ###Subplot 1 TOP LEFT
    plt.subplot(2,2,1)
    plt.plot(nodes, mink0norm, color ='mediumorchid')
    plt.grid(True, linestyle='--', alpha=0.5)
    plt.xlim(100, 3000)
    plt.ylim(0, 0.5)
    plt.xlabel('Voxels ($n^{1/3}$)')
    plt.ylabel('MF$_0$ (-)')
    #plt.title('Nodes vs. Minkowski Functional 0')

    ###Subplot 2 TOP RIGHT
    plt.subplot(2,2,2)
    plt.plot(nodes, mink1norm, color ='lightseagreen')
    plt.grid(True, linestyle='--', alpha=0.5)
    plt.xlim(100, 3000)
    plt.ylim(bottom=0)
    plt.xlabel('Voxels ($n^{1/3}$)')
    plt.ylabel('MF$_1$ $\\left(\\frac{1}{L}\\right)$')
    #plt.title('Nodes vs. Minkowski Functional 1')

    ###Subplot 3 BOTTOM LEFT
    plt.subplot(2,2,3)
    plt.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
    plt.plot(nodes, mink2norm, color ='dodgerblue')
    plt.grid(True, linestyle='--', alpha=0.5)
    plt.xlim(100, 3000)
    plt.ylim(bottom=0)
    plt.xlabel('Voxels ($n^{1/3}$)')
    plt.ylabel('MF$_2$ $\\left(\\frac{1}{L^2}\\right)$')
    #plt.title('Nodes vs. Minkowski Functional 2')

    ###Subplot 4 BOTTOM RIGHT
    plt.subplot(2,2,4)
    plt.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
    plt.plot(nodes, mink3norm, color ='mediumslateblue')
    plt.grid(True, linestyle='--', alpha=0.5)
    plt.xlim(100, 3000)
    plt.xlabel('Voxels ($n^{1/3}$)')
    plt.ylabel('MF$_3$ $\\left(\\frac{1}{L^3}\\right)$')
    #plt.title('Nodes vs. Minkowski Functional 3')

    ###Save figure and prevent clipping
    plt.tight_layout()
    plt.savefig('plots/nodes_vs_mink.png')

###We want to plot relative change (derivative) of the vs. the node amount
#initialize new lists for change values
    
#divide diff by old MF value trailing
dnodes = diff(nodes,axis=0)
dm0dn = diff(mink0,axis=0)/dnodes
dm1dn = diff(mink1,axis=0)/dnodes
dm2dn = diff(mink2,axis=0)/dnodes
dm3dn = diff(mink3,axis=0)/dnodes
nodes_approx = (np.array(nodes)[:-1]+np.array(nodes)[1:])/2

list_lengths = len(mink0)-1

#relative change in m0
i = 0
dm0 = []
dm02 = []
dif = diff(mink0,axis=0)
for i in range(list_lengths):
    value = dif[i]/mink0[i]
    value2 = abs(dif[i]/(nodes[i]-nodes[i+1]))/mink0[i]
    dm0.append(value)
    dm02.append(value2)
dmink0 = [abs(prod(x)*100) for x in dm0]

#relative change in m1
i = 0
dm1 = []
dm12 = []
dif = diff(mink1,axis=0)
for i in range(list_lengths):
    value = dif[i]/mink1[i]
    value2 = abs(dif[i]/(nodes[i]-nodes[i+1]))/mink1[i]
    dm1.append(value)
    dm12.append(value2)
dmink1 = [abs(prod(x)*100) for x in dm1]

#relative change in m2
i = 0
dm2 = []
dm22 = []
dif = diff(mink2,axis=0)
for i in range(list_lengths):
    value = dif[i]/mink2[i]
    value2 = abs(dif[i]/(nodes[i]-nodes[i+1]))/mink2[i]
    dm2.append(value)
    dm22.append(value2)
dmink2 = [abs(prod(x)*100) for x in dm2]

#relative change in m3
i = 0
dm3 = []
dm32 = []
dif = diff(mink3,axis=0)
for i in range(list_lengths):
    value = dif[i]/mink3[i]
    value2 = abs(dif[i]/(nodes[i]-nodes[i+1])/mink3[i])
    dm3.append(value)
    dm32.append(value2)
dmink3 = [abs(prod(x)*100) for x in dm3]

flag_plot_relative = True
if flag_plot_relative:
    plt.style.use('seaborn-v0_8-whitegrid')
    plt.figure(figsize=(10,6)) #size is in inches
    #Relative Change Subplot 1 Top Left
    plt.subplot(2,2,1)
    #plt.title('Relative Change in 0th Minkowski Functional')
    plt.semilogy(nodes_approx, dmink0, color='mediumorchid')
    plt.grid(True, linestyle='--', alpha=0.5)
    plt.xlabel('Voxels ($n^{1/3}$)')
    plt.ylabel('MF$_0$ Percent Change')


    #Relative Change Subplot 2 Top Right
    plt.subplot(2,2,2)
    #plt.title('Relative Change in 1st Minkowski Functional')
    plt.semilogy(nodes_approx, dmink1, color='lightseagreen')
    plt.grid(True, linestyle='--', alpha=0.5)
    plt.xlabel('Voxels ($n^{1/3}$)')
    plt.ylabel('MF$_1$ Percent Change')

    #Relative Change Subplot 3 Bottom Left
    plt.subplot(2,2,3)
    plt.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
    #plt.title('Relative Change in 2nd Minkowski Functional')
    plt.semilogy(nodes_approx, dmink2, color='dodgerblue')
    plt.grid(True, linestyle='--', alpha=0.5)
    plt.xlabel('Voxels ($n^{1/3}$)')
    plt.ylabel('MF$_2$ Percent Change')

    #Relative Change Subplot 4 Bottom Right
    plt.subplot(2,2,4)
    plt.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
    #plt.title('Relative Change in 3rd Minkowski Functional')
    plt.semilogy(nodes_approx, dmink3, color='mediumslateblue')
    plt.grid(True, linestyle='--', alpha=0.5)
    plt.xlabel('Voxels ($n^{1/3}$)')
    plt.ylabel('MF$_3$ Percent Change')

    plt.tight_layout()
    plt.savefig('plots/plotMF_relative.png')

flag_plot_relative2 = True
if flag_plot_relative2:
    plt.style.use('seaborn-v0_8-whitegrid')
    plt.figure(figsize=(10,6)) #size is in inches
    #Relative Change Subplot 1 Top Left
    plt.subplot(2,2,1)
    #plt.title('Relative Change in 0th Minkowski Functional')
    plt.semilogy(nodes_approx, dm02, color='mediumorchid')
    plt.grid(True, linestyle='--', alpha=0.5)
    plt.xlabel('Voxels ($n^{1/3}$)')
    plt.ylabel('MF$_0$ Percent Change')


    #Relative Change Subplot 2 Top Right
    plt.subplot(2,2,2)
    #plt.title('Relative Change in 1st Minkowski Functional')
    plt.semilogy(nodes_approx, dm12, color='lightseagreen')
    plt.grid(True, linestyle='--', alpha=0.5)
    plt.xlabel('Voxels ($n^{1/3}$)')
    plt.ylabel('MF$_1$ Percent Change')

    #Relative Change Subplot 3 Bottom Left
    plt.subplot(2,2,3)
    plt.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
    #plt.title('Relative Change in 2nd Minkowski Functional')
    plt.semilogy(nodes_approx, dm22, color='dodgerblue')
    plt.grid(True, linestyle='--', alpha=0.5)
    plt.xlabel('Voxels ($n^{1/3}$)')
    plt.ylabel('MF$_2$ Percent Change')

    #Relative Change Subplot 4 Bottom Right
    plt.subplot(2,2,4)
    plt.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
    #plt.title('Relative Change in 3rd Minkowski Functional')
    plt.semilogy(nodes_approx, dm32, color='mediumslateblue')
    plt.grid(True, linestyle='--', alpha=0.5)
    plt.xlabel('Voxels ($n^{1/3}$)')
    plt.ylabel('MF$_3$ Percent Change')

    plt.tight_layout()
    plt.savefig('plots/plotMF_relative2.png')