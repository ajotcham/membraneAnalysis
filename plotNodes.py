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
###Want to read in CSV data in a 2x2 shaped array

filename = 'final_results/grid_independence_data_final.csv'

# Define a regular expression pattern to extract numeric values including the identifier
pattern = re.compile(r'(\d+)\s+\(([\d.,\s]+)\)')

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


###Now we have a nodes list from n = n to n = nMax with steps of 100 and lists of each MF corresponding to node size
        
flag_plotMF = True
if flag_plotMF:
    ###Plot nodes vs each MF, we need every MF to reach a rate of 0 in change per node count
    plt.style.use('seaborn-v0_8-whitegrid')
    plt.figure(figsize=(10,6))
    ###Subplot 1 TOP LEFT
    plt.subplot(2,2,1)
    plt.plot(nodes, mink0, color ='mediumorchid')
    plt.grid(True, linestyle='--', alpha=0.5)
    plt.xlim(100, 3000)
    plt.ylim(0,2.0e6)
    plt.xlabel('Total Node Amount ($n^{1/3}$)')
    plt.ylabel('MF0, Volume ($l^3$)')
    plt.title('Nodes vs. Minkowski Functional 0')

    ###Subplot 2 TOP RIGHT
    plt.subplot(2,2,2)
    plt.plot(nodes, mink1, color ='lightseagreen')
    plt.grid(True, linestyle='--', alpha=0.5)
    plt.xlim(100, 3000)
    plt.ylim(bottom=0)
    plt.xlabel('Total Node Amount ($n^{1/3}$)')
    plt.ylabel('MF1, Surface Area ($l^2$)')
    plt.title('Nodes vs. Minkowski Functional 1')

    ###Subplot 3 BOTTOM LEFT
    plt.subplot(2,2,3)
    plt.plot(nodes, mink2, color ='dodgerblue')
    plt.grid(True, linestyle='--', alpha=0.5)
    plt.xlim(100, 3000)
    plt.ylim(bottom=0)
    plt.xlabel('Total Node Amount ($n^{1/3}$)')
    plt.ylabel('MF2, Integral Mean Curvature ($l$)')
    plt.title('Nodes vs. Minkowski Functional 2')

    ###Subplot 4 BOTTOM RIGHT
    plt.subplot(2,2,4)
    plt.plot(nodes, mink3, color ='mediumslateblue')
    plt.grid(True, linestyle='--', alpha=0.5)
    plt.xlim(100, 3000)
    plt.xlabel('Total Node Amount ($n^{1/3}$)')
    plt.ylabel('MF3, Euler Characteristic')
    plt.title('Nodes vs. Minkowski Functional 3')

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
for i in range(list_lengths):
    dif = diff(mink0,axis=0)
    value = dif[i]/mink0[i]
    dm0.append(value)
    i = i+1
dmink0 = [prod(x) for x in dm0]
#relative change in m1
i = 0
dm1 = []
for i in range(list_lengths):
    dif = diff(mink1,axis=0)
    value = dif[i]/mink0[i]
    dm1.append(value)
    i = i+1
dmink1 = [prod(x) for x in dm1]
#relative change in m2
i = 0
dm2 = []
for i in range(list_lengths):
    dif = diff(mink2,axis=0)
    value = dif[i]/mink0[i]
    dm2.append(value)
    i = i+1
dmink2 = [prod(x) for x in dm2]
#relative change in m3
i = 0
dm3 = []
for i in range(list_lengths):
    dif = diff(mink3,axis=0)
    value = dif[i]/mink0[i]
    dm3.append(value)
    i = i+1
dmink3 = [prod(x) for x in dm3]


flag_plot_relative = True
if flag_plot_relative:
    plt.style.use('seaborn-v0_8-whitegrid')
    plt.figure(figsize=(10,6)) #size is in inches
    #Relative Change Subplot 1 Top Left
    plt.subplot(2,2,1)
    plt.title('Relative Change in 0th Minkowski Functional')
    plt.plot(nodes_approx, dmink0, color='mediumorchid')
    plt.grid(True, linestyle='--', alpha=0.5)
    plt.xlabel('Total Nodes ($n^{1/3}$)')
    plt.ylabel('Relative Change')


    #Relative Change Subplot 2 Top Right
    plt.subplot(2,2,2)
    plt.title('Relative Change in 1st Minkowski Functional')
    plt.plot(nodes_approx, dmink1, color='lightseagreen')
    plt.grid(True, linestyle='--', alpha=0.5)
    plt.xlabel('Total Nodes ($n^{1/3}$)')
    plt.ylabel('Relative Change')

    #Relative Change Subplot 3 Bottom Left
    plt.subplot(2,2,3)
    plt.title('Relative Change in 2nd Minkowski Functional')
    plt.plot(nodes_approx, dmink2, color='dodgerblue')
    plt.grid(True, linestyle='--', alpha=0.5)
    plt.xlabel('Total Nodes ($n^{1/3}$)')
    plt.ylabel('Relative Change')

    #Relative Change Subplot 4 Bottom Right
    plt.subplot(2,2,4)
    plt.title('Relative Change in 3rd Minkowski Functional')
    plt.plot(nodes_approx, dmink3, color='mediumslateblue')
    plt.grid(True, linestyle='--', alpha=0.5)
    plt.xlabel('Total Nodes ($n^{1/3}$)')
    plt.ylabel('Relative Change')

    plt.tight_layout()
    plt.savefig('plots/plotMF_relative.png')