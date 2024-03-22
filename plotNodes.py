import numpy as np
from numpy import diff
import time
import pmmoto
import warnings
from matplotlib import pyplot as plt
import glob
import csv
import re

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
        
flag_plotMF = False
if flag_plotMF:
    ###Plot nodes vs each MF, we need every MF to reach a rate of 0 in change per node count
    plt.style.use('seaborn-v0_8-whitegrid')

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
    plt.show()

###We want to plot relative change (derivative) of the vs. the node amount
#initialize new lists for change values
        
dnodes = diff(nodes,axis=0)
dm0dn = diff(mink0,axis=0)/dnodes
dm1dn = diff(mink1,axis=0)/dnodes
dm2dn = diff(mink2,axis=0)/dnodes
dm3dn = diff(mink3,axis=0)/dnodes
nodes_approx = (np.array(nodes)[:-1]+np.array(nodes)[1:])/2

flag_plot_relative = True
if flag_plot_relative:
    plt.style.use('seaborn-v0_8-whitegrid')

    #Relative Change Subplot 1 Top Left
    plt.subplot(2,2,1)
    plt.title('Relative Change in Minkowski Functional 0')
    plt.plot(nodes_approx, dm0dn, color='mediumorchid')
    plt.grid(True, linestyle='--', alpha=0.5)
    plt.xlabel('Total Node Amount ($n^{1/3}$)')
    plt.ylabel(r'Change in MF0 ($\frac{l^3}{n})$')


    #Relative Change Subplot 2 Top Right
    plt.subplot(2,2,2)
    plt.title('Relative Change in Minkowski Functional 1')
    plt.plot(nodes_approx, dm1dn, color='lightseagreen')
    plt.grid(True, linestyle='--', alpha=0.5)
    plt.xlabel('Total Node Amount ($n^{1/3}$)')
    plt.ylabel(r'Change in MF1 ($\frac{l^2}{n})$')

    #Relative Change Subplot 3 Bottom Left
    plt.subplot(2,2,3)
    plt.title('Relative Change in Minkowski Functional 2')
    plt.plot(nodes_approx, dm2dn, color='dodgerblue')
    plt.grid(True, linestyle='--', alpha=0.5)
    plt.xlabel('Total Node Amount ($n^{1/3}$)')
    plt.ylabel(r'Change in MF2 ($\frac{l}{n})$')

    #Relative Change Subplot 4 Bottom Right
    plt.subplot(2,2,4)
    plt.title('Relative Change in Minkowski Functional 3')
    plt.plot(nodes_approx, dm3dn, color='mediumslateblue')
    plt.grid(True, linestyle='--', alpha=0.5)
    plt.xlabel('Total Node Amount ($n^{1/3}$)')
    plt.ylabel(r'Change in MF3 ($\frac{1}{n})$')

    plt.tight_layout()
    plt.show()