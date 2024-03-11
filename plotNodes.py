import numpy as np
import time
import pmmoto
import warnings
from matplotlib import pyplot as plt
import glob
import csv
import re

###Want to read in CSV data in a 2x2 shaped array

filename = 'grid_independence_1500.csv'

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
        
###Plot nodes vs each MF, we need every MF to reach a rate of 0 in change per node count
        
plt.style.use('seaborn-v0_8-whitegrid')

###Subplot 1 TOP LEFT
plt.subplot(2,2,1)
plt.plot(nodes, mink0, color ='mediumorchid')
plt.grid(True, linestyle='--', alpha=0.5)
plt.xlim(0, 1500)
plt.title('Nodes vs. Minkowski Functional 0')

###Subplot 2 TOP RIGHT
plt.subplot(2,2,2)
plt.plot(nodes, mink1, color ='lightseagreen')
plt.grid(True, linestyle='--', alpha=0.5)
plt.xlim(0, 1500)
plt.title('Nodes vs. Minkowski Functional 1')

###Subplot 3 BOTTOM LEFT
plt.subplot(2,2,3)
plt.plot(nodes, mink2, color ='dodgerblue')
plt.grid(True, linestyle='--', alpha=0.5)
plt.xlim(0, 1500)
plt.title('Nodes vs. Minkowski Functional 2')

###Subplot 4 BOTTOM RIGHT
plt.subplot(2,2,4)
plt.plot(nodes, mink3, color ='mediumslateblue')
plt.grid(True, linestyle='--', alpha=0.5)
plt.xlim(0, 1500)
plt.title('Nodes vs. Minkowski Functional 3')

###Save figure and prevent clipping
plt.tight_layout()  
plt.show()
