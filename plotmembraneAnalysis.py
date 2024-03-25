import numpy as np
import matplotlib.pyplot as plt
import re 
import glob
import csv 

plt.rcParams.update({'font.size': 16})

file_in = "final_results/membrane_results.csv"

#Define a regular expression pattern to extract numeric values including the identifier
pattern = re.compile(r'(\d+)\s+\(([\d.,\s]+)\)')

nodes = []
mink0 = []
mink1 = []
mink2 = []
mink3 = []

with open(file_in, 'r') as file:
    for line in file:
        # Split the line into identifier and tuple components
        parts = line.split()

        # Extract identifier and numeric values
        #node = int(parts[0])
        MF0 = [float(x) for x in parts[1][1:-1].split(',')]
        MF1 = [float(x) for x in parts[2][0:-1].split(',')]
        MF2 = [float(x) for x in parts[3][0:-1].split(',')]
        MF3 = [float(x) for x in parts[4][0:-1].split(',')]

        # Append values to the lists
        #nodes.append(node)
        mink0.append(MF0)
        mink1.append(MF1)
        mink2.append(MF2)
        mink3.append(MF3)

###Create sim_time list for plot
sim_time = []
i = 0
value = 0.01
for i in range(0,200):
    sim_time.append(value)
    value = value+0.01
    i = i+1

minkowskis = 4
mink_list = [mink0,mink1,mink2,mink3]

###Now plot 4 graphs of sim_time vs. MF0 - MF3

flag_plot_MFvsSim = False
if flag_plot_MFvsSim:
    colors = ['mediumorchid','lightseagreen','dodgerblue','mediumslateblue']
    plot_num = 1
    plt.style.use('seaborn-v0_8-whitegrid')
    plt.figure(figsize=(10,6))
    plt.grid(True,linestyle='--',alpha=0.75)
    for x in range(0,minkowskis):
        plt.subplot(2,2,plot_num)
        plt.plot(sim_time,mink_list[x],'-',color=colors[x])
        plt.xlabel('Simulation Time (nanoseconds)')
        plt.ylabel(f'MF$_{x}$')
        plt.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
        plt.grid(True,linestyle='--',alpha=0.75)
        if x == 3:
            plt.ylim(0, 1e5)
        elif x == 2:
            plt.ylim(0,1e6)
        else:
            plt.ylim(0, 2.5e6)
        x = x+1
        plot_num = plot_num + 1
    #plt.show()
    plt.tight_layout()
    plt.savefig('plots/MF_vs_simtime')


means = []
 ###Need to have new lists of MF variation from mean
for x in range(0,minkowskis):
        value = np.mean(mink_list[x])
        means.append(value)

###Now subtract each mink list values from their respective means
mink0_diff = mink0 - means[0]
mink1_diff = mink1 - means[1]
mink2_diff = mink2 - means[2]
mink3_diff = mink3 - means[3]


flag_histogram_MF_mean = False
if flag_histogram_MF_mean:
     ###We can now plot 4 separate histograms to find variance of each MF to their mean
    colors = ['mediumorchid','lightseagreen','dodgerblue','mediumslateblue']
    plt.style.use('seaborn-v0_8-whitegrid')
    plt.figure(figsize=(12,8))
    plt.grid(True,linestyle='--',alpha=1)
    ###Top Left Subplot 1
    plt.subplot(2,2,1)
    plt.hist(mink0_diff,bins=15,color=colors[0],edgecolor='black')
    plt.xlim(-2.0e4,2.0e4)
    plt.xlabel('MF$_0$ Difference from Mean')
    plt.ylabel('Frequency')
    plt.ticklabel_format(axis='x', style='sci', scilimits=(0,0))
    
    ###Top Right Subplot 2
    plt.subplot(2,2,2)
    plt.hist(mink1_diff,bins=15,color=colors[1],edgecolor='black')
    plt.xlim(-1.5e4,1.5e4)
    plt.xlabel('MF$_1$ Difference from Mean')
    plt.ylabel('Frequency')
    plt.ticklabel_format(axis='x', style='sci', scilimits=(0,0))

    ###Bot Left Subplot 3
    plt.subplot(2,2,3)
    plt.hist(mink2_diff,bins=15,color=colors[2],edgecolor='black')
    plt.xlim(-6.0e4,6.0e4)
    plt.xlabel('MF$_2$ Difference from Mean')
    plt.ylabel('Frequency')
    plt.ticklabel_format(axis='x', style='sci', scilimits=(0,0))

    ###Bot Right Subplot 4
    plt.subplot(2,2,4)
    plt.hist(mink3_diff,bins=20, color=colors[3],edgecolor='black')
    plt.xlim(-6e3,6e3)
    plt.xlabel('MF$_3$ Difference from Mean')
    plt.ylabel('Frequency')
    plt.ticklabel_format(axis='x', style='sci', scilimits=(0,0))

    #plt.show()
    plt.tight_layout()
    plt.savefig('plots/MF_histograms.png')

###We have means so we can calculate the standard deviation of each MF data
mink0_sigma = np.std(mink0)
mink1_sigma = np.std(mink1)
mink2_sigma = np.std(mink2)
mink3_sigma = np.std(mink3)
    