import numpy as np
import matplotlib.pyplot as plt
import re 
import glob
import csv 

plt.rcParams.update({'font.size': 16})

file_in = "final_results/membrane_results.csv"

#Define a regular expression pattern to extract numeric values including the identifier
pattern = re.compile(r'(\d+)\s+\(([\d.,\s]+)\)')

vol = 4.65454e6
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

mink0 = np.array(mink0)
mink1 = np.array(mink1)
mink2 = np.array(mink2)
mink3 = np.array(mink3)

mink0 = mink0/vol
mink1 = mink1/vol
mink2 = mink2/vol
mink3 = mink3/vol

minkowskis = 4
mink_list = [mink0,mink1,mink2,mink3]

###Create sim_time list for plot
sim_time = []
i = 0
value = 0.01
for i in range(0,200):
    sim_time.append(value)
    value = value+0.01
    i = i+1

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

###Now plot 4 graphs of sim_time vs. MF0 - MF3

flag_plot_MFvsSim = True
if flag_plot_MFvsSim:
    colors = ['mediumorchid','lightseagreen','dodgerblue','mediumslateblue']
    plot_num = 1
    plt.style.use('seaborn-v0_8-whitegrid')
    plt.figure(figsize=(12,8))
    plt.grid(True,linestyle='--',alpha=0.75)
    for x in range(0,minkowskis):
        plt.subplot(2,2,plot_num)
        plt.hlines(y=means[x],xmin=0,xmax=2.0,linestyle='--',color='black',label='mean value')
        plt.plot(sim_time,mink_list[x],'-',color=colors[x])
        plt.xlabel('Simulation Time (nanoseconds)')
        plt.ylabel(f'MF$_{x}$')
        if x < 2:
            plt.legend(prop={'size':14},facecolor='white',framealpha=1,loc='lower right')
        else:
            plt.legend(prop={'size':14},facecolor='white',framealpha=1,loc='upper right')
        plt.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
        plt.grid(True,linestyle='--',alpha=0.75)
        x = x+1
        plot_num = plot_num + 1
    plt.tight_layout()
    #plt.show()
    plt.savefig('plots/MF_vs_simtime')


flag_histogram_MF_mean = True
if flag_histogram_MF_mean:
    ###We can now plot 4 separate histograms to find variance of each MF to their mean
    colors = ["mediumorchid",'lightseagreen','dodgerblue','mediumslateblue']
    plt.style.use('seaborn-v0_8-whitegrid')
    plt.figure(figsize=(12,8))
    plt.grid(True,linestyle='--',alpha=1)
    ###Top Left Subplot 1
    plt.subplot(2,2,1)
    plt.hist(mink0,bins=15, color=colors[0], edgecolor='black')
    #plt.xlim(1.75e6,1.80e6)
    plt.xlabel('MF$_0$')
    plt.ylabel('Frequency')
    plt.ticklabel_format(axis='x', style='sci', scilimits=(0,0))
    
    ###Top Right Subplot 2
    plt.subplot(2,2,2)
    plt.hist(mink1,bins=15,color=colors[1],edgecolor='black')
    #plt.xlim(2.32e6, 2.35e6)
    plt.xlabel('MF$_1$')
    plt.ylabel('Frequency')
    plt.ticklabel_format(axis='x', style='sci', scilimits=(0,0))

    ###Bot Left Subplot 3
    plt.subplot(2,2,3)
    plt.hist(mink2,bins=15,color=colors[2],edgecolor='black')
    #plt.xlim(-6.0e4,6.0e4)
    plt.xlabel('MF$_2$')
    plt.ylabel('Frequency')
    plt.ticklabel_format(axis='x', style='sci', scilimits=(0,0))

    ###Bot Right Subplot 4
    plt.subplot(2,2,4)
    plt.hist(mink3,bins=15,color=colors[3],edgecolor='black')
    #plt.xlim(-6e3,6e3)
    plt.xlabel('MF$_3$')
    plt.ylabel('Frequency')
    plt.ticklabel_format(axis='x', style='sci', scilimits=(0,0))

    
    plt.tight_layout()
    #plt.show()
    plt.savefig('plots/MF_histograms.png')

###We have means so we can calculate the standard deviation of each MF data
mink0_sigma = np.std(mink0)
mink1_sigma = np.std(mink1)
mink2_sigma = np.std(mink2)
mink3_sigma = np.std(mink3)
    