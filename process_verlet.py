import numpy as np
import matplotlib.pyplot as plt
file_in = "final_results/verlet_dependence_4.csv"

data = np.loadtxt(file_in)
num_entries = len(data)
nodes_per_verlet = np.zeros(num_entries)
proc = 4**3
###Want to calculate nodes per verlet
for n_row,row in enumerate(data): ###Loops through rows
    total_nodes = row[0]**3
    total_nodes_proc = total_nodes/proc
    total_nodes_verlet = total_nodes_proc/row[1]**3
    nodes_per_verlet[n_row] = total_nodes_verlet**(1/3)
    #print(n_row,row)

#List of each Node size assuming 24 verlet domain tested
time1000 = data[0:23,2]
time1250 = data[24:47,2]
time1500 = data[48:71,2]
time1750 = data[72:95,2]
time2000 = data[96:119,2]
time2250 = data[120:143,2]

#Initialize Plot
plt.title(f'Nodes in Each Verlet Domain vs. Time Elapsed ({proc} Processors)')
plt.xlabel('Total Nodes in Each Verlet Domain (n)')
plt.ylabel('Elapsed Time (seconds)')
plt.grid(True,linestyle='--',alpha=0.75)

#loop plot
verlet_size = 24 #Change if verlets tested per node step is different
node_step = 250 #Change if node step is different (i.e. by 100 step)
i = 0 #doesn't need to change
n = 1 #doesn't need to change
nodes_sizes_tested = 6
colors = ['dodgerblue','lightseagreen','mediumorchid','coral','gold','forestgreen'] #Add more colors to match nodes_sizes length or else it will be randomized
for x in range(0,nodes_sizes_tested):
    plt.plot(nodes_per_verlet[i:(verlet_size*n)-1],data[i:(verlet_size*n)-1,2],'o',label=f'${1000+(x*node_step)}^3$ nodes', color=colors[x])
    print(i)
    i = i+verlet_size
    n = n+1
    x = x+1
plt.legend()
plt.xlim(0,100)
plt.ylim(0,100)
plt.show()
###Find minimum time for each node size


###Not Looping Plot
# plt.plot(nodes_per_verlet[0:23],data[0:23,2],'o',label='$1000^3$ nodes',color='dodgerblue') #n=1000 
# plt.plot(nodes_per_verlet[24:47],data[24:47,2],'o',label='$1250^3$ nodes',color='lightseagreen') #n=1250
# plt.plot(nodes_per_verlet[48:71],data[48:71,2],'o',label='$1500^3$ nodes',color='mediumorchid') #n=1500
# plt.plot(nodes_per_verlet[72:95],data[72:95,2],'o',label='$1750^3$ nodes',color='coral') #n=1750
# plt.plot(nodes_per_verlet[96:119],data[96:119,2],'o',label='$2000^3$ nodes',color='gold') #n=2000
# plt.plot(nodes_per_verlet[120:143],data[120:143,2],'o',label='$2250^3$ nodes',color='forestgreen')#n=2250
# plt.xlim(0,100)
# plt.ylim(0,100)
# plt.legend()
# plt.show()