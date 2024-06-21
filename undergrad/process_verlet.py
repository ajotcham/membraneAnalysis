import numpy as np
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 16})

file_in = "final_results/verlet_independence_3.csv"

data = np.loadtxt(file_in)
num_entries = len(data)
nodes_per_verlet = np.zeros(num_entries)
verlets = np.zeros(num_entries)
proc = 3**3 ###Must change first integer to match the cube root of the amount of processes being ran

###Want to calculate nodes per verlet
for n_row,row in enumerate(data): ###Loops through rows
    total_nodes = row[0]**3
    verlets[n_row] = row[1]
    total_nodes_proc = total_nodes/proc
    total_nodes_verlet = total_nodes_proc/row[1]**3
    nodes_per_verlet[n_row] = total_nodes_verlet**(1/3)
    #print(n_row,row)

all_data = data[0:143,2]

###Find minimum time for each node size
n = 1
i = 0
d = {}
verlet_size = 24
for x in range(0,6):
    d["minimum_node{0}".format(x+1)] = min(data[i:(verlet_size*n)-1,2])
    i = i+24
    n = n+1

#Convert dictionary into separate list
minimum_times = list(d.values())
###Find index in main data list to find node size associated
indexlist = []
for x in range (0,6):
    index = np.where(all_data == minimum_times[x])[0][0]
    indexlist.append(index)
    x = x+1

### Using list of index use that to find nodes per verlet amount
nodes_min = []
for x in indexlist:
    nodes = nodes_per_verlet[x]
    nodes_min.append(nodes)

###Create list with corresponding verlet values for minimum times
verlet_min = []
for item in indexlist:
    verlet = data[item,1]
    verlet_min.append(verlet) ###use this list to highlight minimum points on plot 2

flag_plot = True
if flag_plot:
    #Initialize Plot 1
    plt.figure(figsize=(8,6))
    #plt.title(f'Nodes in Each Verlet Domain vs. Time Elapsed ({proc} Processors)')
    plt.xlabel('Voxels in Verlet Domain ($n^3$)')
    plt.ylabel('Time Elapsed (seconds)')
    plt.grid(True,linestyle='--',alpha=0.75)

    #loop plot
    verlet_size = 24 #Change if verlets tested per node step is different
    node_step = 250 #Change if node step is different (i.e. by 100 step)
    i = 0 #doesn't need to change
    n = 1 #doesn't need to change
    nodes_sizes_tested = 6
    colors = ['dodgerblue','lightseagreen','mediumorchid','coral','gold','forestgreen'] #Add more colors to match nodes_sizes length or else it will be randomized
    for x in range(0,nodes_sizes_tested):
        plt.plot(nodes_per_verlet[i:(verlet_size*n)-1],data[i:(verlet_size*n)-1,2],'--',label=f'${1000+(x*node_step)}^3$ nodes', color=colors[x])
        i = i+verlet_size
        n = n+1
        x = x+1
    plt.xlim(0,80) ###Limited so we can see the minimum points for each line
    plt.ylim(0,60)
    plt.plot(nodes_min,minimum_times,'o',color='royalblue',alpha=1.0,label='minimum points')
    plt.legend(prop={'size':12},facecolor='white',framealpha=1)
    #plt.show()
    plt.savefig(f'plots/nodes_per_verlet_elapsedtime_{proc}')


    #Initialize Plot 2
    plt.figure(figsize=(8,6))
    #plt.title(f'Nodes in Each Verlet Domain vs. Time Elapsed ({proc} Processors)')
    plt.xlabel('Total Verlet Domains ($v^3$)')
    plt.ylabel('Time Elapsed (seconds)')
    plt.grid(True,linestyle='--',alpha=0.75)

    #loop plot
    verlet_size = 24 #Change if verlets tested per node step is different
    node_step = 250 #Change if node step is different (i.e. by 100 step)
    i = 0 #doesn't need to change
    n = 1 #doesn't need to change
    nodes_sizes_tested = 6
    colors = ['dodgerblue','lightseagreen','mediumorchid','coral','gold','forestgreen'] #Add more colors to match nodes_sizes length or else it will be randomized
    for x in range(0,nodes_sizes_tested):
        plt.semilogy(verlets[i:(verlet_size*n)-1],data[i:(verlet_size*n)-1,2],'--',label=f'${1000+(x*node_step)}^3$ nodes', color=colors[x])
        i = i+verlet_size
        n = n+1
        x = x+1
    #highlight minimum points
    plt.plot(verlet_min,minimum_times,'o',color='royalblue',label='minimum points',alpha=1.0)
    plt.legend(prop={'size':12},facecolor='white',framealpha=1)
    plt.xlim(0,25) ###Limited so we can see the minimum points for each line
    if proc == 3**3:
        plt.ylim(0,10**3)
    elif proc == 4**3:
        plt.ylim(0,10**2)
    #plt.show()
    plt.savefig(f'plots/verlet_elapsedtime_{proc}')
