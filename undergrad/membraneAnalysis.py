import numpy as np
from mpi4py import MPI
import time
import pmmoto
import warnings
import matplotlib.pyplot as plt
import glob
import natsort
from datetime import datetime
warnings.filterwarnings("ignore", category=DeprecationWarning) 

import cProfile


def my_function():

    comm = MPI.COMM_WORLD
    size = comm.Get_size()
    rank = comm.Get_rank()

    subdomains = [2,2,2] # Specifies how Domain is broken among procs
    ###increase number of nodes vs. MF values 

    n = 100
    
    boundaries = [[2,2],[2,2],[0,0]] # 0: Nothing Assumed  1: Walls 2: Periodic
    inlet  = [[0,0],[0,0],[0,0]]
    outlet = [[0,0],[0,0],[0,0]]

    if rank == 0:
        out_filename = f"grid_independence_{datetime.now().strftime('%Y%m%d_%H%M%S')}" + ".csv"
        out_file = open(out_filename,'w',encoding='utf-8')

    nodes = [n,n,n] # Total Number of Nodes in Domain, controls resolution
    # xArray.append(n)
    
    if rank == 0:
        start_time = time.time()

    
    r_lookup_file = './rLookups/PA.rLookup'
    r_lookup = pmmoto.io.read_r_lookup_file(r_lookup_file)
    
    num_files = 200 ###Limits amount of files in test domains

    lammps_files = natsort.natsorted(glob.glob('testDomains/*.gz')) # list of files used to iterate
    lammps_files = lammps_files[0:num_files]

    for lammps_file in lammps_files:
        
        sphere_data,domain_data = pmmoto.io.read_lammps_atoms(lammps_file,r_lookup)
    
        # uses the default boundaries from the input file for x and y
        # trims in the z direction to only grab a central slice of membrane
        domain_data[2] = [-75.329262, 75.329262]

        #Verlet is static for iterating through nodes
        sd = pmmoto.initialize(rank,size,subdomains,nodes,boundaries,inlet,outlet)
        pm = pmmoto.domain_generation.gen_pm_verlet_spheres(sd,sphere_data,domain_data,verlet=[7,7,7],add_periodic=True) 
        #Cpu time vs. # of verlet domains and change nodes if time available, iterate through verlet and node sizes

        mink = pmmoto.analysis.minkowski.functionals(sd,pm.grid)

        if rank == 0:
            out_file.write(f'{lammps_file}\t{mink}\n') 
            end_time = time.time()
            elapsed_time = end_time - start_time
            out_file.flush()
            print("elapsed time: ", elapsed_time,flush=True)

        save_grid = False
        if save_grid:
            pmmoto.io.save_grid_data("dataOut/test_lammps_read_grid_trim",sd,pm.grid)

if __name__ == "__main__":
    my_function()
    MPI.Finalize()
