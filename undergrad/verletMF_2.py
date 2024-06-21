import numpy as np
from mpi4py import MPI
import time
import pmmoto
import warnings
import matplotlib.pyplot as plt
from datetime import datetime
warnings.filterwarnings("ignore", category=DeprecationWarning) 

import cProfile

def my_function():

    comm = MPI.COMM_WORLD
    size = comm.Get_size()
    rank = comm.Get_rank()

    subdomains = [2,2,2] # Specifies how Domain is broken among procs
    ## Ordering for Inlet/Outlet ( (-x,+x) , (-y,+y) , (-z,+z) )
    boundaries = [[2,2],[2,2],[0,0]] # 0: Nothing Assumed  1: Walls 2: Periodic
    inlet  = [[0,0],[0,0],[0,0]]
    outlet = [[0,0],[0,0],[0,0]]

    out_filename = f"verlet_{datetime.now().strftime('%Y%m%d_%H%M%S')}" + ".csv"
    out_file = open(out_filename,'w',encoding='utf-8')
    
    for n in range(1000,2500,250):
        nodes = [n,n,n] # Total Number of Nodes in Domain, controls resolution
    
        i = 1
        verlet_max = 25
        x_array = []
        y_array = []

        ###Iterate through verlet sphere sizes (i) and compare run-times (parabolic shape)
        for i in range(i,verlet_max):
            start_time = time.time()

            r_lookup_file = './rLookups/PA.rLookup'
            r_lookup = pmmoto.io.read_r_lookup_file(r_lookup_file)

            lammps_file = './testDomains/membranedata.71005000.gz'
            sphere_data,domain_data = pmmoto.io.read_lammps_atoms(lammps_file,r_lookup)
            
            # uses the default boundaries from the input file for x and y
            # trims in the z direction to only grab a central slice of membrane
            domain_data[2] = [-75.329262, 75.329262]

            sd = pmmoto.initialize(rank,size,subdomains,nodes,boundaries,inlet,outlet)
            pm = pmmoto.domain_generation.gen_pm_verlet_spheres(sd,sphere_data,domain_data,verlet=[i,i,i],add_periodic=True) 
            #Cpu time vs. # of verlet domains and change nodes if time available, iterate through verlet and node sizes

            end_time = time.time()
            elapsed_time = end_time - start_time
            x_array.append(i)
            y_array.append(elapsed_time)


            if rank == 0:
                out_file.write(f'{n}\t{i}\t{elapsed_time}\n') 


if __name__ == "__main__":
    my_function()
    MPI.Finalize()
