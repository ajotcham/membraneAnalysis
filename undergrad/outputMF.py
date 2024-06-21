import numpy as np
from mpi4py import MPI
import time
import pmmoto
import warnings
import matplotlib.pyplot as plt
import glob
import csv
import natsort
from datetime import datetime
from itertools import zip_longest

warnings.filterwarnings("ignore", category=DeprecationWarning) 

import cProfile


# def profile(filename=None, comm=MPI.COMM_WORLD):
#   def prof_decorator(f):
#     def wrap_f(*args, **kwargs):
#       pr = cProfile.Profile()
#       pr.enable()
#       result = f(*args, **kwargs)
#       pr.disable()

#       if filename is None:
#         pr.print_stats()
#       else:
#         filename_r = filename + ".{}".format(comm.rank)
#         pr.dump_stats(filename_r)

#       return result
#     return wrap_f
#   return prof_decorator

# @profile(filename="profile_out")

def my_function():

    comm = MPI.COMM_WORLD
    size = comm.Get_size()
    rank = comm.Get_rank()

    subdomains = [2,2,2] # Specifies how Domain is broken among procs
    ###increase number of nodes vs. MF values 
    nodes = [20,20,20] # Total Number of Nodes in Domain, controls resolution, static for this script
  
    ## Ordering for Inlet/Outlet ( (-x,+x) , (-y,+y) , (-z,+z) )
    boundaries = [[2,2],[2,2],[0,0]] # 0: Nothing Assumed  1: Walls 2: Periodic
    inlet  = [[0,0],[0,0],[0,0]]
    outlet = [[0,0],[0,0],[0,0]]
    
    if rank == 0:
      start_time = time.time()

    num_files = 200 ###Limits amount of files in test domains

    r_lookup_file = './rLookups/PA.rLookup'
    r_lookup = pmmoto.io.read_r_lookup_file(r_lookup_file)

    lammps_files = natsort.natsorted(glob.glob('testDomains/*.gz')) # list of files used to iterate
    lammps_files = lammps_files[0:num_files]

    out_filename = f"output_{datetime.now().strftime('%Y%m%d_%H%M%S')}" + ".csv"
    out_file = open(out_filename,'w',encoding='utf-8')

    for filename in lammps_files:
        sphere_data,domain_data = pmmoto.io.read_lammps_atoms(filename,r_lookup)
        
        # uses the default boundaries from the input file for x and y
        # trims in the z direction to only grab a central slice of membrane
        domain_data[2] = [-75.329262, 75.329262]

        sd = pmmoto.initialize(rank,size,subdomains,nodes,boundaries,inlet,outlet)
        pm = pmmoto.domain_generation.gen_pm_verlet_spheres(sd,sphere_data,domain_data,verlet=[12,12,12],add_periodic=True) 
        #verlet spheres will stay static for this script

        mink = pmmoto.analysis.minkowski.functionals(sd,pm.grid)


        # if rank == 0:
        #     end_time = time.time()
        #     elapsed_time = end_time - start_time
        #     print("elapsed time: ", elapsed_time)

        save_grid = False
        if save_grid:
            pmmoto.io.save_grid_data("dataOut/test_lammps_read_grid_trim",sd,pm.grid)
        
### Now create csv output file with testdomain and corresponding MF
        if rank == 0:
            out_file.write(f'{filename}\t{mink}\n') 

    out_file.close()

if __name__ == "__main__":
    my_function()
    MPI.Finalize()
