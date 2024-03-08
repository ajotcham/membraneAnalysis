import numpy as np
from mpi4py import MPI
import time
import pmmoto
import warnings
import matplotlib.pyplot as plt
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

    n = 100
    nMax = 501
    # xArray = []
    
    boundaries = [[2,2],[2,2],[0,0]] # 0: Nothing Assumed  1: Walls 2: Periodic
    inlet  = [[0,0],[0,0],[0,0]]
    outlet = [[0,0],[0,0],[0,0]]

    for n in range(n,nMax, 100):
        nodes = [n,n,n] # Total Number of Nodes in Domain, controls resolution
        # xArray.append(n)
        
        if rank == 0:
            start_time = time.time()

        r_lookup_file = './rLookups/PA.rLookup'
        r_lookup = pmmoto.io.read_r_lookup_file(r_lookup_file)

        lammps_file = './testDomains/membranedata.71005000.gz'
        sphere_data,domain_data = pmmoto.io.read_lammps_atoms(lammps_file,r_lookup)
        
        # uses the default boundaries from the input file for x and y
        # trims in the z direction to only grab a central slice of membrane
        domain_data[2] = [-75.329262, 75.329262]

        #Verlet is static for iterating through nodes
        sd = pmmoto.initialize(rank,size,subdomains,nodes,boundaries,inlet,outlet)
        pm = pmmoto.domain_generation.gen_pm_verlet_spheres(sd,sphere_data,domain_data,verlet=[12,12,12],add_periodic=True) 
        #Cpu time vs. # of verlet domains and change nodes if time available, iterate through verlet and node sizes

        mink = pmmoto.analysis.minkowski.functionals(sd,pm.grid)

        if rank == 0:
            print(mink)

        # if rank == 0:
        #     end_time = time.time()
        #     elapsed_time = end_time - start_time
        #     print("elapsed time: ", elapsed_time)

        save_grid = False
        if save_grid:
            pmmoto.io.save_grid_data("dataOut/test_lammps_read_grid_trim",sd,pm.grid)

    # plt.figure()
    # plt.plot(xArray,mink0,color='red')
    # plt.plot(xArray,mink1,color='blue')
    # # plt.plot(xArray,mink2,color='green')
    # # plt.plot(xArray,mink3,color='yellow')
    # plt.title("Node Sizes vs. Minkowski Functional Values")
    # plt.xlabel("Node Sizes")
    # plt.ylabel("Minkowski Functional Values")
    # file = "nodes_vs_MF.png"
    # plt.savefig(file)
if __name__ == "__main__":
    my_function()
    MPI.Finalize()
