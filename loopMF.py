import numpy as np
from mpi4py import MPI
import time
import pmmoto
import warnings
import matplotlib.pyplot as plt
warnings.filterwarnings("ignore", category=DeprecationWarning) 

import cProfile

def profile(filename=None, comm=MPI.COMM_WORLD):
  def prof_decorator(f):
    def wrap_f(*args, **kwargs):
      pr = cProfile.Profile()
      pr.enable()
      result = f(*args, **kwargs)
      pr.disable()

      if filename is None:
        pr.print_stats()
      else:
        filename_r = filename + ".{}".format(comm.rank)
        pr.dump_stats(filename_r)

      return result
    return wrap_f
  return prof_decorator

@profile(filename="profile_out")

def my_function():

    comm = MPI.COMM_WORLD
    size = comm.Get_size()
    rank = comm.Get_rank()

    subdomains = [1,1,1] # Specifies how Domain is broken among procs
    nodes = [200,200,200] # Total Number of Nodes in Domain, controls resolution
  
    ## Ordering for Inlet/Outlet ( (-x,+x) , (-y,+y) , (-z,+z) )
    boundaries = [[2,2],[2,2],[0,0]] # 0: Nothing Assumed  1: Walls 2: Periodic
    inlet  = [[0,0],[0,0],[0,0]]
    outlet = [[0,0],[0,0],[0,0]]
    
    i = 8
    verlet_max = 33
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
      print("elapsed time: ", end_time - start_time)

      pmmoto.io.save_grid_data("dataOut/test_lammps_read_grid_trim",sd,pm.grid)
    print(x_array)
    print(y_array)
    plt.plot(x_array, y_array)
    plt.title("CPU Time vs. Verlet Sphere Count")
    plt.xlabel("Verlet Spheres ($\sqrt[3]{x}$)")
    plt.ylabel("Time Elapsed (seconds)")
    file = "verlet_vs_time_2.png"
    plt.savefig(file)


if __name__ == "__main__":
    my_function()
    MPI.Finalize()
