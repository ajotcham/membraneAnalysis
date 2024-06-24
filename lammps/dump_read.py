import numpy as np
from numpy import diff
import time
import pmmoto
import os
import gzip
import pyvista as pv
from mpi4py import MPI
from matplotlib import pyplot as plt
import itertools as IT
from numpy import prod
comm = MPI.COMM_WORLD

#Read in lammps file and track trajectories over sim time
#Then find the atom with the highest instantaneous velocity and track its location

flag_read = True
if flag_read:
    filename = "dump_md.lammpstrj"
    lammps_input = open(filename, 'r') ###Opens dump file to be read
    skip_key = "TIMESTEP" #Signals where line skips need to be everytime
    new_input = open('new_dump.lammpstrj', 'w') ###New dump file where data lines can be compiled
    lines = lammps_input.readlines()

    for line in lines:
        if skip_key in line:
            c = 0 ###resets counter everytime key is found in one of the lines
            continue 
        elif 0 <= c < 8:
            c = c + 1
            continue
        elif c >= 8:
            new_input.write(line)

    lammps_input.close()
    new_input.close()

### Code now is all data in format of Atom ID, Atom Type, Xpos, Ypos, Zpos

### Need to take data in new_dump and put in np arrays for calculations with data
