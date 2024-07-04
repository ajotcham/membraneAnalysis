import numpy as np
import matplotlib.pyplot as plt
from mpi4py import MPI
import os 
import glob
import natsort
#Read in locations and velocities, find atom with highest instanteous velocity, track that atom's position through simulation

###Group all relevant output files from lammp simulation
files = natsort.natsorted(glob.glob('output/dump_md.*.lammpstrj'))
numAtoms = 1149
numSteps = 301
numTotal = numAtoms * numSteps
###Initialize numpy arrays
atomPositions = np.zeros((numTotal, 5)) ### atom number, atom type, x, y, z positions for all atoms at each timestep
atomVelocities = np.zeros((numTotal, 5)) ### atom number, atom type, x, y, z velocities for all atoms at each timestep
c = 0

for filename in files:
    with open(filename, 'r') as file:
        Lines = file.readlines()[9:] #Skips first 9 rows for headers for every file
        for line in Lines:
            #Atom Positions

            atomPositions[c, 0] = float((line.split(" "))[0])
            atomPositions[c, 1] = float((line.split(" "))[1])
            atomPositions[c, 2] = float((line.split(" "))[2])
            atomPositions[c, 3] = float((line.split(" "))[3])
            atomPositions[c, 4] = float((line.split(" "))[4])
            #Atom Velocities

            atomVelocities[c, 0] = float((line.split(" "))[0])
            atomVelocities[c, 1] = float((line.split(" "))[1])
            atomVelocities[c, 2] = float((line.split(" "))[5])
            atomVelocities[c, 3] = float((line.split(" "))[6])
            atomVelocities[c, 4] = float((line.split(" "))[7])
            c = c+1
    
veloMagnitude = np.zeros((numTotal, 3))
i = 0
for i in range(numTotal):
    ###Need to find sum of velocity components to find magnitude
    veloMagnitude[i, 0] = atomVelocities[i, 0]
    veloMagnitude[i, 1] = atomVelocities[i, 1]
    veloMagnitude[i, 2] = abs(sum([atomVelocities[i, 2], atomVelocities[i, 3], atomVelocities[i, 4]])) ###Only care about greatest magnitude = moving fastest
    i = i + 1

###Find max and min value in veloMagnitude array

maxVelo = np.max(veloMagnitude[:,2])

###Indexes where max and mins are found, can use these to track trajectories since it will yield atom# and id type
maxIndex = np.asfarray(np.where(veloMagnitude[:,2] == maxVelo))
maxIndex = int(maxIndex.item())

###Returns atom number
max_atom = (veloMagnitude[maxIndex, 0])

###Track trajectory
maxTrajectory = np.zeros((numSteps, 4))
max_atomIndexes = np.where(atomPositions[:,0] == max_atom)[0].tolist()
j = 0
for value in max_atomIndexes:
    print(value)
    maxTrajectory[j, 0] = j ###Time step
    maxTrajectory[j, 1] = atomPositions[value, 2] # x position
    maxTrajectory[j, 2] = atomPositions[value, 3] # y position
    maxTrajectory[j, 3] = atomPositions[value, 4] # z position
    j = j + 1

###Returns x, y, z positions of atom that had highest instantaneous velocity magnitude
print(maxTrajectory)
