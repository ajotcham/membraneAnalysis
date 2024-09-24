import mypackage
import numpy as np
import glob as glob
import natsort
import os

# Path to the folder containing the files
folder_path = '../../output/poiseuille'

# Define the file pattern to match specific files
file_pattern = os.path.join(folder_path, '3D_poi.*.lammpstrj')

# List all files matching the pattern
all_files = natsort.natsorted(glob.glob(file_pattern))

# Select specific files (500th to 550th)
selected_files = all_files[500:550]

all_positions = []
all_velocities_mag = []
#Loop through desired timesteps and append timestep specific data to larger dataset
for file in selected_files:
    atoms, num_atoms, boundaries = mypackage.read_atoms(file)

    all_positions.append(atoms["positions"])
    mag_velocity = mypackage.calc_velo(atoms["velocities"])
    all_velocities_mag.append(mag_velocity)

all_positions = np.concatenate(all_positions, axis=0)
all_velocities_mag = np.concatenate(all_velocities_mag, axis=0)
top_atoms, top_indices = mypackage.top_velo(all_velocities_mag, 1)

#create True
PLOT = True
if PLOT:
    xbounds = (-2.5, 12.5)
    ybounds = (-10,10)
    mypackage.plot_trajectory(all_positions, xbounds, ybounds)
    mypackage.plot_individual_location(top_indices, all_positions, xbounds, ybounds)

GRID = False
if GRID:
    discretization = [5,5,5]
    grid = np.zeros((discretization[0], discretization[1]))
    lengths = np.zeros((3))
    deltas = np.zeros((3))
    for n,axis in enumerate(boundaries):
        lengths[n] = axis[1] - axis[0]
        deltas[n] = lengths[n]/discretization[n]

    #For this case delta x and delta y will be the same
    x = np.linspace(deltas[0]/2, lengths[0]- deltas[0]/2, discretization[0])
    y = np.linspace(deltas[1]/2, lengths[1]- deltas[1]/2, discretization[1])

    for index in top_indices:
        ind = mypackage.find_grid_cell(atom_pos[index, :], boundaries, discretization, lengths)
        grid[ind[0], ind[1]] += 1. #increments grid by one 

    print(grid)

HEATMAP = False
if HEATMAP:
    xbins = 150
    ybins = 150
    xbounds = (0, 10.0)
    ybounds = (-8.0, 8.0)
    mypackage.create_heatmap(all_positions, all_velocities_mag, xbins, ybins, xbounds, ybounds)

###Subdomain of interest that only considers atoms in specified region x,y,z
x_sub = (0,5)
y_sub = (-5, 5)
z_sub = (0,5)
sub_atoms = mypackage.poiseuille_subregion_atoms(atoms["positions"], atoms["velocities"], atoms["map"], x_sub, y_sub, z_sub)
