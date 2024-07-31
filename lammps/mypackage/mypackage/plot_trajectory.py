"""plotting and math packages"""
import matplotlib.pyplot as plt
import numpy as np
import natsort
import glob
from mypackage import readlammps

files = files = natsort.natsorted(glob.glob('output/dump_md.*.lammpstrj'))
# boundaries = np.zeros((3, 2))
# atom_info = process_atoms(files)
###To plot trajectories first need to calculate magnitude of velocity of all atoms
###Equation for magnitude of velocity = sqrt(a^2 + b^2 + c^2)

###Logic -> since process_atoms updates the np array every new file then points need to be plotted as it updates
###So a loop for process_atoms should be called that loops through files then calculate velocity magnitudes and plot trajectory

def get_velocity_mag(vx,vy,vz):
    """
    Calculates magnitude given velocity components
    """
    return np.sqrt(vx**2 + vy**2 + vz**2)



def calc_velo(atom_velo):
    """
    calculates magnitudes of atom velocities
    """

    mag_velocity = get_velocity_mag(atom_velo[:,0], atom_velo[:,1], atom_velo[:,2])

    # #Now calculate magnitude of velocity for all lines in atom_info
    # for c,lines in enumerate(atom_info):
    #     v_mag = ((atom_info[c,5])**2 + (atom_info[c,6])**2 + (atom_info[c,7])**2)**0.5
    #     #inserts velocity magnitude into atom_info numpy array
    #     atom_info[c,8] = v_mag
    #     if v_mag > max_velo:
    #         max_velo = v_mag
    #         max_info = atom_info[c] #gives other info about highest velo atom
    #         max_timestep = file.split(".")[1] ###Gives timestep
    #         max_file = file ###Gives file the highest instant velo is found in

    return mag_velocity

def top_velo(mag_velocity, x):
    """
    finds top x% of velocity magnitude given array of magnitudes and x%
    """
    ###Finds the X percentile value for comparison
    #100 - x so the arguement is more intuitive (ex. top 10% pass in 10 instead of 90)
    threshold = np.percentile(mag_velocity, 100 - x)

    ###Finds the top x% of atoms in the passed in array and indices where threshold is true
    top_indices = np.where(mag_velocity >= threshold)[0]
    top_atoms = mag_velocity[top_indices]
    return top_atoms, top_indices

def track_trajectory(atom_pos, top_indices):
    """
    tracks trajectories of top velo atoms given positions and indices of top atoms
    """
    return

def plot_trajectory(atom_pos, boundaries):
    """
    plots trajectories using matplotlib
    """
    ###Data is in the shape of a 2 column array
    # plt.imshow(atom_pos, cmap='viridis', interpolation='none', extent= boundaries)

    # plt.colorbar()
    plt.plot(atom_pos[:,0], atom_pos[:,1], 'o')

    plt.title('Sample Plot')
    plt.xlabel('X Position')
    plt.ylabel('Y Positions')
    return plt.show()

def plot_individual_location(top_indices, atom_pos):
    """
    plot atom positions based on provided index list
    """
    for index in top_indices:
        plt.plot(atom_pos[index,0], atom_pos[index,1], 'o')

    plt.show()

def find_grid_cell(atom_positions, boundaries, discretization, lengths):
    """
    Converts atom's position to grid cell location based on closest node
    """
    dims = atom_positions.shape[0]
    cell_lengths = np.zeros((dims))
    for n in range(dims):
        cell_lengths[n] = lengths[n]/discretization[n]


    ###Calculate cell indice
    grid_indices = np.zeros((dims),dtype=np.int64)
    for n in range(dims):
        grid_indices[n] = int((atom_positions[n] - boundaries[n,0]) / cell_lengths[n])
        
    # # Makes sure indices are inside boundaries
    # column = min(max(column, 0), discretization[0] - 1) #subtract 1 since indice range from 0 - 4
    # row = min(max(row, 0), discretization[1] - 1)

    return grid_indices
