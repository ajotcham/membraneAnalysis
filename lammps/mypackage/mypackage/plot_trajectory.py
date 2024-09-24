"""plotting and math packages"""
import glob
import matplotlib.pyplot as plt
import numpy as np
import natsort
from scipy.ndimage import gaussian_filter


files = files = natsort.natsorted(glob.glob('output/dump_md.*.lammpstrj'))

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


def plot_trajectory(atom_pos, xboundaries, yboundaries):
    """
    plots trajectories using matplotlib
    """
    plt.plot(atom_pos[:,0], atom_pos[:,1], 'o')

    plt.title('Sample Plot')
    plt.xlabel('X Position')
    plt.ylabel('Y Positions')
    plt.xlim(xboundaries)
    plt.ylim(yboundaries)
    return plt.show()


def plot_individual_location(top_indices, atom_pos, xboundaries, yboundaries):
    """
    plot atom positions based on provided index list
    """
    for index in top_indices:
        plt.plot(atom_pos[index,0], atom_pos[index,1], 'o')
    plt.title('Atoms With Greatest Velocity in System')
    plt.xlabel('X Position')
    plt.ylabel('Y Position')
    plt.xlim(xboundaries)
    plt.ylim(yboundaries)
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

    return grid_indices

def create_heatmap(positions, velocities, xbins, ybins, xbounds, ybounds):
    """
    creates heatmap to show hotspots of atom velocity
    """
    x = positions[:,0]
    y = positions[:,1]
    heatmap, xedges, yedges = np.histogram2d(
    x,
    y,
    bins=[xbins, ybins],
    weights = velocities,
    range=[[xbounds[0], xbounds[1]], [ybounds[0], ybounds[1]]]
    )

    counts, _, _ = np.histogram2d(
        x,
        y,
        bins=[xbins, ybins],
        range=[[xbounds[0], xbounds[1]], [ybounds[0], ybounds[1]]]
    )

    counts[counts == 0] = 1
    heatmap /= counts

    # heatmap = gaussian_filter(heatmap, sigma=0.5)

    plt.figure(figsize=(10,8))
    plt.imshow(
        heatmap.T,
        origin='lower',
        extent=[xbounds[0], xbounds[1], ybounds[0], ybounds[1]],
        cmap='magma',
        aspect='auto'
    )
    plt.colorbar(label='Average Velocity Magnitude')
    plt.title('Heatmap of Velocity Magnitudes')
    plt.xlabel('X Position')
    plt.ylabel('Y Position')
    FIGNAME = f"../tests/test_figures/heatmap_{xbins}_{ybins}.png"
    return plt.show()
