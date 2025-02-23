"""Plot aggregate density of water/membrane atoms and total"""

import mypackage
import numpy as np
import glob as glob
import natsort
import os
import matplotlib.pyplot as plt


# Path to the folder containing the files
folder_path = "/Users/ajotcham/Desktop/rileydata_mem_and_water"

# Define the file pattern to match specific files
membrane_file_pattern = os.path.join(folder_path, "membranedata.*.gz")
water_file_pattern = os.path.join(folder_path, "waterdata.*.gz")

# List all files matching the pattern
all_water_files = natsort.natsorted(glob.glob(water_file_pattern))
all_membrane_files = natsort.natsorted(glob.glob(membrane_file_pattern))

# Select specific files (500th to 550th)
selected_membrane_files = all_membrane_files
selected_water_files = all_water_files

all_water_dicts = []
all_membrane_dicts = []
# #Loop through desired timesteps and append timestep specific data to larger dataset
for file1, file2 in zip(all_membrane_files, all_water_files):
    print(file1)
    membrane_atoms, num_membrane_atoms, membrane_boundaries = mypackage.read_atoms(
        file1
    )
    print(file2)
    water_atoms, num_water_atoms, water_boundaries = mypackage.read_atoms(file2)

    all_membrane_dicts.append(membrane_atoms)
    all_water_dicts.append(water_atoms)

# water_atoms and mem_atoms are dictionaries with keys "positions", "velocities", "map", "mass", "peratomvol"
combined_membrane_dict = {}
combined_water_dict = {}
keys_to_combine = ["positions", "velocities", "map", "mass", "peratomvol"]
for key in keys_to_combine:
    combined_membrane_dict[key] = np.concatenate(
        [d[key] for d in all_membrane_dicts], axis=0
    )

    combined_water_dict[key] = np.concatenate([d[key] for d in all_water_dicts])


def get_bin_edges(min_z, max_z, n_bins):
    """Calculate bin_edges based on min, max and number of bins wanted"""
    bin_edges = np.linspace(min_z, max_z, n_bins + 1)
    return bin_edges


def get_bin_centers(bin_edges):
    """
    Calculate bin centers from bin edges.
    """

    return (bin_edges[:-1] + bin_edges[1:]) / 2


def get_volume(x, y, z):
    """Calculate volume"""
    return x * y * z


def bin_atoms_sum_mass(atoms_dict, bin_edges, key="positions", mass_key="mass"):
    """
    Bins atoms based on their z-positions and sums their mass for each bin

    Parameters:
        atoms_dict (dict): Dictionary containing atoms data
        bin_edges (ndarraY): array of bin edges along z-axis
        key (str): key for atoms_dict for positions
        mass_key (str): key in atoms_dict for masses

    Returns:
        atom_counts (ndarray): number of atoms in each bin
        mass_sums (ndarray): mass per bin of atoms from dict
    """

    n_bins = len(bin_edges) - 1
    atom_counts = np.zeros(n_bins, dtype=int)
    mass_sums = np.zeros(n_bins, dtype=float)

    if key in atoms_dict and mass_key in atoms_dict:
        atom_pos = atoms_dict[key]
        atom_masses = atoms_dict[mass_key]

        for i in range(len(atom_pos)):
            z_value = atom_pos[i, 2]  # extract z-coordinate from pos array
            mass = float(atom_masses[i].item())
            bin_index = np.digitize(z_value, bin_edges) - 1
            if 0 <= bin_index < n_bins:
                atom_counts[bin_index] += 1
                mass_sums[bin_index] += mass

    return atom_counts, mass_sums


def calculate_density(bin_edges, mass_sums, xdim, ydim):
    """
    Calculate density of atoms per bin

    Parameters:
        bins_edges (ndarray): bin edges on z-axis
        mass_sums (ndarray): summed mass per bin
        xdim (float): length of box in x direction
        ydim (float): length of box in y direction
        Cross-sectional average of x-y for this case and only plotting against z-axis

    Returns:
        bin_centers (ndarray): center of each bin on z-coordinate for plotting
        density (ndarray): density in each bin (g/cm^3)
    """

    avogadro = 6.022e23  # Avogadro's number
    angstrom_to_cm = 1e-8
    bin_widths = np.diff(bin_edges)  # bin width in angstroms

    bin_volumes = get_volume(xdim, ydim, bin_widths) * (angstrom_to_cm**3)

    # Convert mass to grams (input mass is in g/mol)
    mass_grams = mass_sums / avogadro

    # calculate density (g/cm^3)
    density = mass_grams / bin_volumes

    return density


USE_FUNCTIONS = True
if USE_FUNCTIONS:
    if len(all_water_files) == len(all_membrane_files):
        num_timesteps = len(
            all_water_files
        )  # assuming number of water/membrane files are the same
    else:
        pass
    n_bins = 300
    min_z, max_z = -287, 237  # z dimensions of simualted domain
    xdim = ydim = 175
    # get bin edges and bin centers
    bin_edges = get_bin_edges(min_z, max_z, n_bins)
    bin_centers = get_bin_centers(bin_edges)

    # get water and membrane atom counts and masses per bin
    water_counts, water_mass_sum = bin_atoms_sum_mass(combined_water_dict, bin_edges)
    membrane_counts, membrane_mass_sum = bin_atoms_sum_mass(
        combined_membrane_dict, bin_edges
    )

    # calculate average water and membrane mass per bin over all timesteps
    average_water_mass = water_mass_sum / num_timesteps
    average_membrane_mass = membrane_mass_sum / num_timesteps

    # calclate average density per bin for water and membrane atoms
    water_density = calculate_density(bin_edges, average_water_mass, xdim, ydim)
    membrane_density = calculate_density(bin_edges, average_membrane_mass, xdim, ydim)
    total_density = water_density + membrane_density


# Plot avg density of all timesteps of water/membrane and total.
PLOT_AVG_DENSITY = True
if PLOT_AVG_DENSITY:

    def plot_density(bin_centers, water_density, membrane_density, total_density):
        """
        Plot density of water and membrane atoms in polyamide membrane.
        """

        # Plot Water Atoms Density
        plt.figure(figsize=(10, 8))
        plt.plot(
            bin_centers,
            water_density,
            linestyle="-",
            color="dodgerblue",
            label="Water Density",
        )

        # Plot Membrane Atoms Density
        plt.plot(
            bin_centers,
            membrane_density,
            linestyle="-",
            color="seagreen",
            label="Membrane Density",
        )

        # Plot Total Density
        plt.plot(
            bin_centers,
            total_density,
            linestyle="-",
            color="slategrey",
            label="Total Density",
        )

        plt.title("Density Profile of Polyamide RO Membrane")
        plt.xlabel("Direction of Flow on Z-axis Position (Å)")
        plt.ylabel("Density (g/cm³)")
        plt.legend()
        plt.grid(True)

        return plt.savefig("../../../../../average_density_profile.png")

    plot_density(bin_centers, water_density, membrane_density, total_density)
