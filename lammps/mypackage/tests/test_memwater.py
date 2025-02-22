import mypackage
import numpy as np
import glob as glob
import natsort
import os
import matplotlib.pyplot as plt

# Path to the folder containing the files
folder_path = "/Users/ajotcham/Desktop/rileydata_mem_and_water"

# Define the file pattern to match specific files
membrane_file_pattern = os.path.join(folder_path, "membranedata.90000000.gz")
water_file_pattern = os.path.join(folder_path, "waterdata.90000000.gz")

# List all files matching the pattern
all_water_files = natsort.natsorted(glob.glob(water_file_pattern))
all_membrane_files = natsort.natsorted(glob.glob(membrane_file_pattern))

# print(all_water_files)
# print(all_membrane_files)
# Select specific files (500th to 550th)
selected_membrane_files = all_membrane_files
selected_water_files = all_water_files
all_membrane_positions = []
all_water_positions = []
# all_velocities_mag = []
# #Loop through desired timesteps and append timestep specific data to larger dataset
for file1, file2 in zip(all_membrane_files, all_water_files):
    print(file1)
    mem_atoms, num_mem_atoms, mem_boundaries = mypackage.read_atoms(file1)
    print(file2)
    water_atoms, num_water_atoms, water_boundaries = mypackage.read_atoms(file2)

    # all_positions.append(atoms["positions"])
    # mag_velocity = mypackage.calc_velo(atoms["velocities"])
    # all_velocities_mag.append(mag_velocity)
# print(mem_atoms['mass'])
# print(water_atoms['mass'])
# all_positions = np.concatenate(all_positions, axis=0)
# all_velocities_mag = np.concatenate(all_velocities_mag, axis=0)
# top_atoms, top_indices = mypackage.top_velo(all_velocities_mag, 1)
# np.savetxt('test_membranedata.csv', all_positions, delimiter=' ')

###Calculate density, assume volume = 175^3
# Query water_atoms and based on their positions determine if inside box or outside box.
# Only consider atoms inside box x = 0-175, y = 0-175, z = 0-175 for this case


def bin_atoms(atoms_dict, n_bins, min_z, max_z, key="positions", return_type="all"):
    """
    Filters atoms in water_atoms based on their z position to find if they are inside box or not.
    Parameters
    water_dict(dict) = dictionary for water atoms in system
    key(string) = key for dictionary should be positions
    min z(float) = minimum z value
    max z(float) = maximum z value
    """
    bin_edges = np.linspace(min_z, max_z, n_bins + 1)
    atom_counts = np.zeros(n_bins, dtype=int)
    # indices = []

    def count_atoms_in_bin(atom_positions, bin_edges, counts_array):
        """Bin atoms based on z coordinate"""
        for z_value in atom_positions[:, 2]:  # Extract z-coordinate
            bin_index = np.digitize(z_value, bin_edges) - 1
            if 0 <= bin_index < n_bins:
                counts_array[bin_index] += 1

    if key in atoms_dict:
        count_atoms_in_bin(atoms_dict[key], bin_edges, atom_counts)

    if return_type == "counts":
        return atom_counts
    elif return_type == "edges":
        return bin_edges
    else:
        return bin_edges, atom_counts


def bin_atoms_with_mass(
    atoms_dict, n_bins, min_z, max_z, key="positions", mass_key="mass"
):
    """
    Bins atoms based on their z-position and sums their mass in each bin.

    Parameters:
        atoms_dict (dict): Dictionary containing atomic data (positions, mass).
        n_bins (int): Number of bins.
        min_z (float): Minimum z-coordinate.
        max_z (float): Maximum z-coordinate.
        key (str): Key in atoms_dict for atomic positions.
        mass_key (str): Key in atoms_dict for atomic masses.

    Returns:
        bin_edges (ndarray): Bin edges along the z-axis.
        atom_counts (ndarray): Number of atoms per bin.
        mass_sums (ndarray): Total mass per bin.
    """
    bin_edges = np.linspace(min_z, max_z, n_bins + 1)
    atom_counts = np.zeros(n_bins, dtype=int)
    mass_sums = np.zeros(n_bins, dtype=float)

    if key in atoms_dict and mass_key in atoms_dict:
        atom_positions = atoms_dict[key]
        atom_masses = atoms_dict[mass_key]

        for i in range(len(atom_positions)):
            z_value = atom_positions[i, 2]  # Extract z-coordinate
            mass = float(atom_masses[i])  # Ensure scalar value
            bin_index = np.digitize(z_value, bin_edges) - 1
            if 0 <= bin_index < n_bins:
                atom_counts[bin_index] += 1
                mass_sums[bin_index] += mass  # Sum mass in each bin

    return bin_edges, atom_counts, mass_sums


# selected_indices = query_water_atoms(water_atoms, 0.0, 175.0) #indices where masses are summed
# print(water_atoms['mass'])


def sum_water_mass(water_dict, indices, key="mass"):
    "Sum water masses at given indices"
    # print(len(water_dict[key][0]))
    mass_water = 0
    for i in indices:
        mass_water += water_dict[key][i]
    # mass_water = sum(water_dict[key][i] for i in indices)
    return mass_water


def sum_water_vol(water_dict, indices, key="peratomvol"):
    "Sum water masses at given indices"
    # print(len(water_dict[key][0]))
    vol_water = 0
    for i in indices:
        vol_water += water_dict[key][i]
    # mass_water = sum(water_dict[key][i] for i in indices)
    return vol_water


def calculate_density(bins, mass_sums, xdim, ydim):
    """
    Calculates the density of water atoms per bin.

    Parameters:
        bins (ndarray): Bin edges along the z-axis.
        water_mass_sums (ndarray): Summed mass per bin.
        x_dim (float): Length of the box in the x-direction.
        y_dim (float): Length of the box in the y-direction.

    Returns:
        bin_centers (ndarray): Center of each bin (z-coordinate).
        density (ndarray): Density in each bin (g/cm³).
    """

    avogadro = 6.022e23  # Avogadro's number
    angstrom_to_cm = 1e-8  # Convert angstrom to cm
    bin_widths = np.diff(bins)  # Bin width in angstroms
    bin_centers = (bins[:-1] + bins[1:]) / 2  # Compute bin centers

    bin_volumes = bin_widths * x_dim * y_dim * (angstrom_to_cm**3)

    # Convert mass to grams (assuming input mass is in g/mol, divide by Avogadro)
    mass_grams = mass_sums / avogadro  # Convert to grams

    # Calculate density (g/cm³)
    density = mass_grams / bin_volumes

    return bin_centers, density


PLOT_ATOM_COUNTS = True
if PLOT_ATOM_COUNTS:

    def plot_atom_distribution(bin_edges, water_counts, mem_counts):
        """
        Plot number of atoms (water or membrane) per bin as a function of z-axis
        """
        bin_edge_color = "darkslategray"
        water_color = "dodgerblue"
        mem_color = "mediumpurple"
        bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2

        plt.figure(figsize=(8, 6))
        # Plot water atoms
        plt.plot(
            bin_centers,
            water_counts,
            marker="o",
            linestyle="-",
            color=water_color,
            label="Water",
        )

        # Plot membrane atoms
        plt.plot(
            bin_centers,
            mem_counts,
            marker="D",
            linestyle="-",
            color=mem_color,
            label="Membrane",
        )

        ymax = int(water_counts.max().item())
        plt.vlines(
            bin_edges,
            ymin=0,
            ymax=ymax,
            color=bin_edge_color,
            linestyle="--",
            alpha=0.5,
        )

        plt.xlabel("Z-Position (Å)")
        plt.ylabel("Number of atoms")
        plt.title("Atom Distribution Along Z-axis")
        plt.legend()
        plt.grid(False)

        return plt.savefig("../../../../../atom_distribution.png")


# mass_water_in_box = sum_water_mass(water_atoms, selected_indices) #mass returned in g/mol
# print(mass_water_in_box)
# avogadro = 6.023e23
# angstromCube_cmCube = (1e-8)**3
# volume = (175*175*175) #volume in angstroms
# density = (mass_water_in_box) / (total_vol * angstromCube_cmCube * avogadro)
# print(f"Mass water in box {mass_water_in_box}")
# water_vol = sum_water_vol(water_atoms, selected_indices)
# print(f"Water volume {water_vol}")
# mem_vol = np.sum(mem_atoms['peratomvol'])
# print(f"Membrane volume {mem_vol}")
# total_vol = water_vol + mem_vol
# print(f"Total Peratom Volume {total_vol}")
# print(f"Volume {volume}")

n_bins = 100
min_z, max_z = -287, 237
# bins, water_counts = bin_atoms(water_atoms, n_bins, min_z, max_z)
# mem_counts = bin_atoms(mem_atoms, 100, -287, 237, return_type="counts")

bins, water_counts, sum_water_mass = bin_atoms_with_mass(
    water_atoms, n_bins, min_z, max_z
)
# print(bins)

bins, mem_counts, sum_mem_mass = bin_atoms_with_mass(mem_atoms, n_bins, min_z, max_z)

print(bins)

x_dim = y_dim = 175
bin_centers, water_density = calculate_density(bins, sum_water_mass, x_dim, y_dim)
print(water_density)
bin_centers, membrane_density = calculate_density(bins, sum_mem_mass, x_dim, y_dim)
print(membrane_density)
# print(water_mass)
# plot = plot_atom_distribution(bins, water_counts, mem_counts)

bin_edges = np.linspace(min_z, max_z, n_bins + 1)
PLOT_DENSITY = True
if PLOT_DENSITY:

    def plot_density_profile(bin_centers, bin_edges, water_density, mem_density):
        """
        Plots the density of water atoms along the z-axis.

        Parameters:
            bin_centers (ndarray): Z-coordinate centers of bins.
            density (ndarray): Density of water atoms per bin.
        """
        plt.figure(figsize=(10, 8))
        plt.plot(
            bin_centers,
            water_density,
            marker="o",
            linestyle="-",
            color="dodgerblue",
            label="Water Density",
        )

        plt.plot(
            bin_centers,
            mem_density,
            marker="D",
            linestyle="-",
            color="seagreen",
            label="Membrane Density",
        )

        # plt.vlines(
        #     bin_edges,
        #     ymin=0,
        #     ymax=max(np.max(water_density), np.max(mem_density)),
        #     colors="gray",
        #     linestyles="dashed",
        #     alpha=0.5,
        # )

        plt.xlabel("Z-Position (Å)")
        plt.ylabel("Density (g/cm³)")
        plt.title("Water and Membrane Density Profile Along Z-axis")
        plt.legend()
        plt.grid(True)

        return plt.savefig("../../../../../density_profile_with_grid.png")

    plot_density_profile(bin_centers, bin_edges, water_density, membrane_density)
