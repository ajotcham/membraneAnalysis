import gzip
import glob
import numpy as np
import natsort

files = natsort.natsorted(glob.glob('output/dump_md.*.lammpstrj'))
"""boundaries are in 3D x,y,z"""

def read_atoms(file):
    """
    Passes one file through and returns arrays xyz pos, xyz velo, atom# & type, #atoms in sim, and boundaries
    """
    # for filename in files:
    timestep = 0
    atoms = 0
    boundaries = np.zeros((3, 2))
    split_file = file.split('.')
    abs_path = split_file[0]
    prefix = abs_path.split('/')
    # print(prefix[-1])
    file_extension = split_file[-1]
    if file_extension == 'gz':
        open_file = gzip.open(file, 'rt')
    else:
        open_file = open(file, 'r')

    lines = open_file.readlines()

    # for n_line, line in enumerate(lines):
    #     #Removes '\n' string in found in lines
    #     lines = [line.rstrip() for line in file]
    #     if n_line < 9:
    #         continue
    headers = lines[0:9]
    data = lines[9:]
    # print(headers)
    # print(data)
    
    for c, line in enumerate(headers):
        # print(c)
        if c == 1:
            timestep = int(line)
            #timestep = int(lines[0])
            # print(timestep)

        elif c == 3:
            num_atoms = int(line)
            # print(num_atoms)
 
        elif c >= 5 and c <=7:
            for n,n_l in enumerate([0,1]):
                boundaries[c-5, n_l] = float(line.split(" ")[n_l])

    atoms = {
        "positions": np.zeros((num_atoms, 3), dtype=np.float64),
        "mass": np.zeros((num_atoms, 1), dtype=np.float64),
        "velocities": np.zeros((num_atoms, 3), dtype=np.float64),
        "map": np.zeros((num_atoms, 2), dtype=np.int64),
        "peratomvol": np.zeros((num_atoms, 1), dtype=np.float64)
    }
    # atom_pos = np.zeros((num_atoms, 3),dtype=np.float64)
    # atom_velo = np.zeros((num_atoms, 3),dtype=np.float64)
    # atom_map = np.zeros((num_atoms,2),dtype=np.int64)

    for c, data in enumerate(data):
        for n,n_l in enumerate([0,1]):
            atoms["map"][c,n_l] = int(data.split(" ")[n_l])
        
        atoms['mass'][c] = float(data.split(" ")[3])
        
        # for n,n_l in enumerate([2,3,4]):
        for n, n_l in enumerate([5,6,7]):
            atoms["positions"][c,n] = float(data.split(" ")[n_l])
        
        atoms['peratomvol'][c] = float(data.split(" ")[9])

        if prefix[-1] != 'membranedata':
            # for n,n_l in enumerate([5,6,7]):  
            for n,n_l in enumerate([11,12,13]): 
                atoms["velocities"][c,n] = float(data.split(" ")[n_l])
            
            atoms['peratomvol'][c] = float(data.split(" ")[15])


        
        else:
            continue

    # print(atoms)
    return atoms, num_atoms, boundaries


# def poiseuille_subregion_atoms(atom_pos, atom_velo, atom_map, x_range, y_range, z_range):
#     """
#     Filters data to specific sub_region
#     x_range = tuple, (min_x, max_x) for desired subdomain
#     y_range = tuple, (min_y, max_y) for desired subdomain
#     z_range = tuple, (min_z, max_z) for desired subdomain
#     """
#     mask = (
#         (atom_pos[:, 0] >= x_range[0]) & (atom_pos[:, 0] <= x_range[1]) &
#         (atom_pos[:, 1] >= y_range[0]) & (atom_pos[:, 1] <= y_range[1]) &
#         (atom_pos[:, 2] >= z_range[0]) & (atom_pos[:, 2] <= z_range[1])
#     )

#     sub_atoms = {
#         "positions": atom_pos[mask],
#         "velocities": atom_velo[mask],
#         "map": atom_map[mask]
#     }

#     return sub_atoms