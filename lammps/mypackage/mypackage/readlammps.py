"""glob provides file sorting"""
import glob
import numpy as np
import natsort

files = natsort.natsorted(glob.glob('output/dump_md.*.lammpstrj'))
"""boundaries are in 3D x,y,z"""

def process_atoms(file):
    """
    Passes one file through and returns arrays xyz pos, xyz velo, atom# & type, #atoms in sim, and boundaries
    """
    # for filename in files:
    timestep = 0
    atoms = 0
    boundaries = np.zeros((3, 2))
    with open(file, 'r', encoding="utf-8") as file:
        #Removes '\n' string in found in lines
        lines = [line.rstrip() for line in file]
        headers = lines[0:9]
        data = lines[9:]

    for c, line in enumerate(headers):
        if c == 1:
            timestep = int(line)
            #timestep = int(lines[0])
            # print(timestep)

        elif c == 3:
            num_atoms = int(line)
            # print(atoms)
 
        elif c >= 5 and c <=7:
            for n,n_l in enumerate([0,1]):
                boundaries[c-5, n_l] = float(line.split(" ")[n_l])

    atom_pos = np.zeros((num_atoms, 3),dtype=np.float64)
    atom_velo = np.zeros((num_atoms, 3),dtype=np.float64)
    atom_map = np.zeros((num_atoms,2),dtype=np.int64)

    for c, data in enumerate(data):
        for n,n_l in enumerate([0,1]):
            atom_map[c,n_l] = int(data.split(" ")[n_l])

        for n,n_l in enumerate([2,3,4]):
            atom_pos[c,n] = float(data.split(" ")[n_l])
        
        for n,n_l in enumerate([5,6,7]):
            atom_velo[c,n] = float(data.split(" ")[n_l])

    return atom_pos, atom_velo, atom_map, num_atoms, boundaries


# def test_function(test_path):
#     """
#     test for process_atoms function
#     """
#     with open(test_path, 'r') as file:
#         lines = [line.rstrip() for line in file]
#         headers = lines[0:9]
#         data = lines[9:]

#     for c, lines in enumerate(headers):
#             if c == 2:
#                 timestep = int(headers[1])
#                 print(timestep)

#             if c == 4:
#                 atoms = int(headers[3])
#                 print(atoms)

#             if c == 6:
#                 #Can be rewritten into for loop
#                 continue

#             else:
#                 continue
                
#             test_info = np.zeros((atoms, 8))
#             for c, data in enumerate(data):
#                 for n,n_l in enumerate([0,1,2,3,4,5,6,7]):
#                     test_info[c,n_l] = float(data.split(" ")[n_l])
#     return test_info

# process_atoms(files)
# print(boundaries)

