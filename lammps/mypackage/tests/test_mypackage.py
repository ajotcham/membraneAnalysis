import mypackage
import numpy as np
import matplotlib

file = 'test_dump.lammpstrj'

atom_pos, atom_velo, atom_map, num_atoms, boundaries = mypackage.process_atoms(file)

mag_velocity = mypackage.calc_velo(atom_velo)

top_atoms, top_indices = mypackage.top_velo(mag_velocity, 10)

# print(atom_pos)
# print(mag_velocity)
# print(top_atoms)
# print(top_indices)
# print(boundaries)
# print(num_atoms)

#create plot
bounds = [-20,20,-20,20]
mypackage.plot_trajectory(atom_pos, bounds)
mypackage.plot_individual_location(top_indices, atom_pos)

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

# print(x)
# print(y)
# print(centerpoints)
# print(lengths)
# print(deltas)
# print(bounds)
print(boundaries)
print(top_indices)
for index in top_indices:
    ind = mypackage.find_grid_cell(atom_pos[index, :], boundaries, discretization, lengths)
    print(ind)
    grid[ind[0], ind[1]] += 1. #increments grid by one 

print(grid)
print(lengths)