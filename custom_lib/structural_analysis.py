import numpy as np
import matplotlib.pyplot as plt
import os
import shutil
from scipy.spatial.distance import pdist
import matplotlib.pyplot as plt
import copy


def run_ATK_input(input_file, num_mpi):
	
	''' 
	Runs a ATK melt-quench process 
	
	Args: 
		1. ATKpython input file which can be used to run the process
		2. Number of mpi processes which should be used
		
	Returns:
		(None)
	'''

	os.system('mpiexec -n '+ str(num_mpi) + ' quantumatk-2019.12-kgf atkpython ' + str(input_file) + ' > output.out')		
	
def coord_to_xyz(filename):
	
	''' 
	Reads the coord file produced by atkpython, and converts it to an *.xyz file format
	
	Args: 
		1. *.coord file (atkpython output)
		
	Returns:
		(None)
	'''
	
	coord = []
	atoms = []
	lattice = np.empty([3])
	ind = 0

	with open(filename , "rt") as myfile:
		for line in myfile:

			# Primitive lattice vector lines:
			if len(line.split()) == 6:
				lattice[ind] = line.split()[2+ind]
				ind += 1
				
			# Cartesian coordinates
			if len(line.split()) == 7:
				try:
					coord.append(list(map(float, line.split()[1:4])))
					atoms.append(line.split()[0])
				except ValueError:
					continue
				
	num_atoms = np.shape(coord)[0]

	# Sort by the x-coordinate
	atom_coord = zip(atoms, coord)
	atom_coord = sorted(atom_coord, key=lambda x: x[1])

	# Assemble xyz file
	f = open('HfO2.xyz', 'w')
	f.write(str(num_atoms)+'\n')
	f.write('Cell:  '+str(lattice[0])+'  '+str(lattice[1])+'  '+str(lattice[2])+'\n')
	for x in atom_coord:
		f.write(str(x[0])+'\t'+str(x[1][0])+'\t'+str(x[1][1])+'\t'+str(x[1][2])+'\n')
	f.close()

def read_xyz(filename):
	
	''' 
	Reads data from xyz files
	
	Args:
		1. filename as '*.xyz'
		
	Returns:
		1. A (1x3) numpy array containing the (rectangular) unit cell size
		2. A list containing the atom symbols
		3. A (*x3) numpy array containing the atomic coordinates
	'''
	 
	atoms = []
	coords = []
	 
	with open(filename , "rt") as myfile:
		for line in myfile:
			if len(line.split()) == 1:
				pass
			elif line.split()[0] in ['Cell:', 'cell:']:
				lattice = line.split()[1:4]
			else:
				atoms.append(line.split()[0])
				coords.append(line.split()[1:])
	coords = np.asarray(coords, dtype=np.float64)
	lattice = np.asarray(lattice, dtype=np.float64)
	
	return np.array(lattice), np.array(atoms), np.array(coords)	
	

def write_xyz(lattice, coords, filename):
	
	''' 
	Writes an x-sorted xyz file with the given atomic number, lattice, and coordinate information
	
	Args:
		1. filename as '*.xyz'
		
	Returns:
		(None)
	'''
	
	# Sort the coords!!!
	print('sorting to be written')
	
	# Print the number of atoms on the first line
	num_atoms = np.shape(coords)[0]
	
	# Print the cell size on the second line	
	f = open(filename, 'w')
	f.writelines(str(num_atoms)+'\n')
	f.write('Cell:  '+"{:.10f}".format(lattice[0])\
	+'  '+"{:.10f}".format(lattice[1])\
	+'  '+"{:.10f}".format(lattice[2])+'\n')
	
	# Print the atomic name followed by coordinates on subsequent lines
	for line in coords:
		f.write(line[0]+'\t'\
		+"{:.10f}".format(float(line[1]))+'\t'\
		+"{:.10f}".format(float(line[2]))+'\t'\
		+"{:.10f}".format(float(line[3]))+'\n')

	f.close()
	
	
def update_xyz_coordinates(old_file, new_file):
	
	''' 
	Updates the xyz coordinates in 'old_file' with the last N xyz coordinates in 'new-file', where N is the number of atoms in the structure.
	Can be used to update existing atomic positions with geometry optimization runs from cp2k.
	
	Args:
		1. (str) Name of old xyz file to be replaced
		2. (str) Name of new xyz file to be replaced
		
	Returns:
		(none)
	'''
	
	f = open(old_file, 'r')
	structure_lines = f.readlines()
	num_atoms = int(structure_lines[0][0:]) 
	cell_size = structure_lines[1]
	f.close()
	
	f = open('pre_update.xyz', 'w')
	f.writelines(structure_lines)
	f.close
	
	f = open(new_file, 'r')
	structure_lines = f.readlines()
	new_coords=structure_lines[-num_atoms:]
	f.close()
	
	f = open(old_file, 'w')
	f.writelines(str(num_atoms)+'\n')
	f.writelines(cell_size)
	f.writelines(new_coords)
	f.close
	
def find_coordination(lattice, atoms, coords, from_atom, to_atom, bond_length):
	
	''' 
		!!! TO BE FIXED !!! 
		
	Calculates the atomic coordination numbers of a given periodic lattice
	
	Args: 
		1. A (*x3) numpy array containing the atomic lattice
		2. A list containing the atom symbols
		3. A (*x3) numpy array containing the atomic coordinates
		4. The atom symbol of the atomic centers
		5. The atom symbol of the ligands
		6. The expected bond length (in Angstrom), below which a bond will be counted
		
	Returns:
		1. A (1xX) numpy array containing the coordination numbers of every type of 'from_atom' in the atomic coordinates array 
	'''
	
	from_atom_list = [i for i in range(len(coords)) if atoms[i] == from_atom]
	to_atom_list = [i for i in range(len(coords)) if atoms[i] == to_atom]
	coordination_numbers = np.zeros(len(coords))
	
	for i in from_atom_list:
		for j in to_atom_list:
			
			to_atom_pbc1 = coords[j]
			to_atom_pbc2 = coords[j] - lattice
			
			min_x = np.minimum(np.absolute(coords[i, 0]-to_atom_pbc1[0]),  np.absolute(coords[i, 0]-to_atom_pbc2[0]))
			min_y = np.minimum(np.absolute(coords[i, 1]-to_atom_pbc1[1]),  np.absolute(coords[i, 1]-to_atom_pbc2[1]))
			min_z = np.minimum(np.absolute(coords[i, 2]-to_atom_pbc1[2]),  np.absolute(coords[i, 2]-to_atom_pbc2[2]))
			
			dist = np.sqrt(min_x**2 + min_y**2 + min_z**2)
			
			if dist < bond_length:
				coordination_numbers[i] += 1
		
	coordination_numbers = coordination_numbers[from_atom_list]
	
	plt.hist(coordination_numbers, bins = [3, 4, 5, 6, 7, 8])
	plt.xlabel(from_atom+' Coordination')
	plt.ylabel('Counts (#)')
	plt.ylim(0, len(from_atom_list))
	plt.savefig("coordination_numbers.png")	
		
	return coordination_numbers

def find_voids(lattice, atoms, coords, from_atom, to_atom, bond_length):
	print("To be written!")
	pass


def make_contacts(lattice, metal_cell, atom_positions, atom_names, interface_spacing, thickness):
	
	''' 
	When given information about a 'channel' lattice, attaches metal 'contacts' to either side of it. The contact material can be fully specified
	by the unit cell, atomic positions, and atom names, and the contact is created by tiling it according to the size of the channel and the thickness specified.
	
	Args:
		1. (3x3 numpy array) The rectangular atomic lattice
		2. (3x3 numpy array) The smallest rectangular unit cell of the contact metal
		3. (*x3 numpy array) The coordinates of all the atoms in the channel
		4. (*x3 list) The atom names of all the atoms in the channel
		5. (numeric float) The desired spacing between the contacts and the channel
		6. (numeric float) The desired thickness of the contacts in the transport direction 
		
	Returns:
		1. (*x4) List of atom names (first column) and coordinates (2nd-4th columns) for the left metal contact
		2. (*x4) List of atom names (first column) and coordinates (2nd-4th columns) for the right metal contact
	'''
	
	#Find the size of the metal contacts, constrained by the size of the oxide (in the xy directions), and the specified thickness (in z direction)
	z_repetitions = int(np.ceil(lattice[0]/metal_cell[0, 0]))
	y_repetitions = int(np.ceil(lattice[1]/metal_cell[1, 1]))
	x_repetitions = int(np.ceil(thickness/metal_cell[2, 2]))
	tiling = [x_repetitions, y_repetitions, z_repetitions]
	
	output_positions = copy.deepcopy(atom_positions)
	output_atoms = copy.deepcopy(atom_names)
	metal_supercell = []
	
	# Tile cell in each dimension
	for dim, direction in enumerate(tiling):
	
		for i in range(direction-1):
			output_positions = np.append(output_positions, atom_positions+(i+1)*metal_cell[dim], axis=0)
			output_atoms = np.append(output_atoms, atom_names, axis=0)
		
		atom_positions = copy.deepcopy(output_positions)
		atom_names = copy.deepcopy(output_atoms)
		metal_supercell.append(metal_cell[dim] * direction)
	
	metal_supercell = np.array(metal_supercell)
		
	# Reduce strain (?) -- Ask Marko about this
	# for row in output_positions:
	# 	row = np.divide(row, metal_supercell * np.diag(lattice))
	
	# Make the left contact
	shift_left = -interface_spacing
	output_positions_left = copy.deepcopy(output_positions)
	output_positions_left[:,0] = output_positions_left[:,0]*-1 + shift_left
	output_atoms_left = np.flip(atom_names)
	left_contact = np.concatenate((output_atoms_left[:, None], output_positions_left), axis = 1)
	
	# Make the right contact
	shift_right = lattice[0] + interface_spacing
	output_positions_right = copy.deepcopy(output_positions)
	output_positions_right[:, 0] = output_positions[:, 0] + shift_right
	output_atoms_right = atom_names	
	right_contact = np.concatenate((output_atoms_right[:, None], output_positions_right), axis = 1)
	
	# Sort the atom lines by the x-coordinate (the transport direction)
	left_contact = left_contact[np.argsort(np.asarray(left_contact[:, 1], dtype=np.float64, order='C'))]
	right_contact = right_contact[np.argsort(np.asarray(right_contact[:, 1], dtype=np.float64, order='C'))]

	return(left_contact, right_contact, metal_supercell)
