import os
import struct
import numpy as np
import matplotlib.pyplot as plt
from scipy import sparse
from copy import deepcopy


def extract_basis_set(output_log):
	''' 
	Reads the atomic kinds (Hf, Ti, etc) and # orbitals per atom from a cp2k output log
	
	Args:
		1. cp2k output filename to search
		
	Returns:
		1. A (1x* str) list of atomic symbols, where * is the number of different atomic types present
		2. A (1x* int) list of the # of orbitals corresponding to each atom in the atomic_kinds list
	'''
	
	atomic_kinds = []
	no_orbitals = []
	
	# Look for the atomic kinds 
	with open(output_log, 'r') as f:
		for line in f:
			if 'Atomic kind:' in line:
				atomic_kinds.append(line.split()[3])
			if 'Number of spherical basis functions:' in line:
				no_orbitals.append(int(line.split()[5]))

	return atomic_kinds, no_orbitals
	
def get_value_from_file(filename, keywords):
	''' 
	Reads the file given by filename, and finds the first number on every line where the keyword(s) is/are mentioned
	
	Args:
		1. A human-readable filename to search
		
	Returns:
		1. The first number in every line which contains the keyword(s)
	'''
	
	found_lines = []
	
	# Look for the atomic kinds 
	with open(output_log, 'r') as f:
		for line in f:
			if keyword in line:
				for x in line.split():
					if x.isdigit():
						found_values.append(x)
				
	return found_values

def read_bin(binfile, struct_fmt='<IIIdI'):
	''' 
	Reads the .bin files produced by cp2k for the KS and S matrices
	one data set containes: uint32, uint32, uint32, real*8, uint32,
	corresponding to	  : header, x index, y index, data, header
	NOTE: The bytes are written in little endian in the cp2k output files
	
	Args:
		1. (str) name of the .bin file to read
		2. (str) Format in which the .bin file was written. Use a leading '<' for little endian.
		
	Returns:
		1. (*x*) a N-dim list of the .bin contents, with types corresponding to the 'struct_fmt' arg
			The headers of each line is ommitted in the returned matrix.
	'''
	
	# Get the structure bytelengths
	struct_len_I = struct.calcsize(struct_fmt)
	struct_unpack_I = struct.Struct(struct_fmt).unpack_from
	
	# Initialize the returned matrix according to the file size
	filesize = os.path.getsize(binfile)
	num_interaction = int(filesize/struct_len_I)
	M = np.zeros((num_interaction, 3))
	
	with open(binfile, "rb") as f:
		for ind in range(num_interaction):
			data = f.read(struct_len_I)
			
			# Ignore headers when extracting data:
			M[ind] = struct_unpack_I(data)[1:-1]

	return M
	
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
	
	return lattice, np.array(atoms), coords	
	
def bin_to_csr(M):

	''' 
	Converting from
		 .bin format (x-index, y-index, interaction_parameter)
		  --> matrix format (interaction_parameter, interaction_parameter)
	Args:
		1. (*x3) a N-dim list of the .bin contents, as returned by read_bin() 
		for struct_fmt '<IIIdI'
		
	Returns:
		2. (*x*) A csr matrix created by mapping the third row of the input 
			with the indices specified by the first two rows 
	'''
	
	# The size of the matrices is the largest atomic index
	xmax, ymax, zmax = M.max(axis=0) 
	
	# The first row of the data in the .bin file is the interaction parameters
	interaction_parameters = M[:,2]
	
	# The second and third rows are the (x, y) indices
	x_indices = M.astype(int)[:,0]-1
	y_indices = M.astype(int)[:,1]-1
	
	# Assembly the bin files into sparse matrix representations
	N = sparse.csr_matrix((interaction_parameters, (x_indices, y_indices)), shape=(int(xmax), int(ymax))).toarray()
	#Note: Save as upper triangular matrix later
	
	return N
	
def create_device_matrix(M, coords, num_orb_per_atom, blocks, no_atoms_per_block, no_orb, block_delete, block_repeat, eps):
	
	''' 
	Internal fxn. Takes the KS and S matrices and modifies them accordingly to create the device matrices,
	by adding/removing contact blocks, formatting the channel elements, etc.
	
	@Manasa: clean up this function it looks atrocious.
	
	Args:
		
	Returns:
	'''
	M = np.array(M)
	size_M = np.shape(M)[0]

	# Filter the matrix by removing small elements
	M[np.absolute(M) < eps] = 0
	
	# Remove periodic boundary conditions:
	bnd = int(np.round(size_M/2))		
	M_pbc = np.triu(M[:bnd, bnd:])
	M[:bnd, bnd:] = M[:bnd, bnd:] - M_pbc
	
	# Get the total number of orbital functions per contact block
	num_orb_left = np.sum(num_orb_per_atom[:no_atoms_per_block[0]])
	num_orb_right = np.sum(num_orb_per_atom[len(num_orb_per_atom)-no_atoms_per_block[1]:])
	
	# Delete blocks as required
	if block_delete[0]+block_delete[1] != 0:
		print('Note: deleting blocks is not yet implemented in this code.')
		
	# Get the left and right bandwidths of the input matrix
	no_neigh_left = 0	
	for ind1 in range(1, int(np.floor(size_M/num_orb_left))-1):
		if np.count_nonzero(M[:num_orb_left, ind1*num_orb_left:(ind1+1)*num_orb_left]) > 0:
			no_neigh_left += 1
		else:
			break
			
	no_neigh_right = 0	
	for ind1 in range(1, int(np.floor(size_M/num_orb_right))-1):
		if np.count_nonzero(M[size_M-(ind1+1)*num_orb_right:size_M-ind1*num_orb_right, size_M-num_orb_right:size_M]) > 0:
			no_neigh_right += 1
		else:
			break
	
	# Repeat contact blocks
	# @Manasa: remember to implement the pbc conditions for contact block repeats later...
	
	# Average of the second and second last contact atomic positions - gives us the x-spacing in the contacts
	d = [np.mean(coords[no_atoms_per_block[0]:2*no_atoms_per_block[0], 0] - coords[:no_atoms_per_block[0], 0]),\
	 np.mean(coords[-no_atoms_per_block[1]:, 0] - coords[-2*no_atoms_per_block[1]:-no_atoms_per_block[1], 0])]
	 

	# **************************
	# Repeating the left contact
	# **************************
	
	left_repeats = block_repeat[0]
	left_contact_block = M[:num_orb_left, :num_orb_left]
	
	Mnew = np.concatenate((np.kron(np.identity(left_repeats), left_contact_block),\
		   np.zeros((num_orb_left*left_repeats, num_orb_left*no_neigh_left))), axis=1)

	Medge = M[:num_orb_left, num_orb_left:num_orb_left*(no_neigh_left+1)]
	
	for ind in range(1, no_neigh_left+1):
		Mnew = Mnew + np.concatenate((np.zeros((num_orb_left*left_repeats, num_orb_left*ind)),\
									  np.kron(np.identity(left_repeats), Medge[:, (ind-1)*num_orb_left:ind*num_orb_left]),\
									  np.zeros((num_orb_left*left_repeats, num_orb_left*(no_neigh_left-ind)))), axis=1)
		
	# Add the left contact blocks to the main matrix:
	M_top = np.concatenate((Mnew, np.zeros((left_repeats*num_orb_left, size_M-num_orb_left*no_neigh_left))), axis=1);
	M_bottom = np.concatenate((np.zeros((size_M, left_repeats*num_orb_left)), M), axis = 1)
	M = np.concatenate((M_top, M_bottom), axis = 0)
	
	
	# **************************
	# Repeating the right contact
	# **************************
	
	size_M = np.shape(M)[0]
	right_repeats = block_repeat[1]
	right_contact_block = M[size_M-num_orb_right:size_M, size_M-num_orb_right:size_M]
			
	Mnew = np.concatenate((np.zeros((num_orb_right*no_neigh_right, num_orb_right*right_repeats)),\
						   np.kron(np.identity(right_repeats), right_contact_block)), axis=0)
	
	Medge = M[size_M-num_orb_right*(no_neigh_right+1):size_M-num_orb_right, size_M-num_orb_right:size_M]
	size_Medge = np.shape(Medge)[0]
	
	for ind in range(1, no_neigh_right+1):
		Mnew = Mnew + np.concatenate((np.zeros(((no_neigh_right-ind)*num_orb_right, num_orb_right*right_repeats)),\
									  np.kron(np.identity(right_repeats), Medge[-(ind*num_orb_right):size_Medge-(ind-1)*num_orb_right, :]),\
									  np.zeros((ind*num_orb_right, num_orb_right*right_repeats))), axis=0)
	
	# Add the right contact blocks to the main matrix:
	M_top = np.concatenate((M, np.zeros((right_repeats*num_orb_right, size_M))), axis = 0)
	M_bottom = np.concatenate((np.zeros((size_M-num_orb_right*no_neigh_right, right_repeats*num_orb_right)), Mnew), axis=0)
	M = np.concatenate((M_top, M_bottom), axis = 1)
	
	
	# **************************
	# Preparing the boundary
	# **************************
	
	size_M = np.shape(M)[0]
	no_left = no_neigh_left+1
	no_right = no_neigh_right+1
	
	# Top left block (left boundary)
	Mnew = np.zeros(2*no_left*num_orb_left)
	for ind in range(no_left):
		Mnew = Mnew + np.kron(np.diag(np.ones(2*no_left-ind), ind), M[:num_orb_left, ind*num_orb_left:(ind+1)*num_orb_left])
	
	M[:2*no_left*num_orb_left, :2*no_left*num_orb_left] = Mnew
	
	# Bottom right block (right boundary)
	Mnew = np.zeros(2*no_right*num_orb_right)
	for ind in range(no_right):
		Mnew = Mnew + np.kron(np.diag(np.ones(2*no_right-ind), ind), M[size_M-num_orb_right*(ind+1):size_M-num_orb_right*ind, size_M-num_orb_right:size_M])
		
	M[size_M-num_orb_right*no_right*2:size_M, size_M-num_orb_right*no_right*2:size_M] = Mnew
	
	# Create the full Hamiltonian
	M = M + np.transpose(M) - np.diag(np.diag(M))
	
	#Save plot of sparsity pattern
	plt.spy(M)
	plt.title('Sparsity pattern of the device hamiltonian')
	plt.savefig("H.png")
	
	
	# **************************
	# Updating the atomic coordinates file with the added blocks
	# **************************
	
	# Create extra atomic coordinates corresponding to new left contact blocks
	ext_left_contact_x = np.tile(coords[:no_atoms_per_block[0], 0], (left_repeats, 1)).ravel('C') - np.tile(np.arange(left_repeats,0,-1)*d[0], (no_atoms_per_block[0], 1)).ravel('F')	
	ext_left_contact_yz = np.tile(coords[:no_atoms_per_block[0], 1:], (left_repeats, 1))
	left_atoms = np.concatenate((ext_left_contact_x[:,None], ext_left_contact_yz), axis=1)
	
	# Create extra atomic coordinates corresponding to new right contact blocks
	ext_right_contact_x = np.tile(coords[-no_atoms_per_block[1]:, 0], (right_repeats, 1)).ravel('C') + np.tile(np.arange(1,right_repeats+1,1)*d[1], (no_atoms_per_block[1], 1)).ravel('F')
	ext_right_contact_yz = np.tile(coords[-no_atoms_per_block[1]:, 1:], (right_repeats, 1))
	right_atoms = np.concatenate((ext_right_contact_x[:,None], ext_right_contact_yz), axis=1)
	
	# Stitch the new device coordinates together:
	LM = np.concatenate((left_atoms, coords, right_atoms), axis=0)
	
	# Move the x-positions in the coord array to start at zero
	LM[:,0] = LM[:,0] - np.min(LM[:,0])
	
	return M, LM		

def print_lattice_files(LM, atomic_kinds):
	
	''' 
	Internal fxn, prints a plain text with the following format:
	
	@Manasa: write this docstring

	'''
	
	# Construct the lattice_dat file
	lattice = deepcopy(LM[:, :3])

	# The fourth column of the LM matrix holds the index of the corresponding atom name in atomic_kinds
	atom_names = [atomic_kinds[int(index)-1] for index in LM[:, -1]]
	
	# Print the LM_dat and lattice_dat files:
	with open('LM_dat', 'w') as filehandle:
		for row in LM:
			filehandle.write('{:.7f}\t{:.7f}\t{:.7f}\t{:.7f}\n'.format(row[0], row[1], row[2], row[3]))

	with open('lattice_dat', 'w') as filehandle:
		# @Manasa: Add the cell size and the # atoms at the top
		for atom, row in zip(atom_names, lattice):
			filehandle.write('{}\t{:.7f}\t{:.7f}\t{:.7f}\n'.format(atom, row[0], row[1], row[2]))
			
	
	
def print_Smin_file(LM, no_blocks, no_atoms_first_block):
	''' 
	Internal fxn, prints a plain text with the following format:
	
	@Manasa: write this docstring

	'''
	
	num_atoms = np.shape(LM)[0]
	num_left_contact_atoms = no_blocks[0]*no_atoms_first_block[0]
	num_right_contact_atoms = no_blocks[2]*no_atoms_first_block[1]
	num_channel_atoms = num_atoms - num_left_contact_atoms - num_right_contact_atoms
	channel_atoms_per_block = num_channel_atoms/no_blocks[1]
	
	Smin = [x*no_atoms_first_block[0]+1 for x in range(no_blocks[0]+1)]
	Smin.extend([Smin[-1]+x*channel_atoms_per_block for x in range(1, no_blocks[1]+1)])
	Smin.extend([Smin[-1]+x*no_atoms_first_block[1] for x in range(1, no_blocks[2]+1)])
	Smin[-1] = Smin[-1] -1
	
	with open('Smin_dat', 'w') as filehandle:
		filehandle.write('{:d}\n'.format(len(Smin)-1))
		for index in Smin:
			filehandle.write('{:d}\n'.format(int(index)))
	
def get_warnings(M):
	''' 
	Internal fxn, checks the matrices for common errors, and prints warnings to a file 
	
	@Manasa: write this docstring

	'''

	# @Manasa: write this function

	pass

if __name__ == '__main__':
	
	print('**************** Making inputs for OMEN ********************')
	
	# Input dictionary entries (@Manasa: put this in a json later)
	no_blocks = [5, 4, 5]
	no_atoms_first_block = [128, 128]
	delete_blocks = [0, 0]
	repeat_blocks = [3, 3]
	eps = 1e-6
	xyz_file = 'structure.xyz'
	KS_file = 'geoopt-KS_SPIN_1-1_0.csr'
	S_file = 'geoopt-KS_SPIN_1-1_0.csr'
	
	# Get atomic structure information:
	atomic_kinds, no_orbitals = extract_basis_set(output_log='log_energy.out')
	lattice, atoms, coords = read_xyz(xyz_file)
	num_orb_per_atom = np.array([int(no_orbitals[atomic_kinds.index(atom)]) for atom in atoms])
	coords = np.column_stack((coords, np.array([int(atomic_kinds.index(atom))+1 for atom in atoms])))
	print(f'found {len(atomic_kinds)} atomic kinds: {atomic_kinds}, with corresponding # orbitals: {no_orbitals}')
	
	# Get Kohn-Sham and Overlap matrices from bin files in index-value format
	print('Reading binary files...')
	KS = read_bin(binfile=KS_file, struct_fmt='<IIIdI')
	S = read_bin(binfile=S_file, struct_fmt='<IIIdI')
	
	# Convert to full matrices
	hartree_to_eV = 27.2114
	KS = bin_to_csr(KS)*hartree_to_eV
	S = bin_to_csr(S)
	
	# Create the device matrices by adding/removing contact blocks:
	Hmax = np.amax(np.abs(KS))
	print('Building hamiltonian...')
	H, LM = create_device_matrix(KS, coords, num_orb_per_atom, no_blocks, no_atoms_first_block, no_orbitals, delete_blocks, repeat_blocks, Hmax*eps)
	print('Building overlap matrix...')
	S = create_device_matrix(S, coords, num_orb_per_atom, no_blocks, no_atoms_first_block, no_orbitals, delete_blocks, repeat_blocks, Hmax*eps*0.1) 	

	print('Writing output files...')
	
	# Print the LM_dat and lattice_dat files:
	print_lattice_files(LM, atoms, atomic_kinds, no_atoms_first_block, repeat_blocks)
	
	# Print the Smin file (atomic index of the end of each block)
	print_Smin_file(LM, no_blocks, no_atoms_first_block)
	
	# Checking the matrices for building errors:
	get_warnings(H)
	
	# Write binary files for the hamiltonian and overlap matrices

	print('Finished pre-processing, matrices and input files are ready to use.')	
	
	
	
	
	
	
	
	
	
	
	
	
