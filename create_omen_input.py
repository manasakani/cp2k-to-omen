import os
import time
import json
import numpy as np
import matplotlib.pyplot as plt
from copy import deepcopy
from lib import utils
from lib import input_decorators as check

#Load the User Inputs
with open('user_inputs.json') as filehandle:
		user_inputs = json.load(filehandle)


def extract_basis_set(output_log):
	''' 
	Internal. Reads the atomic kinds and # orbitals per atom from a cp2k output log.

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
	
	
def create_device_matrix(M, coords, num_orb_per_atom, blocks, no_atoms_per_block, no_orb, block_delete, block_repeat, eps):
	
	''' 
	Internal. Takes the KS and S matrices and modifies them accordingly to create the device matrices,
	by adding/removing contact blocks, formatting the channel elements, etc.
	
	@Manasa: clean up this function it looks atrocious.

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
	plt.spy(M, markersize=0.01)
	plt.title(f'H Sparsity, nnz = {np.count_nonzero(M)}')
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
	Internal. Prints the lattice_dat and LM_dat files for OMEN.
	
	'''
	
	# Construct the lattice_dat file
	lattice = deepcopy(LM[:, :3])

	# The fourth column of the LM matrix holds the index of the corresponding atom name in atomic_kinds
	atom_names = [atomic_kinds[int(index)-1] for index in LM[:, -1]]
	
	# Make 'box' with 'f' % empty space around it:
	f = 1.1
	Lx = (np.max(LM[:,0])-np.min(LM[:,0]))*f
	Ly = (np.max(LM[:,1])-np.min(LM[:,1]))*f
	Lz = (np.max(LM[:,2])-np.min(LM[:,2]))*f
	box = np.diag([Lx, Ly, Lz])
	
	# Don't count vacancies as atomic kinds:
	if 'V' in atomic_kinds:
		num_types = np.shape(atomic_kinds)[0] - 1
	else:
		num_types = np.shape(atomic_kinds)[0]
	
	# Print the LM_dat and lattice_dat files:
	with open('LM_dat', 'w') as filehandle:
		for row in LM:
			filehandle.write('{:.7f}\t{:.7f}\t{:.7f}\t{:.7f}\n'.format(row[0], row[1], row[2], row[3]))

	with open('lattice_dat', 'w') as filehandle:
		# Write cell size:
		filehandle.write('{} {} {} {} {}\n\n'.format(np.shape(LM)[0], num_types, 0, 0, 0))
		# Write device size (periodic repeating unit)
		filehandle.write('{} {} {}\n'.format(box[0, 0], box[0, 1], box[0, 2]))
		filehandle.write('{} {} {}\n'.format(box[1, 0], box[1, 1], box[1, 2]))
		filehandle.write('{} {} {}\n\n'.format(box[2, 0], box[2, 1], box[2, 2]))
		
		for atom, row in zip(atom_names, lattice):
			filehandle.write('{}\t{:.7f}\t{:.7f}\t{:.7f}\n'.format(atom, row[0], row[1], row[2]))
		
	
	
def print_Smin_file(LM, no_blocks, no_atoms_first_block):
	''' 
	Internal. Prints the Smin file (indices of the start of each block)
	
	'''
	#@Manasa: why are we removing two channel blocks???
	no_blocks[1] = no_blocks[1] -2
	
	num_atoms = np.shape(LM)[0]
	num_left_contact_atoms = no_blocks[0]*no_atoms_first_block[0]
	num_right_contact_atoms = no_blocks[2]*no_atoms_first_block[1]
	num_channel_atoms = num_atoms - num_left_contact_atoms - num_right_contact_atoms
	channel_atoms_per_block = np.ceil(num_channel_atoms/no_blocks[1])
	
	start = 1
	stop = no_atoms_first_block[0]*no_blocks[0]+1
	step = no_atoms_first_block[0]
	Smin = np.arange(start, stop+step, step)
	
	start = num_left_contact_atoms+1 + channel_atoms_per_block
	stop = num_atoms-num_right_contact_atoms
	step = channel_atoms_per_block
	Smin = np.append(Smin, np.arange(start, stop, step))
	
	start = num_atoms - num_right_contact_atoms+1
	stop =  num_atoms
	step = no_atoms_first_block[1]
	Smin = np.append(Smin, np.arange(start, stop, step))
	
	Smin = np.append(Smin, num_atoms) 
	Smin = Smin.astype(int)
	
	with open('Smin_dat', 'w') as filehandle:
		filehandle.write('{:d}\n'.format(len(Smin)-1))
		for index in Smin:
			filehandle.write('{:d}\n'.format(int(index)))
			
	return Smin
	
def print_E_file(Ef, dE_outer, dE_inner, rE_outer, rE_inner):
	''' 
	Internal. Prints the E_dat file (energy)
	
	''' 
	E_outer_low = np.arange(Ef-rE_outer, Ef-rE_inner+dE_outer, dE_outer)
	E_inner = np.arange(Ef-rE_inner+dE_inner, Ef+rE_inner-dE_inner+dE_inner, dE_inner)
	E_outer_high = np.arange(Ef+rE_inner, Ef+rE_outer, dE_outer)
	E = np.concatenate((E_outer_low, E_inner, E_outer_high), axis=0)
	
	with open('E_dat', 'w') as f:
		f.write('{}\n'.format(len(E)))
		for energy_pt in E:
			f.write('{}\n'.format(energy_pt))
	

def print_matpar_file(atomic_kinds, no_orbitals, Ef):
	''' 
	Internal. Prints the mar_par file
		
	'''  
	no_atoms = np.shape(atomic_kinds)[0]
	with open('mat_par', 'w') as f:
		f.write('{}\t{}\n'.format(int(np.ceil(no_atoms/2)), int(np.floor(no_atoms/2))))
		f.write('{}\t{:.7f}\t{:.7f}\n'.format(0.001, Ef+1e-3, Ef))
		for orbital in no_orbitals:
			f.write(str(orbital)+' ')
	f.close()


def get_warnings(M):
	
	''' 
	Internal. Checks the matrices for common errors, and prints warnings to a file 
	
	'''
	# @Manasa: write this function and return the total # warnings to be printed.

	pass
	
def clean_matrix(M, Smin, num_orb_per_atom):
	
	''' 
	Internal. If the matrix has non-zero interaction parameters beyond the expected # of nearest-neighbor blocks,
	this function manual sets them to zero and returns the thus 'cleaned' matrix
	
	Args:
		(*x* numpy array) matrix to be cleaned
		
	Returns:
		(*x* numpy array) cleaned matrix 
	
	'''
	
	# Find the block sizes (# orbitals per block)
	block_size = np.zeros(len(Smin)-1)
	for ind in range(len(block_size)-1):
		block_size[ind] = sum(num_orb_per_atom[Smin[ind]-1:Smin[ind+1]-1])
	
	block_size[-1] = sum(num_orb_per_atom[Smin[-2]-1:Smin[-1]])
	blocks = np.cumsum(block_size)+1
	blocks = np.concatenate(([1], blocks), axis=0)
	blocks = blocks.astype(int)
		
	neigh = -1*np.ones((2, len(blocks)-1))
			
	# Calculate forward neighbors
	for ind1 in range(np.shape(neigh)[1]):
		for ind2 in range(ind1, np.shape(neigh)[1]):
			if (M[blocks[ind1]-1:blocks[ind1+1]-2 , blocks[ind2]-1:blocks[ind2+1]-2] > 0).any():
				neigh[0, ind1] += 1
		
	# Calculate backward neighbors
	for ind1 in range(np.shape(neigh)[1]):
		for ind2 in range(ind1, -1, -1):
			if (M[blocks[ind2]-1:blocks[ind2+1]-2 , blocks[ind1]-1:blocks[ind1+1]-2] > 0).any():
				neigh[1, ind1] += 1
		
	# Delete the entries that are beyond the intended neighbors
	nn = int(neigh[0, 0])
	for ii in range(np.shape(neigh)[1] - nn):
		for jj in range(ii+nn+1, np.shape(neigh)[1]):
			if np.count_nonzero(M[blocks[ii]:blocks[ii+1], blocks[jj]:blocks[jj+1]]) > 0:
				M[blocks[ii]:blocks[ii+1], blocks[jj]:blocks[jj+1]] = 0
				M[blocks[jj]:blocks[jj+1], blocks[ii]:blocks[ii+1]] = 0	
						
	return M


input_files = [user_inputs['xyz_file'], user_inputs['KS_file'], user_inputs['S_file'], user_inputs['output_log']]
@check.input_file_existence(os.getcwd(), input_files)
def main(input_files):
		
	'''
	Main process to create hamiltonian, overlap matrix, and other necessary input files for OMEN calculations
	
	'''
	
	print('**************** Making inputs for OMEN ********************')

	t0 = time.time()
	
	# Input dictionary entries (@Manasa: put this in a json later)
	no_blocks = user_inputs['no_blocks']
	no_atoms_first_block = user_inputs['no_atoms_first_block']
	delete_blocks = user_inputs['delete_blocks']
	repeat_blocks = user_inputs['repeat_blocks']
	eps = user_inputs['eps']
	
	# Parameters:
	hartree_to_eV = 27.2114
	dE_outer = 10e-3
	dE_inner = 2e-3
	rE_outer = 5
	rE_inner = 0.1
	Vd = user_inputs['Vd']
	
	# Get atomic structure information:
	atomic_kinds, no_orbitals = extract_basis_set(user_inputs['output_log'])
	lattice, atoms, coords = utils.read_xyz(user_inputs['xyz_file'])	
	num_orb_per_atom = np.array([int(no_orbitals[atomic_kinds.index(atom)]) for atom in atoms])
	coords = np.column_stack((coords, np.array([int(atomic_kinds.index(atom))+1 for atom in atoms])))
	print(f'found {len(atomic_kinds)} atomic kinds: {atomic_kinds}, with corresponding # orbitals: {no_orbitals}')
	
	
	# Get Kohn-Sham and Overlap matrices from bin files in index-value format
	if os.path.isfile(os.getcwd()+'/H.dat') and os.path.isfile(os.getcwd()+'/S.dat'):
		H = np.load('H.dat', allow_pickle=True)
		S = np.load('S.dat', allow_pickle=True)
		t1 = time.time()
		print('Read matrices from pickles in '+str(t1-t0)+' s')
	else:
		H = utils.read_bin(binfile=user_inputs['KS_file'], struct_fmt='<IIIdI')
		S = utils.read_bin(binfile=user_inputs['S_file'], struct_fmt='<IIIdI')
		H.dump('H.dat')
		S.dump('S.dat')
		t1 = time.time()
		print('Read matrices from binary files in '+str(t1-t0)+' s')
		
		
	# Convert to full matrices
	H = utils.bin_to_csr(H)*hartree_to_eV
	S = utils.bin_to_csr(S)
		
	# Create the device matrices by adding/removing contact blocks:
	Hmax = np.amax(np.abs(H))
	H, LM = create_device_matrix(H, coords, num_orb_per_atom, no_blocks, no_atoms_first_block, no_orbitals, delete_blocks, repeat_blocks, Hmax*eps)
	S = create_device_matrix(S, coords, num_orb_per_atom, no_blocks, no_atoms_first_block, no_orbitals, delete_blocks, repeat_blocks, Hmax*eps*0.1)[0]
	t2 = time.time()
	print('Built Hamiltonian and overlap matrices in '+str(t2-t1)+' s')

	# Print text files:
	print_lattice_files(LM, atomic_kinds)	
	Smin = print_Smin_file(LM, no_blocks, no_atoms_first_block)
	Ef = utils.get_value_from_file(user_inputs['output_log'], 'Fermi level')*hartree_to_eV
	print_E_file(Ef, dE_outer, dE_inner, rE_outer, rE_inner)
	print_matpar_file(atomic_kinds, no_orbitals, Ef)
	
	# Checking the matrices for building errors:
	get_warnings(H)

	# Cleaning the matrices
	print('Cleaning matrix entries beyond the expected # of nearest neighbors...')
	num_orb_per_atom = np.array([no_orbitals[int(index)-1] for index in LM[:, -1]])
	H = clean_matrix(H, Smin, num_orb_per_atom)
	S = clean_matrix(S, Smin, num_orb_per_atom)	

	# Write binary files for the hamiltonian and overlap matrices
	print('Writing Hamiltonian and Overlap matrices to .bin...')
	utils.write_mat_to_bin('H_4.bin', H, Hmax*eps)
	utils.write_mat_to_bin('S_4.bin', S, Hmax*eps*0.1)
	t3 = time.time()
	print('Binary files written in '+str(t3-t2)+' s')
	
	# Load the command file dictionary and modify it
	dimensions_min = np.min(LM[:, :3], axis =0)/10
	dimensions_max = np.max(LM[:, :3], axis =0)/10
	with open(os.getcwd()+'/lib/input_templates/omen_cmd.json') as cmd_json_file:
		omen_cmds = json.load(cmd_json_file)
			
	#@Manasa: Ask Fabian what this '-6' is
	omen_cmds.update({
	"fermi_level": Ef,
	"Vdmin": Vd,
	"Vdmax": Vd,
	"restart": "[2 0 0 0]",
	"vact_file": 'vact_dat',
	"Lc": dimensions_max[0]-6+1e-3,
	"tc": dimensions_max[1]+1e-3,
	"hc": dimensions_max[2]+1e-3,
	"x0": dimensions_min[1]-1e-3,
	"z0": dimensions_min[2]-1e-3})
	utils.dump_dict_plaintext('omen.cmd', omen_cmds, 'C')
	
	# Edit job script:
	with open(os.getcwd()+'/lib/input_templates/run_transport_sh.json') as job_json_file:
		jobfile = json.load(job_json_file)
	utils.dump_dict_plaintext('run_transport.sh', jobfile, 'bash')
		
	t4 = time.time()
	print('Finished pre-processing. Matrices and input files are ready to use.')
	print('Total runtime '+str(t4-t0)+' s')	

if __name__ == '__main__':
	main()
	
	
	
	
	
	
	
	
	
	
	
	
	
