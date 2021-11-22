import os
import struct
import numpy as np
from scipy import sparse


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