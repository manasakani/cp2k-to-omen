import os
import sys
import shutil
import paramiko
import numpy as np

sys.path.insert(0, '/home/mkaniselvan/Documents/ScriptLibrary/omen_preprocessing')
from lib import remote_connection as rc
from lib import utils

def main(folder, base_structure_location=os.getcwd()):
	
	'''
	# This script automates the process of preparing and submitting an scf calculation, by
	# 1. Converting the 'vacancy+oxygen' file from KMC to a structure.xyz file
	# 2. Preparing cp2k input for an scf calculation
	# 3. Uploading this with a job script to daint and running the calculation 
	
	Args:
		1. The folder in which either the vacancies.xyz or structure.xyz file is located
		2. If starting from a vacancy file, the folder in which the base structure (incl contacts) is located)
	'''
	
	sim_path = os.getcwd()+folder
	
	# if the structure file isn't there, then the vacancy file should be, and structure.xyz is created from it
	if os.path.exists(sim_path+'/vacancies.xyz') and not os.path.exists(sim_path+'/structure.xyz'):
	
		print('Creating structure.xyz file from the vacancy coordinates...')
		# Make the structure file using the vacancy file
		with open(sim_path+'vacancies.xyz', 'r') as f:
			vacancies = f.readlines()
			
		with open(base_structure_location+'/base_structure.xyz', 'r') as f:
			base_structure = f.readlines()
		
		vacancies = vacancies[2:]
		num_atoms = base_structure[0]
		cell_size = base_structure[1]
		contacts = [entry for entry in base_structure if 'Ti' in entry or 'N' in entry]
		Hf_atoms = [entry for entry in base_structure if 'Hf' in entry]
		device = vacancies + contacts + Hf_atoms
		
		# sort the coordinates:
		x_positions = [float(atom.split()[1]) for atom in device]
		y_positions = [float(atom.split()[2]) for atom in device]
		z_positions = [float(atom.split()[3]) for atom in device] 
		sorted_inds = np.lexsort((z_positions, y_positions, x_positions))
		device = [device[i] for i in sorted_inds]
		
		# Write the structure.xyz file
		f = open(sim_path+'structure.xyz', 'w')
		f.writelines(num_atoms)
		f.writelines(cell_size)
		f.writelines(device)
		f.close()
	
	
	# Get the energy.inp and job submission scripts from the lib folder
	shutil.copyfile(os.getcwd()+'/lib/input_templates/energy.inp', sim_path+'/energy.inp')
	shutil.copyfile(os.getcwd()+'/lib/input_templates/run_energy.sh', sim_path+'/run_energy.sh')
	
	# Make sure the energy.inp scf input is consistent with the structure.xyz file provided
	f = open(sim_path+'structure.xyz', 'r')
	structure_lines = f.readlines()
	cell_size = structure_lines[1][5:] # cell size should be on the second line
	
	modify_energy = {
	'ABC' : '\tABC ' + cell_size,
	}
	utils.replace_line(sim_path+'energy.inp', modify_energy)
	
	# Edit the job name
	modify_job = {
	'--job-name' : '#SBATCH --job-name=HfO2-'+folder+'\n',
	}
	utils.replace_line(sim_path+'run_energy.sh', modify_job)
		
	# Setup sftp connection
	username='mkanisel'
	password=''
	daintclient = rc.connect_to_daint(username, password)
	
	source_folder = sim_path	
	remote_location = '/scratch/snx3000/mkanisel/Nov23_Sweep/' + folder
	outbound_files=os.listdir(source_folder)
	
	print('Transporting the current directory to Daint at '+ remote_location)
	rc.push_folder_to_remote(source_folder, remote_location, outbound_files, client=daintclient)
	
	print('Submitting the job...')
	stdout, stderr = rc.remote_command(client=daintclient, cmd='cd '+ remote_location +'; sbatch run_energy.sh')
	
	print('output:')
	print(stdout.read())
	print('any errors?')
	print(stderr.read())
	daintclient.close()		

if __name__ == '__main__':
	folder = '/6/'
	main(folder)
