# 15.11.21 - Manasa Kaniselvan
import os
import sys
import shutil
import paramiko

IP_dict = {
	'ela_IP': '148.187.1.17',
	'daint_IP': '148.187.26.68'
	}

def connect_to_daint(username, password):
	
	''' 
	Returns a paramiko client connected to daint, upon which an sftp connection can be opened
	
	Inputs: 
		1. Username of account
		2. (Optional) Password of account
		
	Outputs:
		1. Paramiko client called 'daintclient'. SSH with it using "sftp = daintclient.open_sftp()", and then use daintclient.exec_command('')
	'''

	# Connect to ela.cscs.ch (Login)
	elaclient = paramiko.SSHClient()
	elaclient.set_missing_host_key_policy(paramiko.AutoAddPolicy())
	elaclient.connect('ela.cscs.ch', username=username, password=password)
	print('connected to ela.cscs.ch')
	
	# Open a channel between ela and daint
	elaclienttransport = elaclient.get_transport()
	dest_addr = (IP_dict.get('daint_IP'), 22)
	local_addr = (IP_dict.get('ela_IP'), 22)
	client_channel = elaclienttransport.open_channel("direct-tcpip", dest_addr, local_addr)
	
	# Connect to daint.cscs.ch (Compute)
	daintclient = paramiko.SSHClient()
	daintclient.set_missing_host_key_policy(paramiko.AutoAddPolicy())
	daintclient.connect('daint.cscs.ch', username=username, password=password, sock=client_channel)
	print('connected to daint.cscs.ch')
	
	return daintclient

def push_folder_to_remote(source_folder, remote_location, outbound_files, client):
	
	''' 
	Uploads a specified folder to the remote location specified by the client
	
	Inputs: 
		1. Path of local folder which should be uploaded
		2. Location on the server that the folder should be placed in
		3. A paramiko client for the server, upon which sftp commands can be used.
		
	Outputs:
		(none)
	'''
	
	sftp = client.open_sftp()
	try:
		sftp.chdir(remote_location)
	except IOError:
		sftp.mkdir(remote_location)
	for file in outbound_files :
		filepath = source_folder+'/'+file
		moveto = remote_location +'/'+file
		sftp.put(filepath, moveto)
	sftp.close()

def pull_folder_from_remote(remote_location, dest_folder, inbound_files, client):
	
	''' 
	Downloads specified files from a remote server onto the local machine
	
	Inputs: 
		1. Path on remote machine which the files can be found in
		2. Path of destination on local machine that the files should be placed in
		3. List of file names which should be downloaded
		4. A paramiko client for the server, upon which sftp commands can be used.
		
	Outputs:
		(none)
	'''
	
	sftp = client.open_sftp()
	for file in inbound_files:
		filepath = remote_location+'/'+file
		moveto = dest_folder +'/'+file
		sftp.get(filepath, moveto)
	sftp.close()
		
def remote_command(client, cmd):
	
	''' 
	Executes a terminal command on a remote server linked to with 'client'. 
	
	Inputs: 
		1. A paramiko client for the server, upon which sftp commands can be used.
		2. A string of terminal commands to be executed. Multiple commands should be linked by ';' in the same function call
		
	Outputs:
		1. The terminal output upon during execution of the command
		2. Any error messages encountered during execution of the command
	'''
	
	stdin, stdout, stderr  = client.exec_command(cmd)
	return stdout, stderr
	
