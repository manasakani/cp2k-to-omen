import os
from functools import wraps

def input_file_existence(working_dir, files):
	def decorate(func):
		@wraps(func)
		def wrapper():
			for name in files:
				if os.path.isfile(os.path.join(working_dir,name)):
					continue
				else:
					raise Exception(f"A '{name}' file should be found in the current folder!")
			return func(files)
		return wrapper
	return decorate
