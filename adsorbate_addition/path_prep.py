from settings import newmofs_path, error_path
import os

def prep_paths():
	
	if not os.path.exists(newmofs_path):
		os.makedirs(newmofs_path)
	if not os.path.exists(error_path):
		os.makedirs(error_path)

def get_refcode(struct_file):

	if '.cif' in struct_file:
		refcode = struct_file.split('.cif')[0]
	elif 'CAR_' in struct_file:
		refcode = struct_file.split('CAR_')[-1]
	elif '_CAR' in struct_file:
		refcode = struct_file.split('_CAR')[0]
	else:
		raise ValueError('Unknown file naming scheme')
		
	return refcode