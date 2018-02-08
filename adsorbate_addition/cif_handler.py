import os
from settings import (coremof_path,newmofs_path,error_path)

def get_cif_files():
#read in the MOF CIF files

	cif_files = []
	for filename in os.listdir(coremof_path):
		filename = filename.strip()
		if len(filename.split('.cif')) == 2:
			cif_files.append(filename)
	if not os.path.exists(newmofs_path):
		os.makedirs(newmofs_path)
	if not os.path.exists(error_path):
		os.makedirs(error_path)

	return cif_files