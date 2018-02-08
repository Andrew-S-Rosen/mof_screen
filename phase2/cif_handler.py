import os
from settings import mofpath, skip_mofs
from writers import pprint
from ase.io import read

def get_cif_files():
#Get CIF files from mofpath
	cif_files = []
	for filename in os.listdir(mofpath):
		filename = filename.strip()
		if len(filename.split('.cif')) == 2:
			refcode = filename.split('.cif')[0]
			if refcode not in skip_mofs:
				cif_files.append(filename)
			else:
				pprint('Skipping '+refcode)
	sorted_cifs = sorted(cif_files)
	return sorted_cifs

def cif_to_mof(cif_file):
#Save MOF file as reduced unit cell and read in ASE
	cifpath = mofpath+cif_file
	mof = read(cifpath)
	return mof