import os
from settings import mofpath, skip_mofs, phase
from writers import pprint
from ase.io import read
from pymatgen.io.cif import CifParser

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
#Read MOF as ASE atoms object
#if running Phase 1, do Niggli-reduction. else, get CIF from prior phase

	cifpath = mofpath+cif_file
	if phase == 1:
		parser = CifParser(cifpath)
		pm_mof = parser.get_structures(primitive=True)[0]
		pm_mof.to(filename='POSCAR')
		mof = read('POSCAR')
	else:
		mof = read(cifpath)

	return mof