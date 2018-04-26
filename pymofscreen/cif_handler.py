import os
from ase.io import read
import pymatgen as pm
from pymatgen.io.cif import CifParser
from pymofscreen.writers import pprint

def get_cif_files(mofpath,skip_mofs=None):
	"""
	Get the list of CIF files
	Args:
		mofpath (string): directory to CIF files 
		skip_mofs (list): list of MOFs to ignore
	Returns:
		sorted_cifs (list): alphabetized list of CIF files
	"""
	cif_files = []
	if skip_mofs is None:
		skip_mofs = []
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

def cif_to_mof(filepath,niggli):
	"""
	Convert file to ASE Atoms object
	Args:
		filepath (string): full path to structure file
		niggli (bool): if Niggli-reduction should be performed
	Returns:
		sorted_cifs (list): alphabetized list of CIF files
	"""
	if niggli == True:
		if '.cif' in os.path.basename(filepath):
			parser = CifParser(filepath)
			pm_mof = parser.get_structures(primitive=True)[0]
		else:
			pm.Structure.from_file(filepath,primitive=True)
		pm_mof.to(filename='POSCAR')
		mof = read('POSCAR')
	else:
		mof = read(filepath)

	return mof
