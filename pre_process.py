from pymatgen.io.cif import CifParser
from ase.io import read
import os
import numpy as np
mofpath = '/projects/p30148/vasp_jobs/structures/CoRE1-DFT-OMS-v2/'
refcodes = []
stoichs = []
for filename in os.listdir(mofpath):
	filename = filename.strip()
	if len(filename.split('.cif')) == 2:
		refcode = filename.split('.cif')[0]
		refcodes.append(refcode)
		parser = CifParser(mofpath+filename)
		pm_mof = parser.get_structures(primitive=True)[0]
		pm_mof.to(filename='POSCAR')
		mof = read('POSCAR')
		atom_numbers = mof.get_atomic_numbers().tolist()
		if 6 not in atom_numbers:
			print('No C in MOF: '+refcode)
			continue
		stoichs.append(mof.get_chemical_formula())
unique_id,id_counts = np.unique(stoichs,return_counts=True)
mult_ids = unique_id[id_counts > 1].tolist()
stoichs = np.array(stoichs)
refcodes = np.array(refcodes)
for mult_id in mult_ids:
	position_mats = []
	idx = []
	for i, stoich in enumerate(stoichs):
		if stoich == mult_id:
			parser = CifParser(mofpath+refcodes[i]+'.cif')
			pm_mof = parser.get_structures(primitive=True)[0]
			pm_mof.to(filename='POSCAR')
			mof = read('POSCAR')
			position_mats.append(mof.get_positions())
			idx.append(i)
	idx = np.array(idx)
	position_mats = np.array(position_mats)
	n = np.shape(position_mats)[0]
	for i in range(n):
		j_vec = np.setdiff1d(np.arange(n),i)
		dups = []
		for j in j_vec:
			diff = (np.sum((position_mats[i,:,:]-position_mats[j,:,:])**2)/n)**(0.5)
			if diff < 0.1:
				dups.append(j)
		if dups:
			print(refcodes[idx[i]]+' ('+stoichs[idx[i]]+') is a duplicate of '+str(refcodes[idx[dups]]))
	print('')
os.remove('POSCAR')
