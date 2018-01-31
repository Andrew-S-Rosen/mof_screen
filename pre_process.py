from pymatgen.io.cif import CifParser
from ase.io import read
import os
import numpy as np
mofpath = 'C:/Users/asros/Desktop/test/'
n_atoms = []
elements = []
counts = []
refcodes = []
for filename in os.listdir(mofpath):
	filename = filename.strip()
	if len(filename.split('.cif')) == 2:
		refcode = filename.split('.cif')[0]
		refcodes.append(refcode)
		parser = CifParser(mofpath+filename)
		pm_mof = parser.get_structures(primitive=True)[0]
		pm_mof.to(filename='POSCAR')
		mof = read('POSCAR')
		if 6 not in mof.get_atomic_numbers():
			print('NOT MOF: '+refcode)
			continue
		n_atoms.append(len(mof))
		[elements_temp,counts_temp] = np.unique(mof.get_atomic_numbers(),return_counts=True)
		elements.append(elements_temp)
		counts.append(counts_temp)
unique_n_atoms,counts_temp = np.unique(n_atoms,return_counts=True)
mult_n_atoms = unique_n_atoms[counts_temp > 1]
overlap_n_atoms = set(n_atoms) & set(mult_n_atoms)
unique_elements,counts_temp = np.unique(elements,return_counts=True)
mult_elements = unique_elements[counts_temp > 1]
unique_counts,counts_temp = np.unique(counts,return_counts=True)
mult_counts = unique_counts[counts_temp > 1]
for i, refcode in enumerate(refcodes):
	n_atoms_i = n_atoms[i]
	if n_atoms_i in mult_n_atoms:
		elements_i = elements[i]
		if elements_i in mult_elements:
			counts_i = counts[i]
			if counts_i in mult_counts:
				print('WARNING: '+refcode+' likely a duplicate ('+str(n_atoms_i)+', '+str(elements_i)+', '+str(counts_i)+')')
