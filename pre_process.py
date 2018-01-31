from pymatgen.io.cif import CifParser
from ase.io import read
import os
import numpy as np
mofpath = '/projects/p30148/vasp_jobs/structures/CoRE1-DFT-OMS-v2/'
n_atoms = []
elements = []
counts = []
refcodes = []
dup_list = []
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
		n_atoms.append(len(mof))
		[elements_temp,counts_temp] = np.unique(atom_numbers,return_counts=True)
		elements.append(elements_temp.tolist())
		counts.append(counts_temp.tolist())
os.remove(mofpath+'POSCAR')
unique_n_atoms,counts_temp = np.unique(n_atoms,return_counts=True)
mult_n_atoms = unique_n_atoms[counts_temp > 1].tolist()
unique_elements,counts_temp = np.unique(elements,return_counts=True)
mult_elements = unique_elements[counts_temp > 1].tolist()
unique_counts,counts_temp = np.unique(counts,return_counts=True)
mult_counts = unique_counts[counts_temp > 1].tolist()
for i, refcode in enumerate(refcodes):
	n_atoms_i = n_atoms[i]
	if n_atoms_i in mult_n_atoms:
		elements_i = elements[i]
		if elements_i in mult_elements:
			counts_i = counts[i]
			if counts_i in mult_counts:
				dup_list.append(refcode)
dup_list = [dup for _,dup in sorted(zip(n_atoms,dup_list))]
for dup in dup_list:
	print('WARNING: '+dup+' likely a duplicate ('+str(n_atoms_i)+', '+str(elements_i)+', '+str(counts_i)+')')
