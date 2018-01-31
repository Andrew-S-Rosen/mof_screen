from pymatgen.io.cif import CifParser
from ase.io import read
import os
import numpy as np
mofpath = '/projects/p30148/vasp_jobs/structures/CoRE1-DFT-OMS-v2/'
n_atoms = []
elements = []
counts = []
refcodes = []
id_list = []
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
		[elements_temp,counts_temp] = np.unique(atom_numbers,return_counts=True)
		n_atoms_temp = len(mof)
		elements_temp = elements_temp.tolist()
		counts_temp = counts_temp.tolist()
		elements.append(elements_temp)
		counts.append(counts_temp)
		n_atoms.append(n_atoms_temp)
		id_list.append([n_atoms_temp]+elements_temp+counts_temp)
os.remove(mofpath+'POSCAR')
full_id_list = list(zip(refcodes,id_list))
full_id_list.sort(key=lambda x: x[1])
unique_id,id_counts = np.unique(id_list,return_counts=True)
mult_id = unique_id[id_counts > 1].tolist()
for i, refcode in enumerate(refcodes):
	n_atoms_i = n_atoms[i]
	elements_i = elements[i]
	counts_i = counts[i]
	id_i = [n_atoms_i]+elements_i+counts_i
	if id_i in mult_id:
		print('WARNING: '+refcode+' likely a duplicate ('+str(n_atoms_i)+', '+str(elements_i)+', '+str(counts_i)+')')
