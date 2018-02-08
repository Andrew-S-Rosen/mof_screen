import numpy as np
from metal_types import mag_list, dblock_metals, fblock_metals, poor_metals
from ase.io import read

def get_incar_magmoms(incarpath,poscarpath):
#get the magnetic moments from the POSCAR
	mof_mag_list = []
	init_mof = read(poscarpath)
	with open(incarpath,'r') as incarfile:
		for line in incarfile:
			line = line.strip()
			if 'MAGMOM' in line:
				mag_line = line.split('= ')[1:][0].split(' ')
				for val in mag_line:
					mag = float(val.split('*')[1])
					num = int(val.split('*')[0])
					mof_mag_list.extend([mag]*num)
	if bool(mof_mag_list) == False:
		mof_mag_list = np.zeros(len(init_mof))
	if len(mof_mag_list) != len(mof_mag_list):
		raise ValueError('Error reading INCAR magnetic moments')
	return mof_mag_list

def get_mag_indices(mof):
#Get indices of d-block, f-block, and semimetal atoms
	mag_indices = []
	for i, atom in enumerate(mof):
		if atom.number in mag_list:
			mag_indices.append(i)
	return mag_indices

def set_initial_magmoms(mof,spin_level):
#Add initial magnetic moments to atoms object
	mag_indices = get_mag_indices(mof)
	mof.set_initial_magnetic_moments(np.zeros(len(mof)))
	for i, atom in enumerate(mof):
		if i in mag_indices:
			mag_number = atom.number
			if spin_level == 'spin1':
				if mag_number in dblock_metals:
					atom.magmom = 5.0
				elif mag_number in fblock_metals:
					atom.magmom = 7.0
				elif mag_number in poor_metals:
					atom.magmom = 0.1
				else:
					raise ValueError('Metal not properly classified')
			elif spin_level == 'spin2':
				atom.magmom = 0.1
			else:
				raise ValueError('Spin iteration out of range')
	return mof

def continue_magmoms(mof,incarpath):
#Read in the old magmoms
	mag_indices = get_mag_indices(mof)
	ispin = False
	with open(incarpath,'r') as incarfile:
		for line in incarfile:
			line = line.strip()
			if 'ISPIN = 2' in line:
				ispin = True
				mof_magmoms = mof.get_magnetic_moments()
				mof.set_initial_magnetic_moments(mof_magmoms)
				abs_magmoms = np.abs(mof_magmoms[mag_indices])
	if ispin == False:
		abs_magmoms = np.zeros(len(mag_indices))
	return mof, abs_magmoms

def continue_failed_magmoms(mof):
	self_resort = []
	file = open('ase-sort.dat', 'r')
	lines = file.readlines()
	file.close()
	for line in lines:
	    data = line.split()
	    self_resort.append(int(data[1]))
	magnetic_moments = np.zeros(len(mof))
	n = 0
	lines = open('OUTCAR', 'r').readlines()
	for line in lines:
	    if line.rfind('magnetization (x)') > -1:
	        for m in range(len(mof)):
	            magnetic_moments[m] = float(lines[n + m + 4].split()[4])
	    n += 1
	sorted_magmoms = np.array(magnetic_moments)[self_resort]
	ispin = False
	with open('INCAR','r') as incarfile:
		for line in incarfile:
			line = line.strip()
			if 'ISPIN = 2' in line:
				ispin = True
	if ispin == True and all(sorted_magmoms == 0.0) == True:
		raise ValueError('Error reading magmoms from failed OUTCAR')
	mof.set_initial_magnetic_moments(sorted_magmoms)
	return mof