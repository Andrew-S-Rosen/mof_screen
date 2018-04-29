import numpy as np
import os
from ase.io import read
from copy import deepcopy
from pymofscreen.metal_types import mag_list, spblock_metals, dblock_metals, fblock_metals, poor_metals
from pymofscreen.writers import pprint

def get_incar_magmoms(incarpath,poscarpath):
#Get the magnetic moments from the POSCAR

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
#Get indices of potentially magnetic atoms

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

	with open(incarpath,'r') as incarfile:
		for line in incarfile:
			line = line.strip()
			if 'ISPIN = 2' in line:
				mof_magmoms = mof.get_magnetic_moments()
				mof.set_initial_magnetic_moments(mof_magmoms)

	return mof

def get_abs_magmoms(mof,incarpath):

	mag_indices = get_mag_indices(mof)
	ispin = False
	with open(incarpath,'r') as incarfile:
		for line in incarfile:
			line = line.strip()
			if 'ISPIN = 2' in line:
				ispin = True
				mof_magmoms = mof.get_magnetic_moments()
				abs_magmoms = np.abs(mof_magmoms[mag_indices])
	if ispin == False:
		abs_magmoms = np.zeros(len(mag_indices))

	return abs_magmoms, mag_indices, ispin

def continue_failed_magmoms(mof):
#If job failed, try to read magmoms from OUTCAR

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

def check_if_new_spin(screener,mof,refcode,acc_level,current_spin):
	basepath = screener.basepath
	spin_levels = screener.spin_levels
	results_partial_path = os.path.join(basepath,'results',refcode,acc_level)
	success_path = os.path.join(results_partial_path,current_spin)
	incarpath = os.path.join(success_path,'INCAR')
	mof = deepcopy(mof)
	mof = continue_magmoms(mof,incarpath)

	for prior_spin in spin_levels:
		if prior_spin == current_spin:
			continue
		old_mof_path = os.path.join(results_partial_path,prior_spin,'OUTCAR')
		old_incar_path = os.path.join(results_partial_path,prior_spin,'INCAR')
		mag_indices = get_mag_indices(mof)
		old_mof = read(old_mof_path)
		old_abs_magmoms, old_mag_indices, old_ispin = get_abs_magmoms(old_mof,old_incar_path)
		mof_mag = mof.get_initial_magnetic_moments()[mag_indices]
		if old_ispin == True:
			old_mof_mag = old_mof.get_magnetic_moments()[mag_indices]
		else:
			old_mof_mag = [0]*len(mag_indices)
		mag_tol = 0.05
		if np.sum(np.abs(mof_mag - old_mof_mag) >= mag_tol) == 0:
			pprint('Skipping rest because '+current_spin+' converged to '+prior_spin)
			return False

	return True

def check_if_skip_low_spin(screener,mof,refcode,spin_level):
	acc_levels = screener.acc_levels
	acc_level = acc_levels[-1]
	basepath = screener.basepath
	success_path = os.path.join(basepath,'results',refcode,acc_level,spin_level)
	incarpath = os.path.join(success_path,'INCAR')
	skip_low_spin = False

	abs_magmoms, mag_indices, ispin = get_abs_magmoms(mof,incarpath)
	mag_nums = mof[mag_indices].get_atomic_numbers()
	if np.sum(abs_magmoms < 0.1) == len(abs_magmoms) or all(num in spblock_metals+poor_metals for num in mag_nums) == True:
		skip_low_spin = True

	return skip_low_spin