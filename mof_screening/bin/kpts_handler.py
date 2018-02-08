import numpy as np
from settings import ads_species, kpts_path
from calculators import defaults

def get_kpts(cif_file,kppa):
#Get kpoint grid at a given KPPA

	cif_split1 = cif_file.split('_'+ads_species+'_OMS')[0]
	old_cif_name = cif_split1.split('_spin')[0]
	infile = open(kpts_path,'r')
	lines = infile.read().splitlines()
	infile.close()
	for i in range(len(lines)):
		if old_cif_name in lines[i]:
			if kppa == defaults.kppa_lo:
				kpts = lines[i+1]
				gamma = lines[i+2]
			elif kppa == defaults.kppa_hi:
				kpts = lines[i+3]
				gamma = lines[i+4]
			else:
				raise ValueError('Incompatible KPPA with prior runs')
			break
	kpts = np.squeeze(np.asarray(np.matrix(kpts))).tolist()
	if gamma == 'True':
		gamma = True
	elif gamma == 'False':
		gamma = False
	else:
		raise ValueError('Error parsing KPOINTS file')
	if len(kpts) != 3:
		raise ValueError('Error parsing KPOINTS file')

	return kpts, gamma

def get_gpt_version(kpts,n_atoms,nprocs,ppn):
#Determine if gamma-point VASP should be used
#Also reduce processors if too many requested

	if sum(kpts) == 3:
		gpt_version = True
	else:
		gpt_version = False
	while n_atoms < nprocs/2:
		nprocs = nprocs - ppn
		if nprocs == ppn:
			break

	return gpt_version, nprocs