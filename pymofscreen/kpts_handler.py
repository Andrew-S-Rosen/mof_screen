import numpy as np
import os
from pymatgen.io.cif import CifParser
from pymatgen.io.vasp.inputs import Kpoints

def get_kpts(screener,cif_file,level):
	"""
	Obtain the number of kpoints
	Args:
		screener (class): pymofscreen.screener class
		cif_file (string): name of CIF file
		level (string): accuracy level
	Returns:
		kpts (list of ints): kpoint grid
		gamma (bool): True for gamma-centered
	"""

	niggli = screener.niggli
	mofpath = screener.mofpath
	kpts_path = screener.kpts_path
	ads_species = screener.ads_species
	kppas = screener.kppas
	kpts = None

	if kpts_path == 'Auto':

		if level == 'low':
			kppa = kppas[0]
		elif level == 'high':
			kppa = kppas[1]
		else:
			raise ValueError('kpoints accuracy level not defined')

		parser = CifParser(os.path.join(mofpath,cif_file))
		pm_mof = parser.get_structures(primitive=niggli)[0]
		pm_kpts = Kpoints.automatic_density(pm_mof,kppa)
		kpts = pm_kpts.kpts[0]

		if pm_kpts.style.name == 'Gamma':
			gamma = True
		else:
			gamma = None

	else:

		cif_split1 = cif_file.split('_'+ads_species+'_OMS')[0]
		old_cif_name = cif_split1.split('_spin')[0]
		infile = open(kpts_path,'r')
		lines = infile.read().splitlines()
		infile.close()

		for i in range(len(lines)):
			if old_cif_name in lines[i]:
				if level == 'low':
					kpts = lines[i+1]
					gamma = lines[i+2]
				elif level == 'high':
					kpts = lines[i+3]
					gamma = lines[i+4]
				else:
					raise ValueError('Incompatible KPPA with prior runs')
				break
		kpts = np.squeeze(np.asarray(np.matrix(kpts))).tolist()
		if len(kpts) != 3:
			raise ValueError('Error parsing KPOINTS file')

		if gamma == 'True':
			gamma = True
		elif gamma == 'False':
			gamma = False
		else:
			raise ValueError('Error parsing KPOINTS file')

	return kpts, gamma