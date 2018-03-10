import os
import numpy as np
from pymatgen.io.cif import CifParser
from pymatgen.analysis.structure_matcher import StructureMatcher

mofpath = '/projects/p30148/vasp_jobs/structures/CoRE1-DFT-OMS-v2/'
refcodes = []
stoichs = []
for filename in os.listdir(mofpath):
	filename = filename.strip()
	if len(filename.split('.cif')) == 2:
		refcode = filename.split('.cif')[0]
		parser = CifParser(mofpath+filename)
		mof = parser.get_structures(primitive=True)[0]
		if 'Element C,' not in str(mof.species):
			print('No C in MOF: '+refcode)
			continue
		refcodes.append(refcode)
		stoichs.append(mof.formula.replace(' ',''))
unique_id,id_counts = np.unique(stoichs,return_counts=True)
mult_ids = unique_id[id_counts > 1].tolist()
stoichs = np.array(stoichs)
refcodes = np.array(refcodes)
sorted_idx = np.argsort(refcodes)
refcodes = refcodes[sorted_idx]
stoichs = stoichs[sorted_idx]
all_dup_refcodes = []
for mult_id in mult_ids:
	position_mats = []
	idx = []
	for i, stoich in enumerate(stoichs):
		if stoich == mult_id:
			idx.append(i)
	idx = np.array(idx)
	n = len(idx)
	for i in range(n):
		refcode = refcodes[idx[i]]
		if refcode in all_dup_refcodes:
			continue
		parser1 = CifParser(mofpath+refcodes[idx[i]]+'.cif')
		mof1 = parser1.get_structures(primitive=True)[0]
		j_vec = np.setdiff1d(np.arange(n),i)
		dups = []
		for j in j_vec:
			parser2 = CifParser(mofpath+refcodes[idx[j]]+'.cif')
			sm = StructureMatcher(primitive_cell=True)
			mof2 = parser2.get_structures(primitive=True)[0]
			rms = sm.get_rms_dist(mof1,mof2)
			if rms and rms[0] < 0.1:
				dups.append(j)
		dup_refcodes = refcodes[idx[dups]]
		if dups:
			print(refcode+' ('+stoichs[idx[i]]+') is a duplicate of '+str(dup_refcodes))
		all_dup_refcodes.extend(dup_refcodes)
