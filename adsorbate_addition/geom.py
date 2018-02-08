import numpy as np
from ase.io import read
from ase import Atoms, Atom
from settings import coremof_path, ads_species, rcut, guess_length

def get_ase_NN_idx(mof,coords,refcode):
#get ASE indices for coordinating atoms
	ase_NN_idx = []
	zeo_tol = 0.1
	for i in range(np.shape(coords)[0]):
		nn_fail = False
		for j, element in enumerate(mof):
			if sum(coords[i,:] >= element.position-zeo_tol) == 3 and sum(coords[i,:] <= element.position+zeo_tol) == 3:
				ase_NN_idx.append(j)
				nn_fail = True
				break
		if nn_fail == False:
			print('WARNING with '+refcode+': a zeo++ NN (#'+str(i)+') is not in same spot as in ASE CIF')
	return ase_NN_idx

def get_dist_planar(normal_vec):
#get distance vector of adsorbate from planar OMS
	unit_normal = normal_vec/np.linalg.norm(normal_vec)
	dist = unit_normal*guess_length
	return dist

def get_NNs(cif_file,ads_site,ase_cus_idx):
#get the number of NNs and the minimum distance to a NN (other than the OMS)
	mof_temp = read(coremof_path+cif_file)
	adsorbate = Atoms([Atom(ads_species,ads_site)])
	mof_temp.extend(adsorbate)
	compare_with = np.arange(0,len(mof_temp)-1).tolist()
	del compare_with[ase_cus_idx]
	neighbor_dist = mof_temp.get_distances(len(mof_temp)-1,compare_with,mic=True)
	NN = sum(neighbor_dist <= rcut)
	mindist = np.min(neighbor_dist)
	return NN, mindist

def get_best_to_worst_idx(cif_file,ads_sites,ase_cus_idx_list):
#sort the OMS by smallest NNs (using mindist as tiebreaker)
	NN = []
	mindist = []
	i_vec = []
	best_to_worst_idx = []
	if len(ase_cus_idx_list) != np.shape(ads_sites)[0]:
		raise ValueError('Incompatible lengths of lists')
	for i, ase_cus_idx in enumerate(ase_cus_idx_list):
		NN_temp, mindist_temp = get_NNs(cif_file,ads_sites[i,:],ase_cus_idx)
		NN.append(NN_temp)
		mindist.append(mindist_temp)
		i_vec.append(i)
	merged_list = list(zip(i_vec,NN,mindist))
	merged_list.sort(key=lambda x: x[2],reverse=True)
	merged_list.sort(key=lambda x: x[1])
	for item in merged_list:
		best_to_worst_idx.append(item[0])
	return best_to_worst_idx