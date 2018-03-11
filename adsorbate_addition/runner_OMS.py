import os
import numpy as np
from ase.io import read
from geom import get_ase_NN_idx, get_best_to_worst_idx
from settings import omsdata_path, coremof_path
from zeo_handler import get_omsex_data, get_CN
from add_adsorbate import write_oms_file
from ads_sites import get_opt_ads_site
from path_prep import prep_paths, get_refcode

prep_paths()
for struct_file in os.listdir(coremof_path):
	if not os.path.isfile(coremof_path+struct_file):
		continue
	refcode = get_refcode(struct_file)
	if os.stat(omsdata_path+refcode+'.omsex').st_size == 0:
		continue
	mof = read(coremof_path+struct_file)
	n_OMS = get_CN(omsdata_path+refcode+'.oms')
	cnums_all, cus_coords_all, cus_idx_all, cus_sym_all, NN_coords_all = get_omsex_data(refcode,n_OMS,mof)
	cluster_sym = []
	for i, cus_idx in enumerate(cus_idx_all):
		if mof[cus_idx].symbol != cus_sym_all[i]:
			raise ValueError('Reading wrong element')
		sum_prior_cnums = sum(cnums_all[0:i])
		NN_coords = NN_coords_all[sum_prior_cnums:sum_prior_cnums+cnums_all[i],:]
		ase_NN_idx = get_ase_NN_idx(mof,NN_coords)
		if len(ase_NN_idx) != cnums_all[i]:
			raise ValueError('Reading wrong indices')
		if len(ase_NN_idx) == 1:
			NN_atnum_temp = [mof[ase_NN_idx].number]
		else:
			NN_atnum_temp = mof[ase_NN_idx].get_atomic_numbers().tolist()
		oms_atnum_temp = [mof[cus_idx].number]
		NN_atnum_temp.sort()
		cluster_sym.append(oms_atnum_temp+NN_atnum_temp)
	if len(cluster_sym) != len(cus_sym_all):
			raise ValueError('Read wrong indices')
	unique_cluster_sym_all = []
	for entry in cluster_sym:
		if entry not in unique_cluster_sym_all:
			unique_cluster_sym_all.append(entry)
	indices_tot = 0
	for unique_cluster_sym in unique_cluster_sym_all:
		omsex_indices = [idx for idx, entry in enumerate(cluster_sym) if entry == unique_cluster_sym]
		ads_sites = np.zeros((len(omsex_indices),3))
		for i, omsex_idx in enumerate(omsex_indices):
			cnum = cnums_all[omsex_idx]
			cus_sym = cus_sym_all[omsex_idx]
			sum_prior_cnums = sum(cnums_all[0:omsex_idx])
			NN_coords = NN_coords_all[sum_prior_cnums:sum_prior_cnums+cnum,:]
			cus_coords = cus_coords_all[omsex_idx,:]
			cus_idx = cus_idx_all[omsex_idx]
			ase_NN_idx = get_ase_NN_idx(mof,NN_coords)
			mic_coords = np.squeeze(mof.get_distances(cus_idx,ase_NN_idx,mic=True,vector=True))
			ads_sites[i,:] = get_opt_ads_site(struct_file,cnum,mic_coords,cus_idx,cus_coords)
		cus_idx_cluster = [cus_idx_all[j] for j in omsex_indices]
		best_to_worst_idx = get_best_to_worst_idx(struct_file,ads_sites,cus_idx_cluster)
		write_oms_file(refcode,struct_file,ads_sites,best_to_worst_idx,unique_cluster_sym)
		indices_tot += len(omsex_indices)
	if indices_tot != len(cus_sym_all):
		raise ValueError('Did not run through all OMS')