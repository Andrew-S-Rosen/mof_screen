from cif_handler import get_cif_files
import os
from ase.io import read
from geom import get_ase_NN_idx, get_dist_planar, get_best_to_worst_idx
from settings import ads_species, omsdata_path, coremof_path, guess_length, sum_cutoff, rmse_tol
from zeo_handler import get_omsex_data, get_CN
import numpy as np
from regression import OLS_fit, TLS_fit
from add_adsorbate import write_files
from ads_sites import get_planar_ads_site, get_nonplanar_ads_site, get_bi_ads_site, get_tri_ads_site

cif_files = get_cif_files()
for cif_file in cif_files:
	refcode = cif_file.split('.cif')[0]
	if os.stat(omsdata_path+refcode+'.omsex').st_size == 0:
		continue
	basename = refcode+'_'+ads_species
	mof = read(coremof_path+cif_file)
	n_OMS = get_CN(omsdata_path+refcode+'.oms')
	cnums_all, cus_coords_all, ase_cus_idx_all, cus_sym_all, NN_coords_all = get_omsex_data(refcode,n_OMS,mof)
	cluster_sym = []
	for i, ase_cus_idx in enumerate(ase_cus_idx_all):
		if mof[ase_cus_idx].symbol != cus_sym_all[i]:
			raise ValueError('Reading wrong element')
		sum_prior_cnums = sum(cnums_all[0:i])
		NN_coords = NN_coords_all[sum_prior_cnums:sum_prior_cnums+cnums_all[i],:]
		ase_NN_idx = get_ase_NN_idx(mof,NN_coords,refcode)
		if len(ase_NN_idx) != cnums_all[i]:
			raise ValueError('Reading wrong indices')
		if len(ase_NN_idx) == 1:
			NN_atnum_temp = [mof[ase_NN_idx].number]
		else:
			NN_atnum_temp = mof[ase_NN_idx].get_atomic_numbers().tolist()
		oms_atnum_temp = [mof[ase_cus_idx].number]
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
			ase_cus_idx = ase_cus_idx_all[omsex_idx]
			ase_NN_idx = get_ase_NN_idx(mof,NN_coords,refcode)
			mic_coords = np.squeeze(mof.get_distances(ase_cus_idx,ase_NN_idx,mic=True,vector=True))
			if cnum == 1:
				normal_vec = mic_coords
			elif cnum == 2:
				normal_vec = OLS_fit(mic_coords)			
			elif cnum >= 3:
				scaled_mic_coords = mic_coords*guess_length/np.linalg.norm(mic_coords,axis=1)[np.newaxis].T
				scaled_sum_dist = sum(scaled_mic_coords)
				sum_dist = sum(mic_coords)
				norm_scaled = np.linalg.norm(scaled_sum_dist)
				rmse, normal_vec = TLS_fit(mic_coords)
			if cnum == 1:
				dist = get_dist_planar(normal_vec)
				ads_sites[i,:] = O_coords-dist
			elif cnum == 2:
				ads_sites[i,:] = get_bi_ads_site(cif_file,normal_vec,cus_coords,ase_cus_idx)
			elif cnum == 3 and np.linalg.norm(scaled_sum_dist) > sum_cutoff:
				ads_sites[i,:] = get_tri_ads_site(cif_file,normal_vec,sum_dist,cus_coords,ase_cus_idx)
			elif norm_scaled <= sum_cutoff or rmse <= rmse_tol:
				dist = get_dist_planar(normal_vec)
				ads_sites[i,:] = get_planar_ads_site(cif_file,cus_coords,dist,ase_cus_idx)
			else:
				ads_sites[i,:] = get_nonplanar_ads_site(sum_dist,cus_coords)
		ase_cus_idx_cluster = [ase_cus_idx_all[j] for j in omsex_indices]
		best_to_worst_idx = get_best_to_worst_idx(cif_file,ads_sites,ase_cus_idx_cluster)
		write_files(refcode,cif_file,ads_sites,best_to_worst_idx,unique_cluster_sym)
		indices_tot += len(omsex_indices)
	if indices_tot != len(cus_sym_all):
		raise ValueError('Did not run through all OMS')
