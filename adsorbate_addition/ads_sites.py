import numpy as np
from geom import get_NNs, get_dist_planar
from regression import OLS_fit, TLS_fit
from settings import guess_length, overlap_tol, rcut, coremof_path, ads_species, sum_cutoff, rmse_tol
from ase import Atoms, Atom
from ase.io import read, write
import os

def get_planar_ads_site(struct_file,center_coord,dist,center_idx):
#Get adsorption site for planar structure

	NN = []
	mindist = []
	for i in range(2):
		if i == 0:
			ads_site_temp = center_coord + dist
		elif i == 1:
			ads_site_temp = center_coord - dist
		NN_temp, mindist_temp = get_NNs(struct_file,ads_site_temp,center_idx)
		NN.append(NN_temp)
		mindist.append(mindist_temp)
	if NN[0] == NN[1]:
		if mindist[0] >= mindist[1]:
			ads_site = center_coord + dist
		else:
			ads_site = center_coord - dist
	elif NN[0] <= NN[1]:
		ads_site = center_coord + dist
	else:
		ads_site = center_coord - dist

	return ads_site

def get_nonplanar_ads_site(sum_dist,center_coord):
#Get adsorption site for nonplanar structure

	dist = guess_length*sum_dist/np.linalg.norm(sum_dist)
	ads_site =  center_coord - dist

	return ads_site

def get_bi_ads_site(struct_file,normal_vec,center_coord,center_idx):
#Get adsorption site for 2-coordinate geometry

	try_angles = np.arange(0,360,10)
	dist = get_dist_planar(normal_vec)
	ads_site_temp_unrotated1 = center_coord + dist
	ads_site_temp_unrotated2 = center_coord - dist
	ads_temp1 = Atoms([Atom(ads_species,ads_site_temp_unrotated1)])
	ads_temp2 = Atoms([Atom(ads_species,ads_site_temp_unrotated2)])
	mof_temp = read(coremof_path+struct_file)
	write('mof_temp.cif',mof_temp)
	for i, angle in enumerate(try_angles):
		mof_temp = read('mof_temp.cif')
		mof_temp.extend(ads_temp1)
		mof_temp.extend(ads_temp2)
		mof_temp.set_distance(center_idx,len(mof_temp)-1,guess_length,fix=0,mic=True)
		mof_temp.set_distance(center_idx,len(mof_temp)-2,guess_length,fix=0,mic=True)	
		mof_temp.set_angle(len(mof_temp)-1,center_idx,len(mof_temp)-2,angle)
		dist_mat = mof_temp.get_distances(len(mof_temp)-2,np.arange(0,len(mof_temp)-2).tolist(),mic=True)
		NNs = sum(dist_mat <= rcut)
		if i == 0:
			ads_site = mof_temp[-2].position
			old_min_NNs = NNs
		elif sum(dist_mat <= overlap_tol) == 0 and NNs < old_min_NNs:
			ads_site = mof_temp[-2].position
			old_min_NNs = NNs
	os.remove('mof_temp.cif')

	return ads_site

def get_tri_ads_site(struct_file,normal_vec,sum_dist,center_coord,center_idx):
#Get adsorption site for 3-coordinate (not trigonal planar)

	dist = get_dist_planar(normal_vec)
	ads_site_planar = get_planar_ads_site(struct_file,center_coord,dist,center_idx)
	NN_planar, mindist_planar = get_NNs(struct_file,ads_site_planar,center_idx)
	ads_site_nonplanar = get_nonplanar_ads_site(sum_dist,center_coord)
	NN_nonplanar, mindist_nonplanar = get_NNs(struct_file,ads_site_nonplanar,center_idx)
	if NN_planar == NN_nonplanar:
		if mindist_planar >= mindist_nonplanar:
			ads_site = ads_site_planar
		else:
			ads_site = ads_site_nonplanar
	elif NN_planar <= NN_nonplanar:
		ads_site = ads_site_planar
	else:
		ads_site = ads_site_nonplanar
		
	return ads_site

def get_opt_ads_site(struct_file,cnum,mic_coords,center_idx,center_coord):
#get the optimal adsorption site

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
			ads_site = center_coord-dist
	elif cnum == 2:
		ads_site = get_bi_ads_site(struct_file,normal_vec,center_coord,center_idx)
	elif cnum == 3 and np.linalg.norm(scaled_sum_dist) > sum_cutoff:
		ads_site = get_tri_ads_site(struct_file,normal_vec,sum_dist,center_coord,center_idx)
	elif norm_scaled <= sum_cutoff or rmse <= rmse_tol:
		dist = get_dist_planar(normal_vec)
		ads_site = get_planar_ads_site(struct_file,center_coord,dist,center_idx)
	else:
		ads_site = get_nonplanar_ads_site(sum_dist,center_coord)

	return ads_site