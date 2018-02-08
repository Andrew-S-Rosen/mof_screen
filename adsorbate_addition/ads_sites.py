import numpy as np
from geom import get_NNs, get_dist_planar
from settings import guess_length, overlap_tol, rcut, coremof_path, ads_species
from ase import Atoms, Atom
from ase.io import read

def get_planar_ads_site(cif_file,cus_coord,dist,ase_cus_idx):
#Get adsorption site for planar structure
	NN = []
	mindist = []
	for i in range(2):
		if i == 0:
			ads_site_temp = cus_coord + dist
		elif i == 1:
			ads_site_temp = cus_coord - dist
		NN_temp, mindist_temp = get_NNs(cif_file,ads_site_temp,ase_cus_idx)
		NN.append(NN_temp)
		mindist.append(mindist_temp)
	if NN[0] == NN[1]:
		if mindist[0] >= mindist[1]:
			ads_site = cus_coord + dist
		else:
			ads_site = cus_coord - dist
	elif NN[0] <= NN[1]:
		ads_site = cus_coord + dist
	else:
		ads_site = cus_coord - dist
	return ads_site

def get_nonplanar_ads_site(sum_dist,cus_coord):
#Get adsorption site for nonplanar structure
	dist = guess_length*sum_dist/np.linalg.norm(sum_dist)
	ads_site =  cus_coord - dist
	return ads_site

def get_bi_ads_site(cif_file,normal_vec,cus_coord,ase_cus_idx):
#Get adsorption site for 2-coordinate geometry
	try_angles = np.arange(0,360,10)
	dist = get_dist_planar(normal_vec)
	ads_site_temp_unrotated1 = cus_coord + dist
	ads_site_temp_unrotated2 = cus_coord - dist
	ads_temp1 = Atoms([Atom(ads_species,ads_site_temp_unrotated1)])
	ads_temp2 = Atoms([Atom(ads_species,ads_site_temp_unrotated2)])
	for i, angle in enumerate(try_angles):
		mof_temp = read(coremof_path+cif_file)
		mof_temp.extend(ads_temp1)
		mof_temp.extend(ads_temp2)
		mof_temp.set_distance(ase_cus_idx,len(mof_temp)-1,guess_length,fix=0,mic=True)
		mof_temp.set_distance(ase_cus_idx,len(mof_temp)-2,guess_length,fix=0,mic=True)	
		mof_temp.set_angle(len(mof_temp)-1,ase_cus_idx,len(mof_temp)-2,angle)
		dist_mat = mof_temp.get_distances(len(mof_temp)-2,np.arange(0,len(mof_temp)-2).tolist(),mic=True)
		NNs = sum(dist_mat <= rcut)
		if i == 0:
			ads_site = mof_temp[-2].position
			old_min_NNs = NNs
		elif sum(dist_mat <= overlap_tol) == 0 and NNs < old_min_NNs:
			ads_site = mof_temp[-2].position
			old_min_NNs = NNs
	return ads_site

def get_tri_ads_site(cif_file,normal_vec,sum_dist,cus_coord,ase_cus_idx):
#Get adsorption site for 3-coordinate (not trigonal planar)
	dist = get_dist_planar(normal_vec)
	ads_site_planar = get_planar_ads_site(cif_file,cus_coord,dist,ase_cus_idx)
	NN_planar, mindist_planar = get_NNs(cif_file,ads_site_planar,ase_cus_idx)
	ads_site_nonplanar = get_nonplanar_ads_site(sum_dist,cus_coord)
	NN_nonplanar, mindist_nonplanar = get_NNs(cif_file,ads_site_nonplanar,ase_cus_idx)
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