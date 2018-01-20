from ase.io import read, write
import numpy as np
import os
from ase import Atoms, Atom

#Paths for files
coremof_path = '/projects/p30148/vasp_jobs/MOFs/reoptimized_core1/results_cifs/reoptimized_cifs/' #path for MOFs to oxygenate
newmofs_path = '/projects/p30148/vasp_jobs/MOFs/reoptimized_core1/results_cifs/oxygenated_reoptimized_oms_cifs/' #path to generate oxygenated MOFs
omsdata = '/projects/p30148/vasp_jobs/MOFs/reoptimized_core1/results_cifs/reoptimized_cifs/OMS_data/' #path to .omsex and .oms files
error_path = newmofs_path+'errors/'

#Parameters
guess_length = 2.0 #M-adsorbate bond distance
ads_species = 'O' #adsorbate
rcut = 2.5 #cutoff for nearest neighbors
zeo_tol = 0.1 #tolerance for difference between zeo++ coordinates and ASE
sum_cutoff = 0.5 #cutoff for sum of vectors to be planar
rmse_tol = 0.25 #rmse tolerance for if a geometry is planar-like
overlap_tol = 1.3 #tolerance for overlapping atoms

def get_cif_files():
#read in the CIF files for the MOFs from coremof_path
	cif_files = []
	for filename in os.listdir(coremof_path):
		filename = filename.strip()
		if len(filename.split('.cif')) == 2:
			cif_files.append(filename)
	if not os.path.exists(newmofs_path):
		os.makedirs(newmofs_path)
	if not os.path.exists(error_path):
		os.makedirs(error_path)
	return cif_files

def fit_line(xyz):
	x = mic_coords[:,0][np.newaxis].T
	y = mic_coords[:,1][np.newaxis].T
	z = mic_coords[:,2][np.newaxis].T
	onevec = np.ones((len(x),1))
	A = np.hstack((x,y,onevec))
	B = z
	fit = np.squeeze(np.dot(np.linalg.pinv(np.dot(A.T,A)),np.dot(A.T,B)))
	z_fit = fit[0]*x+fit[1]*y+fit[2]
	ss_res = sum((z_fit-z)**2)[0]
	ss_tot = sum((z-np.mean(z))**2)[0]
	r2 = 1-ss_res/ss_tot
	normal_vec = np.array([fit[0],fit[1],-1])
	if r2 < 1:
		print('WARNING: Poor linear fit to two points')
	return normal_vec

def fit_plane(xyz):
#orthogonal regression to ax+by+cz+d=0
	xyz_mean = np.mean(xyz,axis=0)
	xyz_sub = xyz-xyz_mean
	[u,s,v] = np.linalg.svd(xyz_sub,full_matrices=False)
	v = v.T
	normal_vec = v[:,-1]
	a = normal_vec[0]
	b = normal_vec[1]
	c = normal_vec[2]
	d = -np.dot(xyz_mean,normal_vec)
	fit = a*xyz[:,0]+b*xyz[:,1]+c*xyz[:,2]+d
	ss_res = sum(fit**2)
	rmse = (ss_res/np.shape(xyz)[0])**0.5
	return rmse, normal_vec

def get_CN(oms_path):
#read OMS file for n_OMS
	f = open(oms_path,'r')
	oms_file = f.read()
	n_OMS = int(oms_file.split('OMS=')[1].split('\n')[0])
	f.close()
	return n_OMS

def get_omsex_line(line):
#read line in OMSEX file
	cus_symbol = line.split(' |')[0]
	cnum = int(line.split('CNUM: ')[1].split('|')[0])
	cus_coord = np.asarray(np.matrix(line.split('COORD: ')[1].split('|')[0][0:-1]))
	coords = np.asarray(np.matrix(line.split('NN: ')[1][0:-4]))
	return cus_symbol, cnum, cus_coord, coords

def get_omsex_data(refcode):
#Get all OMSEX data
	cus_coords_all = np.zeros((n_OMS,3))
	ase_cus_idx_all = []
	cus_sym_all = []
	cnums_all = []
	with open(omsdata+refcode+'.omsex','r') as rf:
		for i, line in enumerate(rf):
			cus_sym_temp, cnum_temp, cus_coords_all[i,:], NN_coords_temp = get_omsex_line(line)
			cus_sym_all.append(cus_sym_temp)
			cnums_all.append(cnum_temp)
			if i == 0:
				NN_coords_all = NN_coords_temp
			else:
				NN_coords_all = np.vstack((NN_coords_all,NN_coords_temp))
			for j, element in enumerate(mof):
				if sum(cus_coords_all[i,:] >= element.position-zeo_tol) == 3 and sum(cus_coords_all[i,:] <= element.position+zeo_tol) == 3:
					ase_cus_idx_all.append(j)
					break
			if len(ase_cus_idx_all) < i+1:
				print('WARNING with '+refcode+': a zeo++ OMS (#'+str(i)+') is not in same spot as in ASE CIF')
	return cnums_all, cus_coords_all, ase_cus_idx_all, cus_sym_all, NN_coords_all

def get_ase_NN_idx(mof,coords):
#get ASE indices for NN
	ase_NN_idx = []
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
#get distance vector of adsorbate from planar CUS site
	unit_normal = normal_vec/np.linalg.norm(normal_vec)
	dist = unit_normal*guess_length
	return dist

def get_NNs(cif_file,ads_site,ase_cus_idx):
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
#sort the OMS by smallest NNs
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

def add_ads_species(cif_file,ads_site):
#add adsorbate to original CIF and save new CIF
	mof_temp = read(coremof_path+cif_file)
	adsorbate = Atoms([Atom(ads_species,ads_site)])
	mof_temp.extend(adsorbate)
	return mof_temp

def write_files(refcode,mof,best_to_worst_idx,cluster):
#Write adsorbed CIF
	basename = refcode+'_'+ads_species
	success = False
	for idx in best_to_worst_idx:
		mof = add_ads_species(cif_file,ads_sites[idx,:])
		dist_mat = mof.get_distances(len(mof)-1,np.arange(0,len(mof)-1).tolist(),mic=True)
		if sum(dist_mat <= overlap_tol) == 0:
			print('SUCCESS: '+refcode+' ('+str(cluster)+')')
			write(newmofs_path+basename+'_OMS'+str(idx)+'.cif',mof)
			success = True
			break
		else:
			del mof[-1]
	if success == False:
		print('ERROR: '+refcode+' ('+str(cluster)+')')
		write(error_path+basename+'_'+str(cluster)+'.cif',mof)

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
	#take the +/- sign of unit normal based on fewest number of NN
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
#Get adsorption site for 2-coordinate
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
	tri_planar = True
	if NN_planar == NN_nonplanar:
		if mindist_planar >= mindist_nonplanar:
			ads_site = ads_site_planar
		else:
			ads_site = ads_site_nonplanar
	elif NN_planar <= NN_nonplanar:
		ads_site = ads_site_planar
	else:
		ads_site = ads_site_nonplanar
		tri_planar = False
	return ads_site, tri_planar

cif_files = get_cif_files()
for cif_file in cif_files:
	refcode = cif_file.split('.cif')[0]
	if os.stat(omsdata+refcode+'.omsex').st_size == 0:
		continue
	basename = refcode+'_'+ads_species
	mof = read(coremof_path+cif_file)
	n_OMS = get_CN(omsdata+refcode+'.oms')
	cnums_all, cus_coords_all, ase_cus_idx_all, cus_sym_all, NN_coords_all = get_omsex_data(refcode)
	cluster_sym = []
	for i, ase_cus_idx in enumerate(ase_cus_idx_all):
		if mof[ase_cus_idx].symbol != cus_sym_all[i]:
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
			ase_NN_idx = get_ase_NN_idx(mof,NN_coords)
			mic_coords = mof.get_distances(ase_cus_idx,ase_NN_idx,mic=True,vector=True)
			scaled_mic_coords = mic_coords*guess_length/np.linalg.norm(mic_coords,axis=1)[np.newaxis].T
			scaled_sum_dist = sum(scaled_mic_coords)
			sum_dist = sum(mic_coords)
			norm_scaled = np.linalg.norm(scaled_sum_dist)
			if cnum >= 3:
				rmse, normal_vec = fit_plane(mic_coords)
			else:
				normal_vec = fit_line(mic_coords)
			if cnum == 1:
				raise ValueError('Not coded!')
			elif cnum == 2:
				ads_sites[i,:] = get_bi_ads_site(cif_file,normal_vec,cus_coords,ase_cus_idx)
				print(refcode+': using angular sweep (cnum='+str(cnum)+')')
			elif cnum == 3 and np.linalg.norm(scaled_sum_dist) > sum_cutoff:
				ads_sites[i,:], tri_planar = get_tri_ads_site(cif_file,normal_vec,sum_dist,cus_coords,ase_cus_idx)
				if tri_planar == True:
					print(refcode+': using least-squares plane (cnum='+str(cnum)+', RMSE='+str(np.round(rmse,2))+', sum(r_i)='+str(np.round(norm_scaled,2))+')')
				else:
					print(refcode+': using sum of vectors (cnum='+str(cnum)+', RMSE='+str(np.round(rmse,2))+', sum(r_i)='+str(np.round(norm_scaled,2))+')')
			elif norm_scaled <= sum_cutoff or rmse <= rmse_tol:
				print(refcode+': using least-squares plane (cnum='+str(cnum)+', RMSE='+str(np.round(rmse,2))+', sum(r_i)='+str(np.round(norm_scaled,2))+')')
				dist = get_dist_planar(normal_vec)
				ads_sites[i,:] = get_planar_ads_site(cif_file,cus_coords,dist,ase_cus_idx)
			else:
				ads_sites[i,:] = get_nonplanar_ads_site(sum_dist,cus_coords)
				print(refcode+': using sum of vectors (cnum='+str(cnum)+', RMSE='+str(np.round(rmse,2))+', sum(r_i)='+str(np.round(norm_scaled,2))+')')
		ase_cus_idx_cluster = [ase_cus_idx_all[j] for j in omsex_indices]
		best_to_worst_idx = get_best_to_worst_idx(cif_file,ads_sites,ase_cus_idx_cluster)
		write_files(refcode,mof,best_to_worst_idx,unique_cluster_sym)
		indices_tot += len(omsex_indices)
	if indices_tot != len(cus_sym_all):
		raise ValueError('Did not run through all OMS')
