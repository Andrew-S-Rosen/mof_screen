from ase.io import read, write
import numpy as np
import os
from ase import Atoms, Atom

#Paths for files
coremof_path = 'C:/Users/asros/OneDrive/Working/coremof_path/' #path for MOFs to oxygenate
newmofs_path = 'C:/Users/asros/OneDrive/Working/newmof_path/' #path to generate oxygenated MOFs
omsdata = 'C:/Users/asros/OneDrive/Working/oms_data/' #path to .omsex and .oms files

#Parameters
guess_length = 2.0 #M-adsorbate bond distance
ads_species = 'O' #adsorbate
rcut = 2.5 #cutoff for nearest neighbors
zeo_tol = 0.1 #tolerance for difference between zeo++ coordinates and ASE
sum_cutoff = 0.5 #cutoff for sum of vectors to be planar
rmse_tol = 0.25 #rmse tolerance for if a geometry is planar-like
r2_tol = 0.95 #r^2 tolerance for if a geometry is planar-like
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
	if not os.path.exists(newmofs_path+'errors/'):
		os.makedirs(newmofs_path+'errors/')
	return cif_files

def fit_plane(mic_coords):
#fit equation of plane and get fit parameters and RMSE
	x = mic_coords[:,0][np.newaxis].T
	y = mic_coords[:,1][np.newaxis].T
	z = mic_coords[:,2][np.newaxis].T
	onevec = np.ones((len(x),1))
	A = np.hstack((x,y,onevec))
	B = z
	fit = np.squeeze(np.dot(np.linalg.pinv(np.dot(A.T,A)),np.dot(A.T,B)))
	z_fit = fit[0]*x+fit[1]*y+fit[2]
	ss_res = sum((z_fit-z)**2)[0]
	rmse = (ss_res/len(z))**0.5
	ss_tot = sum((z-np.mean(z))**2)[0]
	r2 = 1-ss_res/ss_tot
	normal_vec = np.array([fit[0],fit[1],-1])
	return fit, rmse, r2, normal_vec

def get_CN(oms_path):
#read OMS file for n_OMS
	f = open(oms_path,'r')
	oms_file = f.read()
	n_OMS = int(oms_file.split('OMS=')[1].split('\n')[0])
	f.close()
	return n_OMS

def get_omsex_line(line):
#read line in OMSEX file
	oms_symbol = line.split(' |')[0]
	cnum = int(line.split('CNUM: ')[1].split('|')[0])
	cus_coord = np.asarray(np.matrix(line.split('COORD: ')[1].split('|')[0][0:-1]))
	coords = np.asarray(np.matrix(line.split('NN: ')[1][0:-4]))
	return oms_symbol, cnum, cus_coord, coords

def get_omsex_data(refcode):
#Get all OMSEX data
	cus_coords_all = np.zeros((n_OMS,3))
	ase_cus_idx_all = []
	oms_sym_all = []
	cnums_all = []
	with open(omsdata+refcode+'.omsex','r') as rf:
		for i, line in enumerate(rf):
			oms_sym_temp, cnum_temp, cus_coords_all[i,:], NN_coords_temp = get_omsex_line(line)
			oms_sym_all.append(oms_sym_temp)
			cnums_all.append(cnum_temp)
			if i == 0:
				NN_coords_all = NN_coords_temp
			else:
				NN_coords_all = np.vstack((NN_coords_all,NN_coords_temp))
			for j, element in enumerate(mof):
				if sum(cus_coords_all[i,:] > element.position-zeo_tol) == 3 and sum(cus_coords_all[i,:] < element.position+zeo_tol) == 3:
					ase_cus_idx_all.append(j)
					break
			if len(ase_cus_idx_all) < i+1:
				print('WARNING with '+refcode+': a zeo++ OMS (#'+str(i)+') is not in same spot as in ASE CIF')
		ase_cus_idx_all = np.array(ase_cus_idx_all)
		oms_sym_all = np.array(oms_sym_all)
	return cnums_all, cus_coords_all, ase_cus_idx_all, oms_sym_all, NN_coords_all

def get_ase_NN_idx(mof,coords):
#get ASE indices for NN
	ase_NN_idx = []
	for i in range(np.shape(coords)[0]):
		nn_fail = False
		for j, element in enumerate(mof):
			if sum(coords[i,:] > element.position-zeo_tol) == 3 and sum(coords[i,:] < element.position+zeo_tol) == 3:
				ase_NN_idx.append(j)
				nn_fail = True
				break
		if nn_fail == False:
			print('WARNING with '+refcode+': a zeo++ NN (#'+str(i)+') is not in same spot as in ASE CIF')
	ase_idx_asarray = np.asarray(ase_NN_idx)	
	return ase_idx_asarray

def get_dist_planar(normal_vec):
#get distance vector of adsorbate from planar CUS site
	unit_normal = normal_vec/np.linalg.norm(normal_vec)
	dist = unit_normal*guess_length
	return dist

def get_best_idx(cif_file,ads_sites):
#get best OMS with smallest NNs
	NN = np.zeros(np.shape(ads_sites)[0])
	for i in range(len(NN)):
		mof = read(coremof_path+cif_file)
		adsorbate = Atoms([Atom(ads_species,ads_sites[i,:])])
		mof.extend(adsorbate)
		neighbor_dist = mof.get_distances(len(mof)-1,np.arange(0,len(mof)-1).tolist(),mic=True)
		NN[i] = sum(neighbor_dist < rcut)
	best_idx = np.argmin(NN)
	return best_idx

def add_ads_species(cif_file,ads_site):
#add adsorbate to original CIF and save new CIF
	mof = read(coremof_path+cif_file)
	adsorbate = Atoms([Atom(ads_species,ads_site)])
	mof.extend(adsorbate)
	return mof

def write_files(refcode,mof,oms_sym,cnum,best_idx,i):
#Write adsorbed CIF
	dist_mat = mof.get_distances(len(mof)-1,np.arange(0,len(mof)-1).tolist(),mic=True)
	basename = refcode+'_'+ads_species
	if sum(dist_mat <= overlap_tol) > 0:
		error_path = newmofs_path+'errors/'+basename
		if not os.path.exists(error_path):
			os.makedirs(error_path)
		write(error_path+'/'+basename+'_v'+str(i)+'.cif',mof)
		print('ERROR with '+refcode+'_v'+str(i)+' (M = '+oms_sym+', CNUM = '+str(cnum)+'): adsorbate overlaps with NN')
	else:
		result_path = newmofs_path+basename
		if not os.path.exists(result_path):
			os.makedirs(result_path)
		write(result_path+'/'+basename+'_v'+str(i)+'.cif',mof)
		print('SUCCESS: '+refcode +'_v'+str(i)+' (M = '+oms_sym+', CNUM = '+str(cnum)+')')

def get_vert_vec_norm(mic_coords):
#Get normal vector if plane is vertical in z
	vec_norm = np.inf
	for i in range(2,np.shape(mic_coords)[0]):
		vec_temp = np.cross(mic_coords[1,:]-mic_coords[0,:],mic_coords[i,:]-mic_coords[0,:])
		vec_norm_temp = np.linalg.norm(vec_temp)
		if np.abs(vec_temp[-1])/vec_norm_temp < vec_norm:
			normal_vec = vec_temp
			vec_norm = vec_norm_temp
	return normal_vec

def get_NNs(mof):
#Get number of NNs from last atom (presumably adsorbate)
	neighbor_dist_temp = mof.get_distances(len(mof)-1,np.arange(0,len(mof)-1).tolist(),mic=True)
	NN = sum(neighbor_dist_temp < rcut)
	return NN

def get_planar_ads_site(cif_file,cus_coord,dist):
#Get adsorption site for planar structure
	NN_temp = np.zeros(2)
	for i in range(2):
		mof_temp = read(coremof_path+cif_file)
		if i == 0:
			ads_site_temp = cus_coord + dist
		elif i == 1:
			ads_site_temp = cus_coord - dist
		adsorbate_temp = Atoms([Atom(ads_species,ads_site_temp)])
		mof_temp.extend(adsorbate_temp)
		NN_temp[i] = get_NNs(mof_temp)
	#take the +/- sign of unit normal based on fewest number of NN
	if NN_temp[0] <= NN_temp[1]:
		ads_site = cus_coord + dist
	elif NN_temp[0] > NN_temp[1]:
		ads_site = cus_coord - dist
	return ads_site

def get_nonplanar_ads_site(sum_dist,cus_coord):
#Get adsorption site for nonplanar structure
	dist = guess_length*sum_dist/np.linalg.norm(sum_dist)
	ads_site =  cus_coord - dist
	return ads_site

def get_bi_ads_site(cif_file,normal_vec,cus_coord,mic_coords,ase_cus_idx):
#Get adsorption site for 2-coordinate
	try_angles = np.linspace(0,360,37)
	dist = get_dist_planar(normal_vec)
	ads_site_temp_unrotated = cus_coord + dist
	for i, angle in enumerate(try_angles):
		mof_temp = read(coremof_path+cif_file)
		adsorbate_temp = Atoms([Atom(ads_species,ads_site_temp_unrotated)])
		adsorbate_temp.rotate(angle,mic_coords[1,:]-mic_coords[0,:])
		mof_temp.extend(adsorbate_temp)
		mof_temp.set_distance(ase_cus_idx,len(mof_temp)-1,guess_length,fix=0,mic=False)
		dist_mat = mof_temp.get_distances(len(mof_temp)-1,np.arange(0,len(mof_temp)-1).tolist(),mic=True)
		NNs = sum(dist_mat < rcut)
		if i == 0:
			ads_site = mof_temp[-1].position
			old_min_NNs = NNs
		elif sum(dist_mat <= overlap_tol) == 0 and NNs < old_min_NNs:
			ads_site = mof_temp[-1].position
			old_min_NNs = NNs
	return ads_site

def get_tri_ads_site(cif_file,normal_vec,sum_dist,cus_coord):
#Get adsorption site for 3-coordinate (not trigonal planar)
	dist = get_dist_planar(normal_vec)
	ads_site_planar = get_planar_ads_site(cif_file,cus_coord,dist)
	mof_planar = add_ads_species(cif_file,ads_site_planar)
	NN_planar = get_NNs(mof_planar)
	ads_site_nonplanar = get_nonplanar_ads_site(sum_dist,cus_coord)
	mof_nonplanar = add_ads_species(cif_file,ads_site_nonplanar)
	NN_nonplanar = get_NNs(mof_nonplanar)
	if NN_planar <= NN_nonplanar:
		ads_site = ads_site_planar
	else:
		ads_site = ads_site_nonplanar
	return ads_site

cif_files = get_cif_files()
for cif_file in cif_files:
	refcode = cif_file.split('.cif')[0]
	basename = refcode+'_'+ads_species
	mof = read(coremof_path+cif_file)
	n_OMS = get_CN(omsdata+refcode+'.oms')
	cnums_all, cus_coords_all, ase_cus_idx_all, oms_sym_all, NN_coords_all = get_omsex_data(refcode)
	unique_oms_sym = np.unique(oms_sym_all)
	v = 0
	for oms_sym in unique_oms_sym:
		oms_idx = np.where(oms_sym_all == oms_sym)			
		unique_cnums = np.unique(np.array(cnums_all)[oms_idx])
		for cnum in unique_cnums:
			cnum_idx = np.where(cnums_all == cnum)
			intersect_idx = np.intersect1d(oms_idx,cnum_idx)
			ads_sites = np.zeros((len(intersect_idx),3))
			for i, omsex_idx in enumerate(intersect_idx):
				file = basename+'_v'+str(v)+'.cif'
				if os.path.isfile(newmofs_path+file) == True:
					print('Previously completed: '+file)
					continue
				if os.path.isfile(newmofs_path+'errors/'+file) == True:
					print('Previously completed w/ overlapping NNs: '+file)
					continue
				sum_prior_cnums = sum(cnums_all[0:i])
				NN_coords = NN_coords_all[sum_prior_cnums:sum_prior_cnums+cnum,:]
				cus_coords = cus_coords_all[omsex_idx,:]
				ase_cus_idx = ase_cus_idx_all[omsex_idx]
				ase_NN_idx = get_ase_NN_idx(mof,NN_coords)
				mic_coords = mof.get_distances(ase_cus_idx,ase_NN_idx,mic=True,vector=True)
				scaled_mic_coords = mic_coords*guess_length/np.linalg.norm(mic_coords,axis=1)[np.newaxis].T
				scaled_sum_dist = sum(scaled_mic_coords)
				sum_dist = sum(mic_coords)
				fit, rmse, r2, normal_vec = fit_plane(mic_coords)
				if cnum == 2:
					ads_sites[i,:] = get_bi_ads_site(cif_file,normal_vec,cus_coords,mic_coords,ase_cus_idx)
				if cnum == 3 and np.linalg.norm(scaled_sum_dist) >= sum_cutoff:
					ads_sites[i,:] = get_tri_ads_site(cif_file,normal_vec,sum_dist,cus_coords)
				elif np.linalg.norm(scaled_sum_dist) < sum_cutoff or (r2 > r2_tol and rmse < rmse_tol):
					if fit[2] >= 1e10 or rmse >= rmse_tol*2:
						normal_vec = get_vert_vec_norm(mic_coords)
					dist = get_dist_planar(normal_vec)
					ads_sites[i,:] = get_planar_ads_site(cif_file,cus_coords,dist)
				else:
					ads_sites[i,:] = get_nonplanar_ads_site(sum_dist,cus_coords)
			best_idx = get_best_idx(cif_file,ads_sites)
			mof = add_ads_species(cif_file,ads_sites[best_idx,:])
			write_files(refcode,mof,oms_sym,cnum,best_idx,v)
			v += 1
