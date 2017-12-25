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

def get_omsex_data(line):
#read line in OMSEX file
	cnum = int(line.split('CNUM: ')[1].split('|')[0])
	cus_coord = np.asarray(np.matrix(line.split('COORD: ')[1].split('|')[0][0:-1]))
	coords = np.asarray(np.matrix(line.split('NN: ')[1][0:-4]))
	return cnum, cus_coord, coords

def get_ase_idx(mof,coords):
#get ASE indices
	for i in range(np.shape(coords)[0]):
		nn_fail = False
		for j, element in enumerate(mof):
			if sum(coords[i,:] > element.position-zeo_tol) == 3 and sum(coords[i,:] < element.position+zeo_tol) == 3:
				ase_idx.append(j)
				nn_fail = True
				break
		if nn_fail == False:
			print('WARNING with '+refcode+': a zeo++ NN (#'+str(i)+') is not in same spot as in ASE CIF')
	ase_idx_asarray = np.asarray(ase_idx)	
	return ase_idx_asarray

def get_dist_planar(normal_vec):
#get distance vector of adsorbate from planar CUS site
	unit_normal = normal_vec/np.linalg.norm(normal_vec)
	dist = unit_normal*guess_length
	return dist

def get_best_idx(cif_file,n_OMS,ads_site):
#get best OMS with smallest NNs
	NN = np.zeros(n_OMS)
	for i in range(n_OMS):
		mof = read(coremof_path+cif_file)
		adsorbate = Atoms([Atom(ads_species,ads_site[i,:])])
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

def write_files(refcode,mof,cnum,best_idx):
#Write adsorbed CIF
	dist_mat = mof.get_distances(len(mof)-1,np.arange(0,len(mof)-1).tolist(),mic=True)
	if sum(dist_mat <= overlap_tol) > 0:
		write(newmofs_path+'errors/'+refcode+'_'+ads_species+'.cif',mof)
		print('ERROR with '+refcode+' (CNUM = '+str(cnum[best_idx])+'): adsorbate overlaps with NN')
	else:
		write(newmofs_path+refcode+'_'+ads_species+'.cif',mof)
		print('SUCCESS: '+refcode +' (CNUM = '+str(cnum[best_idx])+')')

def get_vert_vec_norm(refcode,mic_coords):
#Get normal vector if plane is vertical in z
	print('NOTE with '+refcode+': NN form vertical plane. Taking cross-product instead of planar fit')
	vec_norm = 0
	for i in range(2,np.shape(mic_coords)[0]):
		vec_temp = np.cross(mic_coords[1,:]-mic_coords[0,:],mic_coords[i,:]-mic_coords[0,:])
		vec_norm_temp = np.linalg.norm(vec_temp)
		if vec_norm_temp > vec_norm:
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

def get_nonplanar_ads_site(dist_orig,cus_coord):
#Get adsorption site for nonplanar structure
	dist = guess_length*dist_orig/np.linalg.norm(dist_orig)
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

def get_tri_ads_site(cif_file,normal_vec,cus_coord):
#Get adsorption site for 3-coordinate (not trigonal planar)
	dist = get_dist_planar(normal_vec)
	ads_site_planar = get_planar_ads_site(cif_file,cus_coord,dist)
	mof_planar = add_ads_species(cif_file,ads_site_planar)
	NN_planar = get_NNs(mof_planar)
	ads_site_nonplanar = get_nonplanar_ads_site(dist_orig,cus_coord)
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
	if os.path.isfile(newmofs_path+refcode+'_'+ads_species+'.cif') == True:
		print('Previously completed: '+refcode)
		continue
	if os.path.isfile(newmofs_path+'errors/'+refcode+'_'+ads_species+'.cif') == True:
		print('Previously completed w/ overlapping NNs: '+refcode)
		continue
	mof = read(coremof_path+cif_file)
	n_OMS = get_CN(omsdata+refcode+'.oms')
	cnum = np.zeros(n_OMS)
	ads_site = np.zeros((n_OMS,3))
	cus_coord = np.zeros((n_OMS,3))
	ase_cus_idx = np.empty(n_OMS)
	ase_cus_idx = []
	with open(omsdata+refcode+'.omsex','r') as rf:
		for i, line in enumerate(rf):
			ase_idx = []
			cnum[i], cus_coord[i,:], coords = get_omsex_data(line)
			ase_idx = get_ase_idx(mof,coords)
			for j, element in enumerate(mof):
				if sum(cus_coord[i,:] > element.position-zeo_tol) == 3 and sum(cus_coord[i,:] < element.position+zeo_tol) == 3:
					ase_cus_idx.append(j)
					break
			if len(ase_cus_idx) < i+1:
				print('WARNING with '+refcode+': a zeo++ OMS (#'+str(i)+') is not in same spot as in ASE CIF')
			mic_coords = mof.get_distances(ase_cus_idx[i],ase_idx,mic=True,vector=True)
			dist_orig = sum(mic_coords)
			fit, rmse, r2, normal_vec = fit_plane(mic_coords)
			if cnum[i] == 2:
				ads_site[i,:] = get_bi_ads_site(cif_file,normal_vec,cus_coord[i,:],mic_coords,ase_cus_idx[i])
			if cnum[i] == 3 and np.linalg.norm(dist_orig) >= sum_cutoff:
				ads_site[i,:] = get_tri_ads_site(cif_file,normal_vec,cus_coord[i,:])
			elif r2 > r2_tol and rmse < rmse_tol:
				if fit[2] >= 1e10:
					print('Doing cross-product for vertical plane')
					vec_norm = get_vert_vec_norm(refcode,mic_coords)
				dist = get_dist_planar(normal_vec)
				ads_site[i,:] = get_planar_ads_site(cif_file,cus_coord[i,:],dist)
			else:
				ads_site[i,:] = get_nonplanar_ads_site(dist_orig,cus_coord[i,:])
	best_idx = get_best_idx(cif_file,n_OMS,ads_site)
	mof = add_ads_species(cif_file,ads_site[best_idx,:])
	write_files(refcode,mof,cnum,best_idx)
