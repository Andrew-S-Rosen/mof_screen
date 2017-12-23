from ase.io import read, write
import numpy as np
import os
from ase import Atoms, Atom

#Paths for files
coremof_path = #path of CIF files to oxygenate
newmofs_path = #path to store oxygenated MOFs
omsdata = #path to .omsex and .oms files from zeo++

#Parameters
guess_length = 2.0 #M-adsorbate bond distance
ads_species = 'O' #adsorbate
rcut = 2.5 #cutoff for nearest neighbors
zeo_tol = 0.1 #tolerance for difference between zeo++ coordinates and ASE
rmse_tol = 0.05 #RMSE tolerance for if a geometry is planar-like
overlap_tol = 1.0

#read in the CIF files for the MOFs from coremof_path
cif_files = []
for filename in os.listdir(coremof_path):
	filename = filename.strip()
	if len(filename.split('.cif')) == 2:
		cif_files.append(filename)

#for each CIF file, add the adsorbate and write the file
for cif_file in cif_files:

	#read in MOF as an ASE Atoms object
	refcode = cif_file.split('.cif')[0]
	mof = read(coremof_path+cif_file)

	#read .oms file to get CN
	f = open(omsdata+refcode+'.oms','r')
	oms_file = f.read()
	n_OMS = int(oms_file.split('OMS=')[1].split('\n')[0])
	f.close()

	#initialize variables
	cnum = np.zeros(n_OMS)
	ads_site = np.zeros((n_OMS,3))
	cus_coord = np.zeros((n_OMS,3))
	ase_cus_idx = np.empty(n_OMS)
	ase_cus_idx = []

	#open .omsex file for MOF
	with open(omsdata+refcode+'.omsex','r') as rf:

		#for each OMS, generate location for adsorbate
		for i, line in enumerate(rf):
			ase_idx = []

			#get the CN from .omsex file
			cnum[i] = int(line.split('CNUM: ')[1].split('|')[0])
			cus_coord[i,:] = np.asarray(np.matrix(line.split('COORD: ')[1].split('|')[0][0:-1]))
			coords = np.asarray(np.matrix(line.split('NN: ')[1][0:-4]))

			#for each NN in .omsex, get the corresponding ASE index
			for j in range(np.shape(coords)[0]):
				nn_pass = False
				for k, element in enumerate(mof):
					if sum(coords[j,:] > element.position-zeo_tol) == 3 and sum(coords[j,:] < element.position+zeo_tol) == 3:
						ase_idx.append(k)
						nn_fail = True
						break
				if nn_fail == False:
					print('WARNING with '+refcode+': a zeo++ NN (#'+str(j)+') is not in same spot as in ASE CIF')	

			ase_idx = np.asarray(ase_idx)

			#for each OMS in .omsex, get the corresponding ASE index
			for j, element in enumerate(mof):
				if sum(cus_coord[i,:] > element.position-zeo_tol) == 3 and sum(cus_coord[i,:] < element.position+zeo_tol) == 3:
					ase_cus_idx.append(j)
					break
			if len(ase_cus_idx) < i+1:
				print('WARNING with '+refcode+': a zeo++ OMS (#'+str(i)+') is not in same spot as in ASE CIF')

			#sum up NN vectors from OMS
			dist_orig = sum(mof.get_distances(ase_cus_idx[i],ase_idx,mic=True,vector=True))

			#fit equation of plane to OMS + NNs
			x = coords[:,0][np.newaxis].T
			y = coords[:,1][np.newaxis].T
			z = coords[:,2][np.newaxis].T
			onevec = np.ones((len(x),1))
			A = np.hstack((x,y,onevec))
			B = z
			try:
				fit = np.squeeze(np.dot(np.linalg.inv(np.dot(A.T,A)),np.dot(A.T,B)))
			except:
				fit = np.squeeze(np.dot(np.linalg.pinv(np.dot(A.T,A)),np.dot(A.T,B)))
			z_fit = fit[0]*x+fit[1]*y+fit[2]
			rmse = (sum((z_fit-z)**2)[0]/len(z))**0.5

			#if M-O bond is very small or planar fit RMSE is small, assume planar
			if np.linalg.norm(dist_orig) < 0.5 or rmse < rmse_tol:

				#calculate unit normal and scale to guess_length
				normal_vec = np.array([fit[0],fit[1],-1])
				unit_normal = normal_vec/np.linalg.norm(normal_vec)
				dist = unit_normal * guess_length
				NN_temp = np.zeros(2)

				#consider both +/- unit normal
				for j in range(2):
					mof_temp = read(coremof_path+cif_file)
					if j == 0:
						ads_site_temp = cus_coord[i,:] + dist
					elif j == 1:
						ads_site_temp = cus_coord[i,:] - dist
					adsorbate_temp = Atoms([Atom(ads_species,ads_site_temp)])
					mof_temp.extend(adsorbate_temp)
					neighbor_dist_temp = mof_temp.get_distances(len(mof_temp)-1,np.arange(0,len(mof_temp)-1).tolist(),mic=True)
					NN_temp[j] = sum(neighbor_dist_temp < rcut)

				#take the +/- sign of unit normal based on fewest number of NN
				if NN_temp[0] <= NN_temp[1]:
					ads_site[i,:] = cus_coord[i,:] + dist
				elif NN_temp[0] > NN_temp[1]:
					ads_site[i,:] = cus_coord[i,:] - dist

			#otherwise scale sum of vectors to 2 
			else:
				dist = 2.0*dist_orig/np.linalg.norm(dist_orig)
				ads_site[i,:] =  cus_coord[i,:] - dist

	#get number of NN for each OMS in MOF
	NN = np.zeros(n_OMS)
	for i in range(n_OMS):
		mof = read(coremof_path+cif_file)
		adsorbate = Atoms([Atom(ads_species,ads_site[i,:])])
		mof.extend(adsorbate)
		neighbor_dist = mof.get_distances(len(mof)-1,np.arange(0,len(mof)-1).tolist(),mic=True)
		NN[i] = sum(neighbor_dist < rcut)

	#keep the OMS that has the smallest NN
	best_idx = np.argmin(NN)

	#add adsorbate to original CIF and save new CIF
	mof = read(coremof_path+cif_file)
	adsorbate = Atoms([Atom(ads_species,ads_site[best_idx,:])])
	mof.extend(adsorbate)
	dist_mat = mof.get_distances(len(mof)-1,np.arange(0,len(mof)-1).tolist(),mic=True)
	if sum(dist_mat <= overlap_tol) > 0:
		print('ERROR with '+refcode+' (CNUM = '+str(cnum[best_idx])+'): adsorbate overlaps with NN')
	else:
		write(newmofs_path+refcode+'_'+ads_species+'.cif',mof)
		print('SUCCESS: '+refcode +' (CNUM = '+str(cnum[best_idx])+')')
