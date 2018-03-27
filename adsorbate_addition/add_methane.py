import pandas as pd
from ase.io import read, write
from ase.geometry import get_distances
import numpy as np
from ase.build import molecule
import os

#Settings
cif_path = 'CIFs/'
grid_path = 'ASCI_Grids/'
new_mof_path = 'CH4_ads/'
max_dist = 3.0
overlap_tol = 1.3
mol_species = 'CH4'

#Get CH4 parameters
CH4 = molecule(mol_species)
CH_length = CH4.get_distance(0,1)
CH_angle = CH4.get_angle(1,0,2)
CH_dihedral = CH4.get_dihedral(2,1,0,4)
n_CH4 = len(CH4)

#Prep paths
error_path = new_mof_path+'errors/'
if not os.path.exists(new_mof_path):
	os.makedirs(new_mof_path)
if not os.path.exists(error_path):
	os.makedirs(error_path)
completed = os.listdir(new_mof_path)

partition = 1e6
for cif_name in os.listdir(cif_path):

	#prep paths
	overlap = False
	mof_name = cif_name.split('.cif')[0]
	if mof_name+'_CH4.cif' in completed:
		continue
	grid_name = mof_name+'.grid'
	if not os.path.isfile(grid_path+grid_name):
		print('ERROR for',mof_name,': no grid')
		continue

	#read grid
	df = pd.read_csv(grid_path+grid_name,delim_whitespace=True,na_values='?',usecols=[0,1,2,3])
	df.columns = ['x','y','z','e']
	df['d'] = ''

	#read MOF
	mof = read(cif_path+cif_name)
	O_idx = [atom.index for atom in mof if atom.symbol == 'O'][-1]
	O_pos = mof[O_idx].position

	#Construct new df based on max_dist
	n_loops = int(np.ceil(len(df)/partition))
	O_df = pd.DataFrame()
	for i in range(n_loops):
		if i == n_loops-1:
			idx = np.arange(i*int(partition),len(df))
		else:
			idx = np.arange(i*int(partition),(i+1)*int(partition))
		D,D_len = get_distances([O_pos],df.loc[idx,['x','y','z']].as_matrix(),cell=mof.cell,pbc=mof.pbc)
		D_len.shape = (-1,)
		df.loc[idx,'d'] = D_len
	O_df = df[df['d'] <= max_dist]

	#Get lowest energy spot for CH4
	best = O_df.loc[O_df.idxmin()['e']]
	ads_site = [best['x'],best['y'],best['z']]

	#Construct CH4 and add to MOF
	CH4 = molecule(mol_species)
	CH4[0].position = ads_site
	D,D_len = get_distances([ads_site],[O_pos],cell=mof.cell,pbc=mof.pbc)
	r_vec = D[0,0]
	r = (r_vec/np.linalg.norm(r_vec))*CH_length
	CH4[1].position = ads_site+r
	CH4.set_distance(0,2,CH_length,fix=0)
	CH4.set_angle(1,0,2,CH_angle)
	CH4.set_distance(0,3,CH_length,fix=0)
	CH4.set_angle(1,0,3,CH_angle)
	CH4.set_dihedral(2,1,0,3,-CH_dihedral)
	CH4.set_distance(0,4,CH_length,fix=0)
	CH4.set_angle(1,0,4,CH_angle)
	CH4.set_dihedral(2,1,0,4,CH_dihedral)
	mof.extend(CH4)

	#Confirm no overlapping atoms
	for i in range(n_CH4):
		dist = mof.get_distances(-(i+1),np.arange(0,len(mof)-n_CH4).tolist(),mic=True)
		if np.sum(dist <= overlap_tol) > 0:
			print('ERROR for',mof_name,': overlap')
			overlap = True
			write(error_path+mof_name+'_CH4.cif',mof)
			break
	if overlap == True:
		continue

	#write new CIF
	write(new_mof_path+mof_name+'_CH4.cif',mof)
