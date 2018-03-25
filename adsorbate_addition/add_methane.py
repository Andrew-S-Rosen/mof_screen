import pandas as pd
from ase.io import read, write
from ase.geometry import get_distances
import numpy as np
from ase.build import molecule
import os

cif_path = 'CIFs/'
grid_path = 'ASCI_grids/'
for cif_name in os.listdir(cif_path):
	mof_name = cif_name.split('.cif')[0]
	CH4 = molecule('CH4')
	cif_name = mof_name+'.cif'
	grid_name = mof_name+'.grid'
	if not os.path.isfile(grid_path+grid_name):
		print('ERROR for ',mof_name,': no grid')
	CH_length = CH4.get_distance(0,1)
	CH_angle = CH4.get_angle(1,0,2)
	CH_dihedral = CH4.get_dihedral(2,1,0,4)
	df = pd.read_csv(grid_path+grid_name,delim_whitespace=True,na_values='?',usecols=[0,1,2,3])
	df.columns = ['x','y','z','e']
	mof = read(cif_path+cif_name)
	O_idx = [atom.index for atom in mof if atom.symbol == 'O'][-1]
	O_pos = mof[O_idx].position
	max_dist = 2.5 #or range between X and Y?
	D,D_len = get_distances([O_pos],df[['x','y','z']].as_matrix(),cell=mof.cell,pbc=mof.pbc)
	D.shape = (-1,3)
	D_len.shape = (-1,)
	df['d'] = D_len
	O_df = df[df['d'] <= max_dist]
	best = O_df.loc[O_df.idxmin()['e']]
	ads_site = [best['x'],best['y'],best['z']]
	CH4[0].position = ads_site
	r_vec = O_pos - ads_site
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
	write(mof_name+'_CH4.cif',mof)
