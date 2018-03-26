import pymatgen as pm
from pymatgen.analysis.local_env import MinimumVIRENN
import numpy as np
import os

results_path = 'results/'
output_name = 'bad_ads_addition.txt'
nonmetals_list = ['H','He','C','N','O','F','Ne','P','S','Cl','Ar','Se','Br','Kr','I','Xe','Rn']
bad_jobs = []

results = os.listdir(results_path)
results.sort
for folder in results:
	spe_path = results_path+folder+'/final_spe/'
	if os.path.isdir(spe_path):
		for subdir in os.listdir(spe_path):
			dist = []
			full_name = folder+'_'+subdir
			contcar_path = spe_path+subdir+'/CONTCAR'
			struct = pm.Structure.from_file(contcar_path,primitive=False,sort=False)
			nn_object = MinimumVIRENN()
			O_ads_idx = [i for i, atom in enumerate(struct) if atom.species_string == 'O'][-1]
			H_ads_idx = [i for i, atom in enumerate(struct) if atom.species_string == 'H'][-1]
			neighbors = nn_object.get_nn_info(struct,H_ads_idx)
			if neighbors:
				for neighbor in neighbors:
					dist.append(struct.get_distance(H_ads_idx,neighbor['site_index']))
				min_neighbor = neighbors[np.argmin(dist)]
				if min_neighbor['site_index'] != O_ads_idx:
					bad_jobs.append(full_name)
					continue
			else:
				bad_jobs.append(full_name)
				continue
			dist = []
			nn_object = MinimumVIRENN()
			del struct[H_ads_idx]
			O_ads_idx = [i for i, atom in enumerate(struct) if atom.species_string == 'O'][-1]
			neighbors = nn_object.get_nn_info(struct,O_ads_idx-1)
			if neighbors:
				for neighbor in neighbors:
					dist.append(struct.get_distance(O_ads_idx,neighbor['site_index']))
				min_neighbor = neighbors[np.argmin(dist)]
				if min_neighbor['site'].species_string in nonmetals_list:
					bad_jobs.append(full_name)
					continue
			else:
				bad_jobs.append(full_name)
				continue
bad_jobs.sort()
with open(results_path+output_name,'w') as wf:
	for bad_job in bad_jobs:
		wf.write(bad_job+'\n')
