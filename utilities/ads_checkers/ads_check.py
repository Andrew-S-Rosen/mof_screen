import pymatgen as pm
from pymatgen.analysis.local_env import MinimumVIRENN
import numpy as np
import os

results_path = '/projects/p30148/vasp_jobs/MOFs/oxidized_oms/results/'
ads_species = 'O'
nonmetals_list = ['H','He','C','N','O','F','Ne','P','S','Cl','Ar','Se','Br','Kr','I','Xe','Rn']
bad_jobs = []
for refcode in os.listdir(results_path):
	spe_path = results_path+refcode+'/final_spe/'
	if os.path.isdir(spe_path):
		for subdir in os.listdir(spe_path):
			dist = []
			full_name = refcode+'_'+subdir
			contcar_path = spe_path+subdir+'/CONTCAR'
			struct = pm.Structure.from_file(contcar_path,primitive=False,sort=False)
			nn_object = MinimumVIRENN()
			ads_idx = [i for i, atom in enumerate(struct) if atom.species_string == ads_species][-1]
			neighbors = nn_object.get_nn_info(struct,ads_idx)
			if neighbors:
				for neighbor in neighbors:
					dist.append(struct.get_distance(ads_idx,neighbor['site_index']))
				if neighbors[np.argmin(dist)]['site'].species_string in nonmetals_list:
					bad_jobs.append(full_name)
bad_jobs.sort()
print(bad_jobs)
