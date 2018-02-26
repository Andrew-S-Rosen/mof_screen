import pymatgen as pm
import os
results_path = '/projects/p30148/vasp_jobs/MOFs/oxidized_oms/results/'
rcut = 2.5
ads_species = 'O'
nonmetals_list = ['H','He','C','N','O','O','F','Ne','P','S','Cl','Ar','Se','Br','Kr','I','Xe','Rn']
for refcode in os.listdir(results_path):
	spe_path = results_path+refcode+'/final_spe/'
	if os.path.isdir(spe_path):
		for subdir in os.listdir(spe_path):
			contcar_path = spe_path+subdir+'/CONTCAR'
			struct = pm.Structure.from_file(contcar_path)
			main_site = [atom for atom in struct if atom.species_string == ads_species][-1]
			neighbors = struct.get_neighbors(main_site,rcut)
			if len(neighbors) == 0:
				print(refcode)
			else:
				neighbors_strs = [specie[0].species_string for specie in neighbors]
				dist = [specie[1] for specie in neighbors]
				orig_env = neighbors_strs+dist
				env = list(zip(neighbors_strs,dist))
				env.sort(key=lambda x: x[1])
				if env[0][0] in nonmetals_list:
					print([refcode+'_'+subdir]+orig_env)
