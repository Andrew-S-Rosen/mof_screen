import numpy as np
import os

coremof_path = '/projects/p30148/vasp_jobs/MOFs/reoptimized_core1/results_cifs/oxygenated_reoptimized_oms_cifs/'
old_mofpath = '/projects/p30148/vasp_jobs/MOFs/reoptimized_core1/results/'
cif_files = []
kppas=[100,1000]
ads_species='O'

for filename in os.listdir(coremof_path):
	filename = filename.strip()
	if len(filename.split('.cif')) == 2:
		cif_files.append(filename)
cif_files.sort()
for cif_file in cif_files:
	for kppa in kppas:
		cif_split1 = cif_file.split('_'+ads_species+'_OMS')[0]
		old_cif_name = cif_split1.split('_spin')[0]
		spin = 'spin'+cif_split1.split('_spin')[1]
		if kppa == 100:
			kpts_path = old_mofpath+old_cif_name+'/isif2/'+spin+'/KPOINTS'
			print(old_cif_name)
		elif kppa == 1000:
			kpts_path = old_mofpath+old_cif_name+'/final_spe/'+spin+'/KPOINTS'
		else:
			raise ValueError('Incompatible KPPA with prior runs')
		with open(kpts_path,'r') as rf:
			for i, line in enumerate(rf):
				line = line.strip()
				if i == 2:
					if 'gamma' in line.lower():
						gamma = True
					else:
						gamma = False
				if i == 3:
					kpts = np.squeeze(np.asarray(np.matrix(line))).tolist()
		if len(kpts) != 3:
			raise ValueError('Error parsing KPOINTS file')
		print(kpts)
		print(gamma)
