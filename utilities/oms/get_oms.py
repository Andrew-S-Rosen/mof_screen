import os
mof_path = '/projects/p30148/vasp_jobs/MOFs/reoptimized_core1/results_cifs/reoptimized_cifs/'
zeo_path = '/home/asr731/software/zeo++-0.3/network'
os.chdir(mof_path)
cif_files = []
for filename in os.listdir(mof_path):
	filename = filename.strip()
	if '.cif' in filename:
		os.system(zeo_path+' -omsex '+filename)
