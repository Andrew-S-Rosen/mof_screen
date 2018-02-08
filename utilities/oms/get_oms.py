import os
mof_path = '/projects/p30148/vasp_jobs/MOFs/reoptimized_core1/results_cifs/reoptimized_cifs/'
os.chdir(mof_path)
cif_files = []
for filename in os.listdir(mof_path):
	filename = filename.strip()
	if len(filename.split('.cif')) == 2:
		os.system('/home/asr731/software/zeo++-0.3/network -omsex '+filename)
