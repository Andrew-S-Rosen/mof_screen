import os
from shutil import copyfile
mof_path = '/projects/p30148/vasp_jobs/MOFs/reoptimized_core1/results_cifs/reoptimized_cifs/'
omsdata = '/projects/p30148/vasp_jobs/MOFs/reoptimized_core1/results_cifs/reoptimized_cifs/OMS_data/'
newmofpath = '/projects/p30148/vasp_jobs/MOFs/reoptimized_core1/results_cifs/reoptimized_oms_cifs/'
os.chdir(mof_path)
cif_files = []
for filename in os.listdir(mof_path):
	filename = filename.strip()
	refcode = filename.split('.cif')[0]
	if len(filename.split('.cif')) == 2:
		f = open(omsdata+refcode+'.omsex')
		if len(f.read()) != 0:
			copyfile(filename,newmofpath+filename)
