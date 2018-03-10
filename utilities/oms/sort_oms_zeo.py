import os
from shutil import copyfile
mof_path = '/projects/p30148/vasp_jobs/MOFs/reoptimized_core1/results_cifs/reoptimized_cifs/'
omsdata = '/projects/p30148/vasp_jobs/MOFs/reoptimized_core1/results_cifs/reoptimized_cifs/OMS_data/'
newmofpath = '/projects/p30148/vasp_jobs/MOFs/reoptimized_core1/results_cifs/reoptimized_oms_cifs/'
cif_files = []
for filename in os.listdir(mof_path):
	filename = filename.strip()
	refcode = filename.split('.cif')[0]
	if '.cif' in filename:
		if os.stat(omsdata+refcode+'.omsex').st_size > 0:
			copyfile(mof_path+filename,newmofpath+filename)
