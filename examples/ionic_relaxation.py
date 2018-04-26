from pymofscreen.cif_handler import get_cif_files
from pymofscreen.screen import screener

#Run screening analysis
mofpath = '/projects/p30148/vasp_jobs/MOFs/testing/structures/'
basepath = '/projects/p30148/vasp_jobs/MOFs/testing/'
kpts_path = '/projects/p30148/vasp_jobs/MOFs/reoptimized_core1/results_cleaned/data_for_phase2/kpts.txt'
submit_script = 'sub_screen4_temp.job'
stdout_file = 'driver.out'
cif_files = get_cif_files(mofpath)
s = screener(mofpath,basepath,submit_script=submit_script,
	stdout_file=stdout_file,kpts_path=kpts_path,niggli=False)
for cif_file in cif_files:
	mof = s.run_screen(cif_file,'ionic')