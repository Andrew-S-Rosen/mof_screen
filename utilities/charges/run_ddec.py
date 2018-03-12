import os
from shutil import copyfile
import time

results_path = '/projects/p30148/vasp_jobs/MOFs/reoptimized_core1/results/'
job_control_path = '/home/asr731/software/chargemol_09_26_2017/scripts/job_control.txt'
submit_script_path = '/home/asr731/software/chargemol_09_26_2017/scripts/sub_ddec.job'
sub_command = 'msub'
refcodes = os.listdir(results_path)
refcodes.sort()
for refcode in refcodes:
	spe_path = results_path+refcode+'/final_spe/'
	if os.path.isdir(spe_path):
		print(refcode)
		for subdir in os.listdir(spe_path):
			ddec_path = spe_path+subdir+'/ddec/'
			if not os.path.exists(ddec_path):
				os.makedirs(ddec_path)
			if 'VASP_DDEC_analysis.output' in os.listdir(ddec_path):
				continue
			files = ['AECCAR0','AECCAR2','CHGCAR','POTCAR']
			for file in files:
				copyfile(spe_path+subdir+'/'+file,ddec_path+file)
			copyfile(job_control_path,ddec_path+'job_control.txt')
			copyfile(submit_script_path,ddec_path+'sub_ddec.job')
			os.chdir(ddec_path)
			time.sleep(2)
			os.system(sub_command +' sub_ddec.job')
