import os
from shutil import copyfile

results_path = '/projects/p30148/vasp_jobs/MOFs/reoptimized_core1/results/'
submit_script_path = '/home/asr731/software/chargemol_09_26_2017/scripts/sub_bader.job'
sub_command = 'msub'
refcodes = os.listdir(results_path)
refcodes.sort()
for refcode in refcodes:
	spe_path = results_path+refcode+'/final_spe/'
	if os.path.isdir(spe_path):
		print(refcode)
		for subdir in os.listdir(spe_path):
			bader_path = spe_path+subdir+'/bader/'
			if not os.path.exists(bader_path):
				os.makedirs(bader_path)
			files = ['AECCAR0','AECCAR2','CHGCAR']
			for file in files:
				copyfile(spe_path+subdir+'/'+file,bader_path+file)
			copyfile(submit_script_path,bader_path+'sub_bader.job')
			os.chdir(bader_path)
			os.system(sub_command +' sub_bader.job')
