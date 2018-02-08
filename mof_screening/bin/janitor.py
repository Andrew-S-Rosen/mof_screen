import os
from settings import basepath
from shutil import copyfile

def prep_paths():
#Folder and file cleanup

	error_path = basepath+'errors'
	results_path = basepath+'results'
	working_path = basepath+'working'
	screen_results_path = basepath+'results/screen_results.dat'
	log_file = 'screening.log'
	if not os.path.exists(error_path):
		os.makedirs(error_path)
	if not os.path.exists(results_path):
		os.makedirs(results_path)
	if not os.path.exists(working_path):
		os.makedirs(working_path)
	if os.path.isfile(screen_results_path) == True:
		open(screen_results_path, 'w').close()
	if os.path.isfile(log_file) == True:
		open(log_file, 'w').close()
	if os.path.isfile('run_vasp.py') == True:
		os.remove('run_vasp.py')

def clean_files(remove_files):
#Clean files

	for file in remove_files:
		if os.path.isfile(file) == True:
			os.remove(file)

def manage_restart_files(file_path):
#Make sure the restart files is copied

	files = ['WAVECAR']
	for file in files:
		if os.path.isfile(file) != True or os.stat(file).st_size == 0:
			if os.path.isfile(file_path+'/'+file) == True:
				copyfile(file_path+'/'+file,os.getcwd()+'/'+file)