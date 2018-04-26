import os
from shutil import copyfile

def prep_paths(basepath):
	"""
	Prepare the necessary file paths
	Args:
		basepath (string): full path to base
	"""

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
	if os.path.isfile(screen_results_path):
		open(screen_results_path,'w').close()
	if os.path.isfile(log_file):
		open(log_file,'w').close()
	if os.path.isfile('run_vasp.py'):
		os.remove('run_vasp.py')

def clean_files(remove_files):
	"""
	Remove specified files
	Args:
		remove_files (list of strings): filenames to remove
	"""

	for file in remove_files:
		if os.path.isfile(file):
			os.remove(file)

def manage_restart_files(file_path):
	"""
	Copy restart files to current directory
	Args:
		file_path (string): path restart files
	"""

	files = ['WAVECAR','CHGCAR']
	for file in files:
		full_path = os.path.join(file_path,file)
		if not os.path.isfile(file) or os.stat(file).st_size == 0:
			if os.path.isfile(full_path) and os.stat(full_path).st_size > 0:
				copyfile(full_path,os.path.join(os.getcwd(),file))