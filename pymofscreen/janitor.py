import os
from shutil import copyfile, rmtree

def clean_files(remove_files):
	"""
	Remove specified files
	Args:
		remove_files (list of strings): filenames to remove
	"""

	for file in remove_files:
		if os.path.isfile(file):
			os.remove(file)

def vtst_cleanup():
	if os.path.exists('neb'):
		rmtree('neb')
	if os.path.exists('dim'):
		rmtree('dim')
	if os.path.exists('neb.tar.gz'):
		os.remove('neb.tar.gz')
	clean_files(['MODECAR'])

def prep_paths(basepath):
	"""
	Prepare the necessary file paths
	Args:
		basepath (string): full path to base
	"""

	error_path = os.path.join(basepath,'errors')
	results_path = os.path.join(basepath,'results')
	working_path = os.path.join(basepath,'working')
	screen_results_path = os.path.join(basepath,'results','screen_results.dat')
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
	clean_files(['run_vasp.py','neb.tar.gz','AECCAR0','AECCAR1','AECCAR2','CENTCAR','CHG','ase-sort.dat','DIMCAR','DOSCAR','EIGENVAL','IBZKPT','OSZICAR','PCDAT','PROCAR','REPORT','vasprun.xml','XDATCAR','WAVECAR','CHGCAR'])
	vtst_cleanup()

def manage_restart_files(file_path,dimer=False,neb=False,wavechg=True):
	"""
	Copy restart files to current directory
	Args:
		file_path (string): path restart files
	"""

	if wavechg:
		files = ['WAVECAR','CHGCAR']
	else:
		files = []
	if dimer == True and neb == True:
		raise ValueError('Cannot be both NEB and dimer')
	if dimer == True:
		files += ['NEWMODECAR']
	if neb == True:
		files = ['neb.tar.gz']
	for file in files:
		full_path = os.path.join(file_path,file)
		if file == 'NEWMODECAR':
			if os.path.isfile('MODECAR'):
				os.remove('MODECAR')
			copyfile(full_path,'MODECAR')
		else:
			if not os.path.isfile(file) or os.stat(file).st_size == 0:
				if os.path.isfile(full_path) and os.stat(full_path).st_size > 0:
					copyfile(full_path,file)
					if '.tar.gz' in file:
						os.system('tar -zxvf '+file)
					elif '.gz' in file:
						os.system('gunzip '+file)