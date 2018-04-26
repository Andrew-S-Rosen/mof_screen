import os
from shutil import copyfile

def pprint(printstr):
	"""
	Redirects pprint to stdout
	Args:
		printstr (string): string to print to stdout
	"""

	with open('screening.log','a') as txtfile:
		txtfile.write(printstr+'\n')

def write_success(workflow):
	"""
	Write out the successful job files
	Args:
		workflow (class): pymofscreen.screen_phases.worfklow class
		acc_level (string): accuracy level
		vasp_files (list of strings): relevant VASP files
	"""
	spin_level = workflow.spin_level
	refcode = workflow.refcode
	basepath = workflow.basepath
	cif_file = workflow.cif_file
	acc_level = workflow.acc_levels[workflow.run_i]
	vasp_files = workflow.vasp_files

	pprint('SUCCESS: '+spin_level+', '+acc_level)
	success_path = os.path.join(basepath,'results',refcode,acc_level,spin_level)
	if acc_level == 'final_spe':
		vasp_success_files = vasp_files+['DOSCAR','AECCAR0','AECCAR2']
	else:
		vasp_success_files = vasp_files
	if not os.path.exists(success_path):
		os.makedirs(success_path)
	for file in vasp_success_files:
		if os.path.isfile(file) and os.stat(file).st_size > 0:
			write_to_path = os.path.join(success_path,file)
			copyfile(file,write_to_path)
	os.remove(os.path.join(basepath,'working',cif_file))
	if os.path.exists('STOPCAR'):
		os.remove('STOPCAR')

def write_errors(workflow):
	"""
	Write out the unsuccesful job files
	Args:
		workflow (class): pymofscreen.screen_phases.worfklow class
		acc_level (string): accuracy level
		vasp_files (list of strings): relevant VASP files
	"""
	spin_level = workflow.spin_level
	refcode = workflow.refcode
	basepath = workflow.basepath
	cif_file = workflow.cif_file
	acc_level = workflow.acc_levels[workflow.run_i]
	vasp_files = workflow.vasp_files

	pprint('ERROR: '+spin_level+', '+acc_level+' failed')
	error_path = os.path.join(basepath,'errors',refcode,acc_level,spin_level)
	if not os.path.exists(error_path):
		os.makedirs(error_path)
	for file in vasp_files:
		if os.path.isfile(file) and os.stat(file).st_size > 0:
			write_to_path = os.path.join(error_path,file)
			copyfile(file,write_to_path)
	os.remove(os.path.join(basepath,'working',cif_file))
	if os.path.exists('STOPCAR'):
		os.remove('STOPCAR')