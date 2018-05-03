import os
from shutil import copyfile

def pprint(printstr):
	"""
	Redirects pprint to stdout
	Args:
		printstr (string): string to print to stdout
	"""
	print(printstr)
	with open('screening.log','a') as txtfile:
		txtfile.write(printstr+'\n')

def write_success(workflow,neb=False):
	"""
	Write out the successful job files
	Args:
		workflow (class): pymofscreen.screen_phases.worfklow class
	"""
	spin_level = workflow.spin_level
	acc_level = workflow.acc_levels[workflow.run_i]
	pprint('SUCCESS: '+spin_level+', '+acc_level)
	refcode = workflow.refcode
	basepath = workflow.basepath
	cif_file = workflow.cif_file
	vasp_files = workflow.vasp_files
	success_path = os.path.join(basepath,'results',refcode,acc_level,spin_level)
	if not os.path.exists(success_path):
		os.makedirs(success_path)
	if not neb:
		if acc_level == 'final_spe':
			files_to_copy = vasp_files+['DOSCAR','AECCAR0','AECCAR2']
		elif 'dimer' in acc_level:
			files_to_copy = vasp_files+['DIMCAR','MODECAR','NEWMODECAR']
		else:
			files_to_copy = vasp_files
		for file in files_to_copy:
			if os.path.isfile(file) and os.stat(file).st_size > 0:
				write_to_path = os.path.join(success_path,file)
				copyfile(file,write_to_path)
	elif neb:
		tar_file = 'neb.tar.gz'
		os.system('tar -zcvf '+tar_file+' neb')
		if os.path.isfile(tar_file) and os.stat(tar_file).st_size > 0:
			write_to_path = os.path.join(success_path,tar_file)
			copyfile(tar_file,write_to_path)
		os.remove('neb.tar.gz')
	os.remove(os.path.join(basepath,'working',cif_file))
	if os.path.exists('STOPCAR'):
		os.remove('STOPCAR')

def write_errors(workflow,mof,neb=False):
	"""
	Write out the unsuccesful job files
	Args:
		workflow (class): pymofscreen.screen_phases.worfklow class
		mof (ASE Atoms object): ASE Atoms object
	"""
	spin_level = workflow.spin_level
	acc_level = workflow.acc_levels[workflow.run_i]
	pprint('ERROR: '+spin_level+', '+acc_level+' failed')

	if acc_level != 'scf_test' and 'neb' not in acc_level:
		if mof == None:
			pprint('^ VASP crashed')
		elif mof.calc.scf_converged == False:
			pprint('^ SCF did not converge')
		elif mof.calc.converged == False:
			pprint('^ Convergence not reached')
	refcode = workflow.refcode
	basepath = workflow.basepath
	cif_file = workflow.cif_file
	vasp_files = workflow.vasp_files
	stdout_file = workflow.stdout_file
	error_path = os.path.join(basepath,'errors',refcode,acc_level,spin_level)
	if not os.path.exists(error_path):
		os.makedirs(error_path)
	if not neb:
		if 'dimer' in acc_level:
			files_to_copy = vasp_files+[stdout_file,'DIMCAR','MODECAR','NEWMODECAR']
		else:
			files_to_copy = vasp_files+[stdout_file]
		for file in files_to_copy:
			if os.path.isfile(file) and os.stat(file).st_size > 0:
				write_to_path = os.path.join(error_path,file)
				copyfile(file,write_to_path)
	elif neb:
		tar_file = 'neb.tar.gz'
		os.system('tar -zcvf '+tar_file+' neb')
		if os.path.isfile(tar_file) and os.stat(tar_file).st_size > 0:
			write_to_path = os.path.join(error_path,tar_file)
			copyfile(tar_file,write_to_path)		
	os.remove(os.path.join(basepath,'working',cif_file))
	if os.path.exists('STOPCAR'):
		os.remove('STOPCAR')