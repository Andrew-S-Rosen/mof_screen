import os
from shutil import copyfile

def pprint(printstr):
#Redirect print commands to log file

	print(printstr)
	with open('screening.log','a') as txtfile:
		txtfile.write(printstr+'\n')

def write_success(workflow,acc_level,vasp_files):
#Write success files

	spin_level = workflow.spin_level
	refcode = workflow.refcode
	basepath = workflow.basepath
	cif_file = workflow.cif_file

	pprint('SUCCESS: '+spin_level+', '+acc_level)
	success_path = os.path.join(basepath,'results',refcode,acc_level,spin_level)
	if acc_level == 'final_spe':
		vasp_success_files = vasp_files+['DOSCAR','AECCAR0','AECCAR2']
	else:
		vasp_success_files = vasp_files
	if not os.path.exists(success_path):
		os.makedirs(success_path)
	for file in vasp_success_files:
		if os.path.isfile(file) == True and os.stat(file).st_size > 0:
			write_to_path = os.path.join(success_path,file)
			copyfile(file,write_to_path)
	os.remove(os.path.join(basepath,'working',cif_file))

def write_errors(workflow,acc_level,vasp_files):
#Write error files

	spin_level = workflow.spin_level
	refcode = workflow.refcode
	basepath = workflow.basepath
	cif_file = workflow.cif_file

	pprint('ERROR: '+spin_level+', '+acc_level+' failed')
	error_path = os.path.join(basepath,'errors',refcode,acc_level,spin_level)
	if not os.path.exists(error_path):
		os.makedirs(error_path)
	for file in vasp_files:
		if os.path.isfile(file) == True and os.stat(file).st_size > 0:
			write_to_path = os.path.join(error_path,file)
			copyfile(file,write_to_path)
	os.remove(os.path.join(basepath,'working',cif_file))
	if os.path.exists('STOPCAR'):
		os.remove('STOPCAR')