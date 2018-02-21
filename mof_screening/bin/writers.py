import os
from settings import basepath
from shutil import copyfile
from ase.io import read

def pprint(printstr):
#Redirect print commands to log file

	print(printstr)
	with open('screening.log','a') as txtfile:
		txtfile.write(printstr+'\n')

def write_success(refcode,spin_level,acc_level,vasp_files,cif_file):
#Write success files

	pprint('SUCCESS: '+spin_level+', '+acc_level)
	success_path = basepath+'results/'+refcode+'/'+acc_level+'/'+spin_level
	if not os.path.exists(success_path):
		os.makedirs(success_path)
	for file in vasp_files:
		if os.path.isfile(file) == True:
			write_to_path = success_path+'/'+file
			copyfile(file,write_to_path)
	os.remove(basepath+'working/'+cif_file)

def write_errors(refcode,spin_level,acc_level,vasp_files,cif_file):
#Write error files

	pprint('ERROR: '+spin_level+', '+acc_level+' failed')
	error_path = basepath+'errors/'+refcode+'/'+acc_level+'/'+spin_level
	if not os.path.exists(error_path):
		os.makedirs(error_path)
	for file in vasp_files:
		if os.path.isfile(file) == True:
			write_to_path = error_path+'/'+file
			copyfile(file,write_to_path)
	os.remove(basepath+'working/'+cif_file)

def write_energy(refcode,acc_level,spin_level):
#Write energy to results file

	outcarpath = basepath+'results/'+refcode+'/'+acc_level+'/'+spin_level+'/OUTCAR'
	final_mof = read(outcarpath)
	E = final_mof.get_potential_energy()
	str_to_write = refcode+'_'+spin_level+': '+str(E)
	results_file = basepath+'results/screen_results.dat'
	with open(results_file,'a') as txtfile:
		txtfile.write(str_to_write+'\n')