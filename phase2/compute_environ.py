import os
import settings

def get_nprocs():
#Get the number of CPUs
	submit_script = settings.submit_script
	with open(submit_script,'r') as rf:
		for line in rf:
			if 'nodes' in line or 'ppn' in line:
				line = line.strip().replace(' ','')
				nodes = int(line.split('nodes=')[1].split(':ppn=')[0])
				ppn = int(line.split('nodes=')[1].split(':ppn=')[1])
	nprocs = nodes*ppn
	return nprocs, ppn

def choose_vasp_version(gpt_version,nprocs,calc_swaps,*fix_version):
#Run a given VASP version
	vasp_ex = ['vasp_std','vasp_gam']
	vtst_ex = ['vasp_std_vtst','vasp_gam_vtst']
	vasp_cmd =  'mpirun -n '+str(nprocs)+' '+vasp_ex[0]
	vtst_cmd = 'mpirun -n '+str(nprocs)+' '+vtst_ex[0]
	gamvasp_cmd = 'mpirun -n '+str(nprocs)+' '+vasp_ex[1]
	gamvtst_cmd = 'mpirun -n '+str(nprocs)+' '+vtst_ex[1]
	vasp_module = settings.vasp_module
	vtst_module = settings.vtst_module
	if 'vasp' in fix_version:
		current_module = vasp_module
	elif 'vtst' in fix_version:
		current_module = vtst_module
	elif 'vtst' in calc_swaps:
		current_module = vtst_module
	else:
		current_module = vasp_module
	if os.path.isfile('run_vasp.py'):
		with open('run_vasp.py','r') as rf:
			for line in rf:
				line = line.strip()
				if 'module load '+vtst_module in line:
					old_module = vtst_module
					break
				if 'module load '+vasp_module in line:
					old_module = vasp_module
					break
	else:
		old_module = current_module
	if old_module != current_module:
		module_cmd = 'module unload '+old_module+'; module load '+current_module+' && '
	else:
		module_cmd = 'module load '+current_module+' && '
	runvasp_file = open('run_vasp.py','w')
	if current_module == vtst_module:
		if gpt_version == True:
			runvasp_file.write("import os\nexitcode = os.system("+"'"+module_cmd+gamvtst_cmd+"'"+')')
		else:
			runvasp_file.write("import os\nexitcode = os.system("+"'"+module_cmd+vtst_cmd+"'"+')')
	else:
		if gpt_version == True:
			runvasp_file.write("import os\nexitcode = os.system("+"'"+module_cmd+gamvasp_cmd+"'"+')')
		else:
			runvasp_file.write("import os\nexitcode = os.system("+"'"+module_cmd+vasp_cmd+"'"+')')
	runvasp_file.close()
