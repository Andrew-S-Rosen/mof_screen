import os
from settings import (submit_script, vasp_module, vtst_module, vasp_ex,
	vtst_ex, parallel_cmd)

def get_nprocs():
#Get the number of CPUs from the job submit script

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

	base = parallel_cmd+' '+str(nprocs)+' '
	vasp_cmd = base+vasp_ex[0]
	vtst_cmd = base+vtst_ex[0]
	gamvasp_cmd = base+vasp_ex[1]
	gamvtst_cmd = base+vtst_ex[1]
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
		module_cmd = ('module unload '+old_module+'; module load '+
			current_module+' && ')
	else:
		module_cmd = 'module load '+current_module+' && '
	runvasp_file = open('run_vasp.py','w')
	if current_module == vtst_module:
		if gpt_version == True:
			runvasp_file.write("import os\nexitcode = os.system("
				+"'"+module_cmd+gamvtst_cmd+"'"+')')
		else:
			runvasp_file.write("import os\nexitcode = os.system("
				+"'"+module_cmd+vtst_cmd+"'"+')')
	else:
		if gpt_version == True:
			runvasp_file.write("import os\nexitcode = os.system("
				+"'"+module_cmd+gamvasp_cmd+"'"+')')
		else:
			runvasp_file.write("import os\nexitcode = os.system("
				+"'"+module_cmd+vasp_cmd+"'"+')')

	runvasp_file.close()
