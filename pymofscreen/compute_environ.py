def get_nprocs(submit_script):
#Get the number of CPUs from the job submit script

	with open(submit_script,'r') as rf:
		for line in rf:
			if 'nodes' in line or 'ppn' in line:
				line = line.strip().replace(' ','')
				nodes = int(line.split('nodes=')[1].split(':ppn=')[0])
				ppn = int(line.split('nodes=')[1].split(':ppn=')[1])
	nprocs = nodes*ppn

	return nprocs, ppn

def choose_vasp_version(gpt_version,nprocs,calc_swaps):
#Run a given VASP version
	parallel_cmd = 'mpirun -n'
	vasp_path = '/home/asr731/software/vasp_builds/bin/'
	vasp_ex = [vasp_path+'vasp_std',vasp_path+'vasp_gam']
	module_cmd = 'module load mpi/openmpi-1.8.3-intel2013.2'
	base = parallel_cmd+' '+str(nprocs)+' '
	vasp_cmd = base+vasp_ex[0]
	gamvasp_cmd = base+vasp_ex[1]
	runvasp_file = open('run_vasp.py','w')
	if gpt_version == True:
		runvasp_file.write("import os\nexitcode = os.system("
			+"'"+module_cmd+' && '+gamvasp_cmd+"'"+')')
	else:
		runvasp_file.write("import os\nexitcode = os.system("
			+"'"+module_cmd+' && '+vasp_cmd+"'"+')')

	runvasp_file.close()
