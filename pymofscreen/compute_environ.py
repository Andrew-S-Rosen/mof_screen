def get_nprocs(submit_script):
	"""
	Get the number of processors from submit script
	Args:
		submit_script (string): name of submission script
	Returns:
		nprocs (int): number of total processors
		ppn (int): number of processors per node
	"""

	#Setup for MOAB
	with open(submit_script,'r') as rf:
		for line in rf:
			if 'nodes' in line or 'ppn' in line:
				line = line.strip().replace(' ','')
				nodes = int(line.split('nodes=')[1].split(':ppn=')[0])
				ppn = int(line.split('nodes=')[1].split(':ppn=')[1])
	nprocs = nodes*ppn

	#Setup for NERSC
	# with open(submit_script,'r') as rf:
	# 	for line in rf:
	# 		if 'SBATCH -N' in line:
	# 			nodes = int(line.split('-N ')[1])
	# ppn = 32
	# nprocs = nodes*ppn

	#Setup for Stampede2
	# with open(submit_script,'r') as rf:
	# 	for line in rf:
	# 		if '-N' in line:
	# 			line = line.strip().replace(' ','')
	# 			nodes = int(line.split('-N')[1])
	# 		if '--ntasks-per-node' in line:
	# 			line = line.strip().replace(' ','')
	# 			ppn = int(line.split('=')[1])
	# nprocs = nodes*ppn

	return nprocs, ppn

def choose_vasp_version(gpt_version,nprocs):
	"""
	Choose the appropriate VASP version (std or gam)
	Args:
		gpt_version (bool): True if gamma-point only or False
		if standard version
		nprocs (int): total number of processors
	"""

	runvasp_file = open('run_vasp.py','w')

	#Setup for A.S. Rosen on Quest
	parallel_cmd = 'mpirun -n'
	vasp_path = '/home/asr731/software/vasp_builds/bin/'
	vasp_ex = [vasp_path+'vasp_std',vasp_path+'vasp_gam']
	module_cmd = 'module load mpi/openmpi-1.8.3-intel2013.2'

	#Setup for NERSC
	# parallel_cmd = 'srun -n'
	# vasp_path = ''
	# vasp_ex = [vasp_path+'vasp_std',vasp_path+'vasp_gam']
	# module_cmd = 'module load vasp-tpc/5.4.1'

	#Setup for Stampede2
	# parallel_cmd = 'ibrun -n'
	# vasp_path = ''
	# vasp_ex = [vasp_path+'vasp_std_vtst',vasp_path+'vasp_gam_vtst']
	# module_cmd = 'module load vasp/5.4.4'

	#Setting up run_vasp.py
	base = parallel_cmd+' '+str(nprocs)+' '
	vasp_cmd = base+vasp_ex[0]
	gamvasp_cmd = base+vasp_ex[1]
	if gpt_version == True:
		runvasp_file.write("import os\nexitcode = os.system("
			+"'"+module_cmd+' && '+gamvasp_cmd+"'"+')')
	else:
		runvasp_file.write("import os\nexitcode = os.system("
			+"'"+module_cmd+' && '+vasp_cmd+"'"+')')

	runvasp_file.close()
