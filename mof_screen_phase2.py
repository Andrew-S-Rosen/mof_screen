import os
from shutil import copyfile
import numpy as np
from ase.calculators.vasp import Vasp
from ase.optimize import BFGSLineSearch
from ase.io import read

#-------------SET PATHS-------------
old_mofpath = '/projects/p30148/vasp_jobs/MOFs/reoptimized_core1/results/'
mofpath = '/projects/p30148/vasp_jobs/MOFs/reoptimized_core1/results_cifs/oxygenated_reoptimized_oms_cifs/'
basepath = '/projects/p30148/vasp_jobs/MOFs/oxidized_oms/'
submit_script = 'sub_asevasp_screening2_temp.job'
stdout_file = 'mof_screen_phase2.out'
ads_species = 'O'
skip_mofs = []

#-------------DEFAULT PARAMETERS-------------
defaults = {
	'xc': 'PBE',
	'ivdw': 12,
	'encut': 520,
	'prec': 'Accurate',
	'algo': 'All',
	'ediff': 1e-4,
	'nelm': 250,
	'nelmin': 6,
	'lreal': False,
	'ncore': 24,
	'ismear': 0,
	'sigma': 0.01,
	'nsw': 500,
	'isif': 2,
	'ediffg': -0.03,
	'lorbit': 11,
	'kppa_lo': 100,
	'kppa_hi': 1000,
	'isym': 0
	}

#Define blocks for elements for spin polarization assignments
#Note: Zn, Cd, Hg treated as non-TMs
sblock_metals = [3,4,11,12,19,20,37,38,55,56,87,88]
dblock_metals = np.concatenate((np.arange(21,30,1),np.arange(39,48,1),np.arange(71,80,1),np.arange(103,112,1)),axis=0).tolist()
fblock_metals = np.concatenate((np.arange(57,71,1),np.arange(89,103,1)),axis=0).tolist()
nonmetal_list = [1,2,6,7,8,9,10,15,16,17,18,34,35,36,53,54,86]
metal_list = [val for val in np.arange(1,119,1) if val not in nonmetal_list]
mag_list = [metal for metal in metal_list if metal not in sblock_metals]
nomag_list = [val for val in np.arange(1,119,1) if val not in mag_list]
poor_metals = [metal for metal in metal_list if metal not in dblock_metals+fblock_metals+sblock_metals]

#-------------FUNCTION DECLARATION-------------
def get_nprocs():
#Get the number of CPUs (modify for each computing environment)
	with open(submit_script,'r') as rf:
		for line in rf:
			if 'nodes' in line or 'ppn' in line:
				line = line.strip().replace(' ','')
				nodes = int(line.split('nodes=')[1].split(':ppn=')[0])
				ppn = int(line.split('nodes=')[1].split(':ppn=')[1])
	nprocs = nodes*ppn
	return nprocs, ppn

def get_kpts(cif_file,kppa):
#Get kpoint grid at a given KPPA
	cif_split1 = cif_file.split('_'+ads_species+'_OMS')[0]
	old_cif_name = cif_split1.split('_spin')[0]
	spin = 'spin'+cif_split1.split('_spin')[1]
	if kppa == 100:
		kpts_path = old_mofpath+old_cif_name+'/isif2/'+spin+'/KPOINTS'
	elif kppa == 1000:
		kpts_path = old_mofpath+old_cif_name+'/final/'+spin+'/KPOINTS'
	else:
		raise ValueError('Incompatible KPPA with prior runs')
	with open(kpts_path,'r') as rf:
		for i, line in enumerate(rf):
			line = line.strip()
			if i == 2:
				if 'gamma' in line.lower():
					gamma = True
				else:
					gamma = False
			if i == 3:
				kpts = np.squeeze(np.asarray(np.matrix(line))).tolist()
	if len(kpts) != 3:
		raise ValueError('Error parsing KPOINTS file')
	return kpts, gamma

def choose_vasp_version(kpts,n_atoms,nprocs,ppn):
#Run the gamma pt only or regular VASP version
	if sum(kpts) == 3:
		gamma_version = True
	else:
		gamma_version = False
	while n_atoms < nprocs/2:
		nprocs = nprocs - ppn
		if nprocs == ppn:
			break
	gamvasp_cmd = 'mpirun -n '+str(nprocs)+' vasp_gam'
	vasp_cmd =  'mpirun -n '+str(nprocs)+' vasp_std'
	runvasp_file = open('run_vasp.py','w')
	if gamma_version == True:
		runvasp_file.write("import os\nexitcode = os.system("+"'"+gamvasp_cmd+"'"+')')
	else:
		runvasp_file.write("import os\nexitcode = os.system("+"'"+vasp_cmd+"'"+')')
	runvasp_file.close()

def pprint(printstr):
#Redirect print commands to log file
	print(printstr)
	with open('screening.log','a') as txtfile:
		txtfile.write(printstr+'\n')

def get_cif_files():
#Get CIF files from mofpath
	cif_files = []
	for filename in os.listdir(mofpath):
		filename = filename.strip()
		if len(filename.split('.cif')) == 2:
			refcode = filename.split('.cif')[0]
			if refcode not in skip_mofs:
				cif_files.append(filename)
			else:
				pprint('Skipping '+refcode)
	sorted_cifs = sorted(cif_files)
	return sorted_cifs

def cif_to_mof(cif_file):
#Save MOF file as reduced unit cell and read in ASE
	cifpath = mofpath+cif_file
	mof = read(cifpath)
	return mof

def prep_paths():
#Folder and file cleanup
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
	if os.path.isfile(screen_results_path) == True:
		open(screen_results_path, 'w').close()
	if os.path.isfile(log_file) == True:
		open(log_file, 'w').close()

def clean_files(remove_files):
#clean files
	for file in remove_files:
		if os.path.isfile(file) == True:
			os.remove(file)

def get_mag_indices(mof):
#Get indices of d-block, f-block, and semimetal atoms
	mag_indices = []
	for i, atom in enumerate(mof):
		if atom.number in mag_list:
			mag_indices.append(i)
	return mag_indices

def set_initial_magmoms(mof,spin_level):
#Add initial magnetic moments to atoms object
	mag_indices = get_mag_indices(mof)
	mof.set_initial_magnetic_moments(np.zeros(len(mof)))
	for i, atom in enumerate(mof):
		if i in mag_indices:
			mag_number = atom.number
			if spin_level == 'spin1':
				if mag_number in dblock_metals:
					atom.magmom = 5.0
				elif mag_number in fblock_metals:
					atom.magmom = 7.0
				elif mag_number in poor_metals:
					atom.magmom = 0.1
				else:
					raise ValueError('Metal not properly classified')
			elif spin_level == 'spin2':
				atom.magmom = 0.1
			else:
				raise ValueError('Spin iteration out of range')
	return mof

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

def continue_magmoms(mof,incarpath):
#Read in the old magmoms
	mag_indices = get_mag_indices(mof)
	ispin = False
	with open(incarpath,'r') as incarfile:
		for line in incarfile:
			line = line.strip()
			if 'ISPIN = 2' in line:
				ispin = True
				mof_magmoms = mof.get_magnetic_moments()
				mof.set_initial_magnetic_moments(mof_magmoms)
				abs_magmoms = np.abs(mof_magmoms[mag_indices])
	if ispin == False:
		abs_magmoms = np.zeros(len(mag_indices))
	return mof, abs_magmoms

def continue_failed_magmoms(mof):
	self_resort = []
	file = open('ase-sort.dat', 'r')
	lines = file.readlines()
	file.close()
	for line in lines:
	    data = line.split()
	    self_resort.append(int(data[1]))
	magnetic_moments = np.zeros(len(mof))
	n = 0
	lines = open('OUTCAR', 'r').readlines()
	for line in lines:
	    if line.rfind('magnetization (x)') > -1:
	        for m in range(len(mof)):
	            magnetic_moments[m] = float(lines[n + m + 4].split()[4])
	    n += 1
	sorted_magmoms = np.array(magnetic_moments)[self_resort]
	ispin = False
	with open('INCAR','r') as incarfile:
		for line in incarfile:
			line = line.strip()
			if 'ISPIN = 2' in line:
				ispin = True
	if ispin == True and all(sorted_magmoms == 0.0) == True:
		raise ValueError('Error reading magmoms from failed OUTCAR')
	mof.set_initial_magnetic_moments(sorted_magmoms)
	return mof

def get_incar_magmoms(incarpath,poscarpath):
#get the magnetic moments from the POSCAR
	mof_mag_list = []
	init_mof = read(poscarpath)
	with open(incarpath,'r') as incarfile:
		for line in incarfile:
			line = line.strip()
			if 'MAGMOM' in line:
				mag_line = line.split('= ')[1:][0].split(' ')
				for val in mag_line:
					mag = float(val.split('*')[1])
					num = int(val.split('*')[0])
					mof_mag_list.extend([mag]*num)
	if bool(mof_mag_list) == False:
		mof_mag_list = np.zeros(len(init_mof))
	if len(mof_mag_list) != len(mof_mag_list):
		raise ValueError('Error reading INCAR magnetic moments')
	return mof_mag_list

def read_outcar(outcarpath):
#read OUTCAR and fixes weird fortran I/O errors
	try:
		mof = read(outcarpath,format='vasp-out')
	except ValueError:
		os.system('sed -i -e "s/\([[:digit:]]\)-\([[:digit:]]\)/\\1 -\\2/g" '+str(outcarpath))
		mof = read(outcarpath,format='vasp-out')
	return mof

def reset_mof():
#reset the MOF object to the POSCAR/INCAR
	mof = read('POSCAR')
	mof.set_initial_magnetic_moments(get_incar_magmoms('INCAR','POSCAR'))	
	return mof

def continue_mof():
#decide if the MOF object should be continued or restarted
	try:
		mof = read('CONTCAR')
		try:
			mof, abs_magmoms = continue_magmoms(mof,'INCAR')
		except:
			mof = continue_failed_magmoms(mof)
	except:
		mof = reset_mof()
	return mof

def check_line_for_error(line,errormsg):
	if 'inverse of rotation matrix was not found (increase SYMPREC)' in line:
		errormsg.append('inv_rot_mat')
	if 'WARNING: Sub-Space-Matrix is not hermitian in DAV' in line:
		errormsg.append('subspacematrix')
	if 'Routine TETIRR needs special values' in line:
		errormsg.append('tetirr')
	if 'Could not get correct shifts' in line:
		errormsg.append('incorrect_shift')
	if 'REAL_OPTLAY: internal error' in line or 'REAL_OPT: internal ERROR' in line:
		errormsg.append('real_optlay')
	if 'ERROR RSPHER' in line:
		errormsg.append('rspher')
	if 'DENTET' in line:
		errormsg.append('dentet')
	if 'TOO FEW BANDS' in line:
		errormsg.append('too_few_bands')
	if 'Found some non-integer element in rotation matrix' in line:
		errormsg.append('rot_matrix')
	if 'BRIONS problems: POTIM should be increased' in line:
		errormsg.append('brions')
	if 'internal error in subroutine PRICEL' in line:
		errormsg.append('pricel')
	if 'One of the lattice vectors is very long (>50 A), but AMIN' in line:
		errormsg.append('amin')
	if 'ZBRENT: fatal internal in' in line or 'ZBRENT: fatal error in bracketing' in line:
		errormsg.append('zbrent')
	if 'ERROR in subspace rotation PSSYEVX' in line:
		errormsg.append('pssyevx')
	if 'WARNING in EDDRMM: call to ZHEGV failed' in line:
		errormsg.append('eddrmm')
	if 'Error EDDDAV: Call to ZHEGV failed' in line:
		errormsg.append('edddav')
	if 'EDWAV: internal error, the gradient is not orthogonal' in line:
		errormsg.append('grad_not_orth')
	if 'ERROR: SBESSELITER : nicht konvergent' in line:
		errormsg.append('nicht_konv')
	if 'ERROR EDDIAG: Call to routine ZHEEV failed!' in line:
		errormsg.append('zheev')
	if 'ELF: KPAR>1 not implemented' in line:
		errormsg.append('elf_kpar')
	if 'RHOSYG internal error' in line:
		errormsg.append('rhosyg')
	if 'POSMAP internal error: symmetry equivalent atom not found' in line:
		errormsg.append('posmap')
	return errormsg

def get_error_msgs(outcarfile,refcode):
#read in any error messages
	errormsg = []
	start = False
	with open(outcarfile,'r') as rf:
		for line in rf:
			errormsg = check_line_for_error(line,errormsg)
	with open(stdout_file,'r') as rf:
		for line in rf:
			if 'STARTING '+refcode in line:
				start = True
			if start == True:
				errormsg = check_line_for_error(line,errormsg)
	return errormsg

def get_warning_msgs(outcarfile):
#read in any warning messages
	warningmsg = []
	with open(outcarfile,'r') as rf:
		for line in rf:
			if 'You have a (more or less)' in line:
				warningmsg.append('large_supercell')
	return warningmsg

def update_calc(calc,calc_swaps):
#update calculator based on calc swaps
	for swap in calc_swaps:
		swap.replace(' ','')
		if swap == 'large_supercell':
			calc.special_params['lreal'] = 'Auto'
		elif 'sigma=' in swap:
			calc.float_params['sigma'] = float(swap.split('=')[-1])
		elif 'nbands=' in swap:
			calc.int_params['nbands'] = int(swap.split('=')[-1])
		elif 'potim=' in swap:
			calc.float_params['potim'] = float(swap.split('=')[-1])
		elif 'nsw=' in swap:
			calc.int_params['nsw'] = int(swap.split('=')[-1])
		elif 'nelm=' in swap:
			calc.int_params['nelm'] = int(swap.split('=')[-1])
		elif 'ibrion=' in swap:
			calc.int_params['ibrion'] = int(swap.split('=')[-1])
		elif 'istart=' in swap:
			calc.int_params['istart'] = int(swap.split('=')[-1])
		elif 'algo=' in swap:
			calc.string_params['algo'] = swap.split('=')[-1]
		elif 'isif=' in swap:
			calc.int_params['isif'] = int(swap.split('=')[-1])
		elif swap == 'edddav':
			calc.string_params['algo'] = 'All'
		elif swap == 'inv_rot_mat':
			calc.exp_params['symprec'] = 1e-8
		elif swap == 'subspacematrix' or swap == 'real_optlay' or swap == 'rspher' or swap == 'nicht_konv':
			calc.special_params['lreal'] = False
			calc.string_params['prec'] = 'Accurate'
		elif swap == 'tetirr' or swap == 'incorrect_shift':
			calc.input_params['gamma'] = True
		elif swap == 'dentet' or swap == 'grad_not_orth':
			calc.int_params['ismear'] = 0
		elif swap == 'rot_matrix':
			calc.input_params['gamma'] = True
			calc.int_params['isym'] = 0
		elif swap == 'pricel':
			calc.exp_params['symprec'] = 1e-8
			calc.int_params['isym'] = 0
		elif swap == 'amin':
			calc.float_params['amin'] = 0.01
		elif swap == 'zbrent':
			calc.int_params['ibrion'] = 1
			calc.exp_params['ediff'] = 1e-6
			calc.int_params['nelmin'] = 8
		elif swap == 'pssyevx' or swap == 'eddrmm':
			calc.string_params['algo'] = 'Normal'
		elif swap == 'zheev':
			calc.string_params['algo'] = 'Exact'
		elif swap == 'elf_kpar':
			calc.int_params['kpar'] = 1
		elif swap == 'rhosyg':
			calc.exp_params['symprec'] = 1e-4
			calc.int_params['isym'] = 0
		elif swap == 'posmap':
			calc.exp_params['symprec'] = 1e-6
	return calc

def update_calc_after_errors(calc,calc_swaps,errormsg):
#make a calc swap based on errors
	for msg in errormsg:
		if msg not in calc_swaps:
			calc_swaps.append(msg)
	calc = update_calc(calc,calc_swaps)
	for swap in calc_swaps:
		if msg == 'too_few_bands':
			with open('OUTCAR','r') as outcarfile:
				for line in outcarfile:
					if "NBANDS" in line:
						try:
							d = line.split("=")
							nbands = int(d[-1].strip())
							break
						except (IndexError, ValueError):
							pass
			nbands = int(1.1*nbands)
			calc_swaps.append('nbands='+nbands)
			calc.int_params['nbands'] = nbands
		elif msg == 'brions':
			with open('OUTCAR','r') as outcarfile:
				for line in outcarfile:
					if 'POTIM' in line:
						try:
							potim = float(d.split('=')[-1].split('time-step')[0])
							break
						except (IndexError, ValueError):
							pass
			calc_swaps.append('potim='+potim)
			calc.float_params['potim'] = potim
	return calc, calc_swaps

def get_niter(outcarfile):
	with open(outcarfile,'r') as rf:
		for line in rf:
			if '- Iteration' in line:
				niter = line.split('(')[0].split('n')[-1].strip()
	niter = int(niter)
	return niter

def mof_run(mof,calc,cif_file,calc_swaps):
#Get the optimized structure of the MOF
	copyfile(mofpath+cif_file,basepath+'working/'+cif_file)
	success = False
	calc = update_calc(calc,calc_swaps)
	mof.set_calculator(calc)
	try:
		mof.get_potential_energy()
		niter = get_niter('OUTCAR')
		if niter < mof.calc.int_params['nsw'] and mof.calc.converged != True:
			raise SystemError('VASP stopped but did not crash and burn')
		success = True
	except:
		pprint('Original run failed. Trying to auto-solve issue.')
		old_error_len = 0
		refcode = cif_file.split('.cif')[0]
		if os.path.isfile('WAVECAR'):
			os.remove('WAVECAR')
		while True:
			errormsg = get_error_msgs('OUTCAR',refcode)
			calc, calc_swaps = update_calc_after_errors(calc,calc_swaps,errormsg)
			error_len = len(errormsg)
			if error_len == old_error_len:
				break
			print('VASP failed with error(s). Trying to auto-solve issue.',errormsg)
			mof = continue_mof()
			try:				
				mof.set_calculator(calc)
				mof.get_potential_energy()
				niter = get_niter('OUTCAR')
				if niter < mof.calc.int_params['nsw'] and mof.calc.converged != True:
					raise SystemError('VASP stopped but did not crash and burn')
				success = True
			except:
				pass
			old_error_len = error_len
	if success == False:
		mof = None
	return mof, calc_swaps

def mof_bfgs_run(mof,calc,cif_file,calc_swaps,steps,fmax):
#Optimize with BFGSLineSearch
	copyfile(mofpath+cif_file,basepath+'working/'+cif_file)
	success = False
	calc = update_calc(calc,calc_swaps)
	mof.set_calculator(calc)
	try:
		dyn = BFGSLineSearch(mof,trajectory='opt.traj')
		dyn.run(fmax=fmax,steps=steps)
		success = True
	except:
		pprint('Original run failed. Trying to auto-solve issue.')
		old_error_len = 0
		refcode = cif_file.split('.cif')[0]
		if os.path.isfile('WAVECAR'):
			os.remove('WAVECAR')
		while True:
			errormsg = get_error_msgs('OUTCAR',refcode)
			calc, calc_swaps = update_calc_after_errors(calc,calc_swaps,errormsg)
			error_len = len(errormsg)
			if error_len == old_error_len:
				break
			print('VASP failed with error(s). Trying to auto-solve issue.',errormsg)
			mof = continue_mof()
			try:				
				mof.set_calculator(calc)
				dyn = BFGSLineSearch(mof,trajectory='opt.traj')
				dyn.run(fmax=fmax,steps=steps)
				success = True
			except:
				pass
			old_error_len = error_len
	if success == False:
		mof = None
	return mof, dyn, calc_swaps

def prep_next_run(acc_level,run_i,refcode,spin_level):
#Update counter and decide if next job should be skipped
	skip_spin2 = False
	success_path = basepath+'results/'+refcode+'/'+acc_level+'/'+spin_level
	incarpath = success_path+'/INCAR'
	outcarpath = success_path+'/OUTCAR'
	errorpath = basepath+'errors/'+refcode+'/'+acc_level+'/'+spin_level
	if os.path.exists(errorpath) == True:
		mof = None
	else:
		mof = read_outcar(outcarpath)
		if acc_level != 'scf_test':
			mof, abs_magmoms = continue_magmoms(mof,incarpath)
			mag_indices = get_mag_indices(mof)
			mag_nums = mof[mag_indices].get_atomic_numbers()
			if np.sum(abs_magmoms < 0.1) == len(abs_magmoms) or all(num in sblock_metals+poor_metals for num in mag_nums) == True:
				skip_spin2 = True
	run_i += 1
	return mof, run_i, skip_spin2

def manage_restart_files(file_path):
#Make sure WAVECAR is copied
	files = ['WAVECAR']
	for file in files:
		if os.path.isfile(file) != True or os.stat(file).st_size == 0:
			if os.path.isfile(file_path+'/'+file) == True:
				copyfile(file_path+'/'+file,os.getcwd()+'/'+file)
	return

def write_energy(refcode,acc_level,spin_level):
#Write energy to results file
	outcarpath = basepath+'results/'+refcode+'/'+acc_level+'/'+spin_level+'/OUTCAR'
	final_mof = read_outcar(outcarpath)
	E = final_mof.get_potential_energy()
	str_to_write = refcode+'_'+spin_level+': '+str(E)
	results_file = basepath+'results/screen_results.dat'
	with open(results_file,'a') as txtfile:
		txtfile.write(str_to_write+'\n')

def calcs(run_i):
#Get calculator
	if run_i == 0:
		calc = Vasp(
			xc=defaults['xc'],
			kpts=defaults['kpts_lo'],
			gamma=defaults['gamma'],
			ivdw=defaults['ivdw'],
			prec=defaults['prec'],
			algo=defaults['algo'],
			nelm=1,
			lreal=defaults['lreal'],
			ncore=defaults['ncore'],
			ismear=defaults['ismear'],
			sigma=defaults['sigma'],
			istart=0,
			lcharg=False,
			lwave=False,
			isym=defaults['isym']
			)
	elif run_i == 1:
		calc = Vasp(
			xc=defaults['xc'],
			kpts=defaults['kpts_lo'],
			gamma=defaults['gamma'],
			ivdw=defaults['ivdw'],
			prec=defaults['prec'],
			algo=defaults['algo'],
			ediff=defaults['ediff'],
			nelm=defaults['nelm'],
			nelmin=defaults['nelmin'],
			lreal=defaults['lreal'],
			ncore=defaults['ncore'],
			ismear=defaults['ismear'],
			sigma=defaults['sigma'],
			lcharg=False,
			lwave=True,
			lorbit=defaults['lorbit'],
			isym=defaults['isym']
			)
	elif run_i == 1.5:
		calc = Vasp(
			xc=defaults['xc'],
			kpts=defaults['kpts_lo'],
			gamma=defaults['gamma'],
			ivdw=defaults['ivdw'],
			prec=defaults['prec'],
			algo=defaults['algo'],
			ediff=defaults['ediff'],
			nelm=defaults['nelm'],
			nelmin=defaults['nelmin'],
			lreal=defaults['lreal'],
			ncore=defaults['ncore'],
			ismear=defaults['ismear'],
			sigma=defaults['sigma'],
			lcharg=False,
			lwave=True,
			ibrion=2,
			isif=defaults['isif'],
			nsw=200,
			ediffg=-0.05,
			lorbit=defaults['lorbit'],
			isym=defaults['isym']
			)
	elif run_i == 2:
		calc = Vasp(
			xc=defaults['xc'],
			kpts=defaults['kpts_hi'],
			gamma=defaults['gamma'],
			ivdw=defaults['ivdw'],
			prec=defaults['prec'],
			algo=defaults['algo'],
			ediff=defaults['ediff'],
			nelm=defaults['nelm'],
			nelmin=defaults['nelmin'],
			lreal=defaults['lreal'],
			ncore=defaults['ncore'],
			ismear=defaults['ismear'],
			sigma=defaults['sigma'],
			lcharg=False,
			lwave=True,
			ibrion=2,
			isif=defaults['isif'],
			nsw=defaults['nsw'],
			ediffg=defaults['ediffg'],
			lorbit=defaults['lorbit'],
			isym=defaults['isym']
			)
	elif run_i == 3:
		calc = Vasp(
			xc=defaults['xc'],
			encut=defaults['encut'],
			kpts=defaults['kpts_hi'],
			gamma=defaults['gamma'],
			ivdw=defaults['ivdw'],
			prec=defaults['prec'],
			algo=defaults['algo'],
			ediff=defaults['ediff'],
			nelm=defaults['nelm'],
			nelmin=defaults['nelmin'],
			lreal=defaults['lreal'],
			ncore=defaults['ncore'],
			ismear=defaults['ismear'],
			sigma=defaults['sigma'],
			lcharg=True,
			laechg=True,
			lwave=True,
			ibrion=2,
			isif=defaults['isif'],
			nsw=defaults['nsw'],
			ediffg=defaults['ediffg'],
			lorbit=defaults['lorbit'],
			isym=defaults['isym']
			)
	else:
		raise ValueError('Out of range for calculators')
	if 'Li' in calc.setups_defaults:
		del calc.setups_defaults['Li']
	return calc

def run_screen(cif_files):
#Run high-throughput screening

	#Files, spin levels, and accuracy levels to iterate over
	vasp_files = ['INCAR','POSCAR','KPOINTS','POTCAR','OUTCAR',
	'CONTCAR','CHGCAR','AECCAR0','AECCAR2','WAVECAR','opt.traj',
	'vasprun.xml']
	spin_levels = ['spin1','spin2']
	acc_levels = ['scf_test','isif2_lowacc','isif2_medacc','final']
	nprocs, ppn = get_nprocs()

	#for each CIF file, optimize the structure
	for cif_file in cif_files:

		refcode = cif_file.split('.cif')[0]
		pprint('***STARTING '+refcode+'***')

		#Make sure MOF isn't running on other process
		working_cif_path = basepath+'working/'+cif_file
		if os.path.isfile(working_cif_path) == True:
			pprint('SKIPPED: Running on another process')
			continue

		#Partial paths to write the OUTCARs
		results_partial_paths = []
		error_outcar_partial_paths = []
		for acc_level in acc_levels:
			results_partial_paths.append(basepath+'results/'+refcode+'/'+acc_level)
			error_outcar_partial_paths.append(basepath+'errors/'+refcode+'/'+acc_level)
		spin1_final_mof_path = results_partial_paths[-1]+'/'+spin_levels[0]+'/OUTCAR'

		#Get the kpoints
		kpts_lo, gamma = get_kpts(cif_file,defaults['kppa_lo'])
		kpts_hi, gamma = get_kpts(cif_file,defaults['kppa_hi'])
		defaults['gamma'] = gamma
		defaults['kpts_lo'] = kpts_lo
		defaults['kpts_hi'] = kpts_hi

		#for each spin level, optimize the structure
		for spin_level in spin_levels:

			#***********PREP FOR RUN***********
			calc_swaps = []
			outcar_paths = []
			error_outcar_paths = []
			run_i = 0
			clean_files(vasp_files)
			for results_partial_path in results_partial_paths:
				outcar_paths.append(results_partial_path+'/'+spin_level+'/OUTCAR')
			for error_outcar_partial_path in error_outcar_partial_paths:
				error_outcar_paths.append(error_outcar_partial_path+'/'+spin_level+'/OUTCAR')

			#***********SCF TEST************
			acc_level = acc_levels[run_i]
			if os.path.isfile(outcar_paths[run_i]) != True and os.path.isfile(error_outcar_paths[run_i]) != True:
				if os.path.isfile(spin1_final_mof_path):
					mof = read(spin1_final_mof_path)
				else:
					mof = cif_to_mof(cif_file)
				mof = set_initial_magmoms(mof,spin_level)
				choose_vasp_version(kpts_lo,len(mof),nprocs,ppn)
				pprint('Running '+spin_level+', '+acc_level)
				mof, calc_swaps = mof_run(mof,calcs(run_i),cif_file,calc_swaps)
				if mof != None:
					write_success(refcode,spin_level,acc_level,vasp_files,cif_file)
				else:
					pprint('^ VASP crashed')
					write_errors(refcode,spin_level,acc_level,vasp_files,cif_file)
			elif os.path.isfile(outcar_paths[run_i]) == True:
				pprint('COMPLETED: '+spin_level+', '+acc_level)
			mof, run_i, skip_spin2 = prep_next_run(acc_level,run_i,refcode,spin_level)
			if mof == None:
				pprint('Skipping rest because of errors')
				break
			warnings = get_warning_msgs(outcar_paths[run_i-1])
			calc_swaps.extend(warnings)

			#***********ISIF 2 (lowacc)************
			acc_level = acc_levels[run_i]
			if os.path.isfile(outcar_paths[run_i-1]) == True and os.path.isfile(outcar_paths[run_i]) != True and os.path.isfile(error_outcar_paths[run_i]) != True:
				if os.path.isfile(spin1_final_mof_path):
					mof = read(spin1_final_mof_path)
				else:
					mof = cif_to_mof(cif_file)
				mof = set_initial_magmoms(mof,spin_level)
				choose_vasp_version(kpts_lo,len(mof),nprocs,ppn)
				pprint('Running '+spin_level+', '+acc_level)
				steps = 100
				fmax = 5.0
				mof, dyn, calc_swaps = mof_bfgs_run(mof,calcs(run_i),cif_file,calc_swaps,steps,fmax)
				if mof != None and dyn and mof.calc.scf_converged == True:
					loop_i = 0
					converged = False
					while mof != None and loop_i < 5 and converged == False and mof.calc.scf_converged == True:
						mof = read_outcar('OUTCAR')
						mof, abs_magmoms = continue_magmoms(mof,'INCAR')
						mof, calc_swaps = mof_run(mof,calcs(1.5),cif_file,calc_swaps)
						if mof == None:
							break
						converged = mof.calc.converged
						loop_i += 1
				if mof != None and mof.calc.scf_converged == True and mof.calc.converged == True:
					write_success(refcode,spin_level,acc_level,vasp_files,cif_file)
				else:
					write_errors(refcode,spin_level,acc_level,vasp_files,cif_file)
					if mof == None:
						pprint('^ VASP crashed')
					elif mof.calc.scf_converged == False:
						pprint('^ SCF did not converge')
					elif mof.calc.converged == False:
						pprint('^ Convergence not reached')
			elif os.path.isfile(outcar_paths[run_i]) == True:
				pprint('COMPLETED: '+spin_level+', '+acc_level)
			mof, run_i, skip_spin2 = prep_next_run(acc_level,run_i,refcode,spin_level)
			if mof == None:
				pprint('Skipping rest because of errors')
				break
			if spin_level == 'spin2':
				mag_indices = get_mag_indices(mof)
				old_mof = read(spin1_final_mof_path)
				if np.sum(np.abs(mof.get_initial_magnetic_moments()[mag_indices] - old_mof.get_magnetic_moments()[mag_indices]) >= 0.05) == 0:
					pprint('Skipping rest because SPIN2 converged to SPIN1')
					continue

			#***********ISIF 2 (medacc)************
			acc_level = acc_levels[run_i]
			if os.path.isfile(outcar_paths[run_i-1]) == True and os.path.isfile(outcar_paths[run_i]) != True and os.path.isfile(error_outcar_paths[run_i]) != True:
				choose_vasp_version(kpts_hi,len(mof),nprocs,ppn)
				if sum(kpts_lo) == 3 and sum(kpts_hi) > 3:
					files_to_clean = ['WAVECAR']
					clean_files(files_to_clean)
				else:
					manage_restart_files(results_partial_paths[run_i-1]+'/'+spin_level)
				if kpts_lo != kpts_hi:
					if 'zbrent' in calc_swaps:
						calc_swaps.remove('zbrent')
				pprint('Running '+spin_level+', '+acc_level)
				mof,calc_swaps = mof_run(mof,calcs(run_i),cif_file,calc_swaps)
				if mof != None and mof.calc.scf_converged == True and mof.calc.converged == True:
					write_success(refcode,spin_level,acc_level,vasp_files,cif_file)
				else:
					write_errors(refcode,spin_level,acc_level,vasp_files,cif_file)
					if mof == None:
						pprint('^ VASP crashed')
					elif mof.calc.scf_converged == False:
						pprint('^ SCF did not converge')
					elif mof.calc.converged == False:
						pprint('^ Convergence not reached')
			elif os.path.isfile(outcar_paths[run_i]) == True:
				pprint('COMPLETED: '+spin_level+', '+acc_level)
			mof, run_i, skip_spin2 = prep_next_run(acc_level,run_i,refcode,spin_level)
			if mof == None:
				pprint('Skipping rest because of errors')
				break

			#***********ISIF 2 (final)************
			acc_level = acc_levels[run_i]
			if os.path.isfile(outcar_paths[run_i-1]) == True and os.path.isfile(outcar_paths[run_i]) != True and os.path.isfile(error_outcar_paths[run_i]) != True:
				choose_vasp_version(kpts_hi,len(mof),nprocs,ppn)
				manage_restart_files(results_partial_paths[run_i-1]+'/'+spin_level)
				if 'zbrent' in calc_swaps:
					calc_swaps.remove('zbrent')
				pprint('Running '+spin_level+', '+acc_level)
				mof,calc_swaps = mof_run(mof,calcs(run_i),cif_file,calc_swaps)
				if mof != None and mof.calc.scf_converged == True and mof.calc.converged == True:
					if 'large_supercell' in calc_swaps:
						calc_swaps.remove('large_supercell')
						mof = read_outcar('OUTCAR')
						mof, abs_magmoms = continue_magmoms(mof,'INCAR')
						mof, calc_swaps = mof_run(mof,calcs(run_i),cif_file,calc_swaps)
						if mof != None and mof.calc.scf_converged == True and mof.calc.converged == True:
							write_success(refcode,spin_level,acc_level,vasp_files,cif_file)
						else:
							write_errors(refcode,spin_level,acc_level,vasp_files,cif_file)
							if mof == None:
								pprint('^ VASP crashed')
							elif mof.calc.scf_converged == False:
								pprint('^ SCF did not converge')
							elif mof.calc.converged == False:
								pprint('^ Convergence not reached')
					else:
						write_success(refcode,spin_level,acc_level,vasp_files,cif_file)
				else:
					write_errors(refcode,spin_level,acc_level,vasp_files,cif_file)
					if mof == None:
						pprint('^ VASP crashed')
					elif mof.calc.scf_converged == False:
						pprint('^ SCF did not converge')
					elif mof.calc.converged == False:
						pprint('^ Convergence not reached')
			elif os.path.isfile(outcar_paths[run_i]) == True:
				pprint('COMPLETED: '+spin_level+', '+acc_level)
			mof, run_i, skip_spin2 = prep_next_run(acc_level,run_i,refcode,spin_level)
			if mof == None:
				pprint('Skipping rest because of errors')
				break

			#***********SAVE and CONTINUE***********
			if os.path.isfile(outcar_paths[-1]) == True:
				write_energy(refcode,acc_level,spin_level)
			if skip_spin2 == True:
				pprint('Skipping '+spin_levels[1]+' run')
				break

#Run screen
prep_paths()
cif_files = get_cif_files()
run_screen(cif_files)
