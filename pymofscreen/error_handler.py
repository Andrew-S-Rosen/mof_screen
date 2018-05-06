from ase.io import read
from pymofscreen.magmom_handler import get_incar_magmoms, continue_failed_magmoms
from pymofscreen.calc_swaps import update_calc

def get_error_msgs(outcarfile,refcode,stdout_file):
	"""
	Parse error messages from VASP
	Args:
		outcarfile (string): parth to OUTCAR file
		refcode (string): name of MOF
		stdout_file (string): path to stdout file
	Returns:
		errormsg (list of strings): error messages in OUTCAR/stdout
	"""

	errormsg = []
	start = False
	with open(outcarfile,'r') as rf:
		for line in rf:
			errormsg = check_line_for_error(line,errormsg)
	with open(stdout_file,'r') as rf:
		for line in rf:
			if 'STARTING '+refcode in line:
				start = True
			if start:
				errormsg = check_line_for_error(line,errormsg)
	errormsg = list(set(errormsg))

	return errormsg

def get_warning_msgs(outcarfile):
	"""
	Parse warning messages from VASP
	Args:
		outcarfile (string): parth to OUTCAR file
	Returns:
		warningmsg (list of strings): warning messages in OUTCAR
	"""

	warningmsg = []
	with open(outcarfile,'r') as rf:
		for line in rf:
			if 'You have a (more or less)' in line:
				warningmsg.append('large_supercell')
	warningmsg = list(set(warningmsg))

	return warningmsg

def check_line_for_error(line,errormsg):
	"""
	Parse given line for VASP error code
	Args:
		line (string): error statement
	Returns:
		errormsg (string): VASP error code
	"""

	if 'inverse of rotation matrix was not found (increase SYMPREC)' in line:
		errormsg.append('inv_rot_mat')
	if 'WARNING: Sub-Space-Matrix is not hermitian in DAV' in line:
		errormsg.append('subspacematrix')
	if 'Routine TETIRR needs special values' in line:
		errormsg.append('tetirr')
	if 'Could not get correct shifts' in line:
		errormsg.append('incorrect_shift')
	if ('REAL_OPTLAY: internal error' in line
		or 'REAL_OPT: internal ERROR' in line):
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
	if ('ZBRENT: fatal internal in' in line
		or 'ZBRENT: fatal error in bracketing' in line):
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

def update_calc_after_errors(calc,calc_swaps,errormsg):
	"""
	Update an ASE Vasp calculators object based on error messages
	Args:
		calc (dict): ASE Vasp calculator dictionary
		calc_swaps (list of strings): list of calc swaps
		errormsg (list of strings): list of error messages
	Returns:
		calc (dict): updated ASE Vasp calculator
		errormsg (list of strings): list of error messages
	"""

	for msg in errormsg:
		if msg not in calc_swaps:
			calc_swaps.append(msg)

	calc, calc_swaps = update_calc(calc,calc_swaps)

	for swap in calc_swaps:

		if swap == 'too_few_bands':
			with open('OUTCAR','r') as outcarfile:
				for line in outcarfile:
					if 'NBANDS' in line:
						try:
							d = line.split("=")
							nbands = int(d[-1].strip())
							break
						except (IndexError, ValueError):
							pass
			nbands = int(1.1*nbands)
			calc_swaps.append('nbands='+nbands)
			calc.int_params['nbands'] = nbands
			calc_swaps.remove('too_few_bands')

		elif swap == 'brions':
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
			calc_swaps.remove('brions')

	return calc, calc_swaps

def reset_mof():
	"""
	Reset the ASE atoms object to the POSCAR/INCAR settings
	Returns:
		mof (ASE Atoms object): reset ASE Atoms object
	"""

	mof = read('POSCAR')
	mof.set_initial_magnetic_moments(get_incar_magmoms('INCAR','POSCAR'))	

	return mof

def continue_mof():
	"""
	Update ASE Atoms object after failed job
	Returns:
		mof (ASE Atoms object): reset ASE Atoms object
	"""

	try:
		mof = read('CONTCAR')
		mof = continue_failed_magmoms(mof)
	except:
		mof = reset_mof()

	return mof

def get_niter(outcarfile):
	"""
	Get the number of ionic steps that were run
	Args:
		outcarfile (string): full path to OUTCAR file
	Returns:
		niter (int): number of ionic iterations
	"""

	with open(outcarfile,'r') as rf:
		for line in rf:
			if '- Iteration' in line:
				niter = line.split('(')[0].split('n')[-1].strip()
	niter = int(niter)
	return niter