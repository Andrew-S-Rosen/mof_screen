from settings import stdout_file
from ase.io import read
from magmom_handler import get_incar_magmoms, continue_failed_magmoms
from calc_swaps import update_calc

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
	errormsg = list(set(errormsg))
	return errormsg

def get_warning_msgs(outcarfile):
#read in any warning messages
	warningmsg = []
	with open(outcarfile,'r') as rf:
		for line in rf:
			if 'You have a (more or less)' in line:
				warningmsg.append('large_supercell')
	warningmsg = list(set(warningmsg))
	return warningmsg

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

def reset_mof():
#reset the MOF object to the POSCAR/INCAR
	mof = read('POSCAR')
	mof.set_initial_magnetic_moments(get_incar_magmoms('INCAR','POSCAR'))	
	return mof

def continue_mof():
#decide if the MOF object should be continued or restarted
	try:
		mof = read('CONTCAR')
		mof = continue_failed_magmoms(mof)
	except:
		mof = reset_mof()
	return mof