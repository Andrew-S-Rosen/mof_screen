from ase.calculators.vasp import Vasp

#default parameters for calculators
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
	'ediffg': -0.03,
	'lorbit': 11,
	'isym': 0,
	'setups':{'base':'recommended','Li':''}
	}

def calcs_ads(run_i):
	"""
	Define the default calculators for ionic relaxations
	Note: it should not include the kpts or gamma keywords!
	Args:
		run_i (int): iteration number
	Returns:
		calc (dict): ASE Vasp calculator dictionary
	"""
	if run_i == 'scf_test':
		calc = Vasp(
			xc=defaults['xc'],
			setups=defaults['setups'],
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
			isym=defaults['isym'],
			nsw=0,
			istart=0
			)
	elif run_i == 'ase_bfgs':
		calc = Vasp(
			xc=defaults['xc'],
			setups=defaults['setups'],
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
	elif run_i == 'isif2_lowacc':
		calc = Vasp(
			xc=defaults['xc'],
			setups=defaults['setups'],
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
			isif=2,
			nsw=200,
			ediffg=-0.05,
			lorbit=defaults['lorbit'],
			isym=defaults['isym']
			)
	elif run_i == 'isif2_medacc':
		calc = Vasp(
			xc=defaults['xc'],
			setups=defaults['setups'],
			ivdw=defaults['ivdw'],
			prec=defaults['prec'],
			algo=defaults['algo'],
			ediff=1e-6,
			nelm=defaults['nelm'],
			nelmin=8,
			lreal=defaults['lreal'],
			ncore=defaults['ncore'],
			ismear=defaults['ismear'],
			sigma=defaults['sigma'],
			lcharg=False,
			lwave=True,
			ibrion=3,
			iopt=7,
			potim=0,
			isif=2,
			nsw=defaults['nsw'],
			ediffg=defaults['ediffg'],
			lorbit=defaults['lorbit'],
			isym=defaults['isym']
			)
	elif run_i == 'isif2_highacc':
		calc = Vasp(
			xc=defaults['xc'],
			setups=defaults['setups'],
			encut=defaults['encut'],
			ivdw=defaults['ivdw'],
			prec=defaults['prec'],
			algo=defaults['algo'],
			ediff=1e-6,
			nelm=defaults['nelm'],
			nelmin=8,
			lreal=defaults['lreal'],
			ncore=defaults['ncore'],
			ismear=defaults['ismear'],
			sigma=defaults['sigma'],
			lcharg=False,
			lwave=True,
			ibrion=3,
			iopt=7,
			potim=0,
			isif=2,
			nsw=defaults['nsw'],
			ediffg=defaults['ediffg'],
			lorbit=defaults['lorbit'],
			isym=defaults['isym']
			)
	elif run_i == 'final_spe':
		calc = Vasp(
			xc=defaults['xc'],
			setups=defaults['setups'],
			encut=defaults['encut'],
			ivdw=defaults['ivdw'],
			prec=defaults['prec'],
			algo=defaults['algo'],
			ediff=1e-6,
			nelm=defaults['nelm'],
			lreal=False,
			ncore=defaults['ncore'],
			ismear=defaults['ismear'],
			sigma=defaults['sigma'],
			lcharg=True,
			laechg=True,
			lwave=True,
			nsw=0,
			lorbit=defaults['lorbit'],
			isym=defaults['isym'],
			addgrid=False
			)
	else:
		raise ValueError('Out of range for calculators')

	return calc

def calcs_vol(run_i):
#calculator definitions for each run
	if run_i == 0:
		calc = Vasp(
			xc=defaults['xc'],
			setups=defaults['setups'],
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
			isym=defaults['isym'],
			nsw=0,
			istart=0
			)
	elif run_i == 1:
		calc = Vasp(
			xc=defaults['xc'],
			setups=defaults['setups'],
			ivdw=defaults['ivdw'],
			prec=defaults['prec'],
			algo=defaults['algo'],
			ediff=defaults['ediff'],
			nelm=defaults['nelm'],
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
			setups=defaults['setups'],
			ivdw=defaults['ivdw'],
			prec=defaults['prec'],
			algo=defaults['algo'],
			ediff=defaults['ediff'],
			nelm=defaults['nelm'],
			lreal=defaults['lreal'],
			ncore=defaults['ncore'],
			ismear=defaults['ismear'],
			sigma=defaults['sigma'],
			lcharg=False,
			lwave=True,
			ibrion=2,
			isif=2,
			nsw=defaults['nsw'],
			ediffg=-0.05,
			lorbit=defaults['lorbit'],
			isym=defaults['isym']
			)
	elif run_i == 2:
		calc = Vasp(
			xc=defaults['xc'],
			setups=defaults['setups'],
			encut=defaults['encut'],
			ivdw=defaults['ivdw'],
			prec=defaults['prec'],
			algo=defaults['algo'],
			ediff=1e-6,
			nelm=defaults['nelm'],
			nelmin=defaults['nelmin'],
			lreal=defaults['lreal'],
			ncore=defaults['ncore'],
			ismear=defaults['ismear'],
			sigma=defaults['sigma'],
			lcharg=False,
			lwave=True,
			ibrion=2,
			isif=3,
			nsw=30,
			ediffg=defaults['ediffg'],
			lorbit=defaults['lorbit'],
			isym=defaults['isym']
			)
	elif run_i == 3:
		calc = Vasp(
			xc=defaults['xc'],
			setups=defaults['setups'],
			encut=defaults['encut'],
			ivdw=defaults['ivdw'],
			prec=defaults['prec'],
			algo=defaults['algo'],
			ediff=1e-6,
			nelm=defaults['nelm'],
			nelmin=defaults['nelmin'],
			lreal=defaults['lreal'],
			ncore=defaults['ncore'],
			ismear=defaults['ismear'],
			sigma=defaults['sigma'],
			lcharg=False,
			lwave=True,
			ibrion=2,
			isif=3,
			nsw=30,
			ediffg=defaults['ediffg'],
			lorbit=defaults['lorbit'],
			isym=defaults['isym']
			)
	elif run_i == 4:
		calc = Vasp(
			xc=defaults['xc'],
			setups=defaults['setups'],
			encut=defaults['encut'],
			ivdw=defaults['ivdw'],
			prec=defaults['prec'],
			algo=defaults['algo'],
			ediff=1e-4,
			nelm=defaults['nelm'],
			nelmin=defaults['nelmin'],
			lreal=False,
			ncore=defaults['ncore'],
			ismear=defaults['ismear'],
			sigma=defaults['sigma'],
			lcharg=False,
			lwave=True,
			ibrion=2,
			isif=2,
			nsw=30,
			ediffg=defaults['ediffg'],
			lorbit=defaults['lorbit'],
			isym=defaults['isym']
			)
	elif run_i == 5:
		calc = Vasp(
			xc=defaults['xc'],
			setups=defaults['setups'],
			encut=defaults['encut'],
			ivdw=defaults['ivdw'],
			prec='Accurate',
			algo=defaults['algo'],
			ediff=1e-6,
			nelm=defaults['nelm'],
			lreal=False,
			ncore=defaults['ncore'],
			ismear=defaults['ismear'],
			sigma=defaults['sigma'],
			lcharg=True,
			laechg=True,
			lwave=True,
			nsw=0,
			lorbit=defaults['lorbit'],
			isym=defaults['isym'],
			addgrid=False
			)
	else:
		raise ValueError('Out of range for calculators')

	return calc
