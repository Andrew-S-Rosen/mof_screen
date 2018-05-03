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
	'ismear': 0,
	'sigma': 0.01,
	'nsw': 500,
	'ediffg': -0.03,
	'lorbit': 11,
	'isym': 0,
	'setups':{'base':'recommended','Li':''}
	}

def calcs(calc_name):
	"""
	Define the default calculators for relaxations
	Note: it should not include the kpts, gamma, or images keywords!
	Args:
		calc_name (string): name of calculator
	Returns:
		calc (dict): ASE Vasp calculator dictionary
	"""
	if calc_name == 'scf_test':
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
			ismear=defaults['ismear'],
			sigma=defaults['sigma'],
			lcharg=False,
			lwave=True,
			lorbit=defaults['lorbit'],
			isym=defaults['isym'],
			nsw=0,
			istart=0
			)
	elif calc_name == 'ase_bfgs':
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
			ismear=defaults['ismear'],
			sigma=defaults['sigma'],
			lcharg=False,
			lwave=True,
			lorbit=defaults['lorbit'],
			isym=defaults['isym']
			)
	elif calc_name == 'isif2_lowacc':
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
			ismear=defaults['ismear'],
			sigma=defaults['sigma'],
			lcharg=False,
			lwave=True,
			ibrion=2,
			isif=2,
			nsw=250,
			ediffg=-0.05,
			lorbit=defaults['lorbit'],
			isym=defaults['isym']
			)
	elif calc_name == 'isif2_medacc':
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
	elif calc_name == 'isif2_highacc':
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
	elif calc_name == 'isif3_lowacc':
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
	elif calc_name == 'isif3_highacc':
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
	elif calc_name == 'final_spe':
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
	elif calc_name == 'cineb_lowacc':
		calc = Vasp(
			xc=defaults['xc'],
			setups=defaults['setups'],
			ivdw=defaults['ivdw'],
			prec=defaults['prec'],
			algo=defaults['algo'],
			ediff=1e-6,
			nelm=defaults['nelm'],
			nelmin=defaults['nelmin'],
			lreal=defaults['lreal'],
			ismear=defaults['ismear'],
			sigma=defaults['sigma'],
			lcharg=False,
			lwave=False,
			ibrion=3,
			potim=0,
			iopt=1,
			nsw=defaults['nsw'],
			ediffg=-0.5,
			lclimb=True,
			lorbit=defaults['lorbit'],
			isym=defaults['isym'],
			ichain=0
			)
	elif calc_name == 'dimer_lowacc':
		calc = Vasp(
			xc=defaults['xc'],
			setups=defaults['setups'],
			ivdw=defaults['ivdw'],
			prec=defaults['prec'],
			algo=defaults['algo'],
			ediff=1e-8,
			nelm=defaults['nelm'],
			nelmin=defaults['nelmin'],
			lreal=defaults['lreal'],
			ismear=defaults['ismear'],
			sigma=defaults['sigma'],
			lcharg=False,
			lwave=True,
			ibrion=3,
			potim=0,
			iopt=2,
			nsw=defaults['nsw']*4,
			ediffg=-0.05,
			lorbit=defaults['lorbit'],
			isym=defaults['isym'],
			ichain=2
			)
	elif calc_name == 'dimer_medacc':
		calc = Vasp(
			xc=defaults['xc'],
			setups=defaults['setups'],
			ivdw=defaults['ivdw'],
			prec=defaults['prec'],
			algo=defaults['algo'],
			ediff=1e-8,
			nelm=defaults['nelm'],
			nelmin=defaults['nelmin'],
			lreal=defaults['lreal'],
			ismear=defaults['ismear'],
			sigma=defaults['sigma'],
			lcharg=False,
			lwave=True,
			ibrion=3,
			potim=0,
			iopt=2,
			nsw=defaults['nsw']*4,
			ediffg=defaults['ediffg'],
			lorbit=defaults['lorbit'],
			isym=defaults['isym'],
			ichain=2
			)
	elif calc_name == 'dimer_highacc':
		calc = Vasp(
			xc=defaults['xc'],
			encut=defaults['encut'],
			setups=defaults['setups'],
			ivdw=defaults['ivdw'],
			prec=defaults['prec'],
			algo=defaults['algo'],
			ediff=1e-8,
			nelm=defaults['nelm'],
			nelmin=defaults['nelmin'],
			lreal=defaults['lreal'],
			ismear=defaults['ismear'],
			sigma=defaults['sigma'],
			lcharg=False,
			lwave=True,
			ibrion=3,
			potim=0,
			iopt=2,
			nsw=defaults['nsw']*4,
			ediffg=defaults['ediffg'],
			lorbit=defaults['lorbit'],
			isym=defaults['isym'],
			ichain=2
			)
	else:
		raise ValueError('Out of range for calculators')

	return calc