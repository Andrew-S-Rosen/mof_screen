import os
from ase.io import read, write
mofpath = '/projects/p30148/vasp_jobs/structures/CoRE1-DFT-OMS-v2/'
results = '/projects/p30148/vasp_jobs/MOFs/reoptimized_core1/results/'

startcifs = []
for filename in os.listdir(mofpath):
	filename = filename.strip()
	if len(filename.split('.cif')) == 2:
		refcode = filename.split('.cif')[0]
		startcifs.append(filename)
startcifs = sorted(startcifs)
for cif in startcifs:
	refcode = cif.split('.cif')[0]
	finalpath = results+refcode+'/final_spe/'
	if os.path.isdir(finalpath):
		spin1path = finalpath+'spin1/'
		spin2path = finalpath+'spin2/'
		mof1 = read(spin1path+'CONTCAR')
		if os.path.isdir(spin2path):
			mof2 = read(spin2path+'CONTCAR')
			E_mof1 = read(spin1path+'OUTCAR')
			E_mof1 = E_mof1.get_potential_energy()
			E_mof2 = read(spin2path+'OUTCAR')
			E_mof2 = E_mof2.get_potential_energy()
			if E_mof2 < E_mof1:
				write(refcode+'_spin2.cif',mof2)
			else:
				write(refcode+'_spin1.cif',mof1)
		else:
			write(refcode+'_spin1.cif',mof1)
	else:
		print(refcode+': not completed')
