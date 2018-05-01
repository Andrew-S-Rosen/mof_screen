import os
import numpy as np
from shutil import copyfile
from ase.io import read

def nebmake(POSCAR1,POSCAR2,n_images):
	os.system('nebmake.pl '+POSCAR1+' '+POSCAR2+' '+str(n_images))

def neb2dim():
	os.system('vfin.pl neb_temp')
	os.chdir('neb_temp')
	os.system('nebresults.pl')
	copyfile('exts.data','../exts.dat')
	os.chdir('../')
	os.system('neb2dim.pl')
	os.chdir('dim')
	mof = read('POSCAR')	
	return mof

def dimmins(dis):
	os.system('vfin.pl dim_temp')
	os.system('dimmins.pl POSCAR MODECAR '+str(dis))

def nebef(ediffg):
	open('nebef.dat','w').close()
	os.system('nebef.pl > nebef.dat')
	max_F = 0
	if os.stat('nebef.dat').st_size == 0:
		raise ValueError('nebef.dat not written')
	with open('nebef.dat','r') as rf:
		for line in rf:
			line = line.strip()
			max_F_temp = np.fromstring(line,dtype=float,sep=' ')[1]
			if max_F_temp > max_F:
				max_F = max_F_temp
	if max_F <= ediffg:
		neb_conv = True
	else:
		neb_conv = False

	return neb_conv