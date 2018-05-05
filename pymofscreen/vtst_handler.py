import os
import numpy as np
import time
from shutil import copyfile, rmtree, move
from ase.io import read, write
from pymofscreen.janitor import clean_files

def nebmake(initial_atoms,final_atoms,n_images):
	pwd = os.getcwd()
	neb_path = os.path.join(pwd,'neb')
	if os.path.exists(neb_path):
		rmtree(neb_path)
	os.makedirs(neb_path)
	os.chdir(neb_path)
	if n_images < 10:
		last_image = '0'+str(n_images+1)
	else:
		last_image = str(n_images+1)
	write(os.path.join(neb_path,'POSCAR1'),initial_atoms,format='vasp')
	write(os.path.join(neb_path,'POSCAR2'),final_atoms,format='vasp')
	os.system('nebmake.pl POSCAR1 POSCAR2 '+str(n_images))
	write_dummy_outcar(os.path.join(neb_path,'00','OUTCAR'),initial_atoms.get_potential_energy())
	write_dummy_outcar(os.path.join(neb_path,last_image,'OUTCAR'),final_atoms.get_potential_energy())

def write_dummy_outcar(name,E):
	with open(name,'w') as wf:
		wf.write('  energy  without entropy=                   energy(sigma->0) =     '+str(E)+'\n')

def neb2dim():
	pwd = os.getcwd()
	neb_path = os.path.join(pwd,'neb')
	os.chdir(neb_path)
	os.system('vfin.pl neb_fin')
	time.sleep(5)
	neb_fin_path = os.path.join(neb_path,'neb_fin')
	os.chdir(neb_fin_path)
	os.system('nebresults.pl')
	copyfile(os.path.join(neb_fin_path,'exts.dat'),os.path.join(neb_path,'exts.dat'))
	os.chdir(neb_path)
	if os.stat(os.path.join(neb_path,'exts.dat')).st_size == 0:
		raise ValueError('Error with exts.dat file')
	os.system('neb2dim.pl')
	old_dim_path = os.path.join(neb_path,'dim')
	new_dim_path = os.path.join(pwd,'dim')
	move(old_dim_path,new_dim_path)
	os.chdir(new_dim_path)
	mof = read('POSCAR')	
	return mof

def dimmins(dis):
	#probably wont work without all files
	os.system('vfin.pl dim_fin')
	rmtree('dim_fin')
	os.system('dimmins.pl POSCAR MODECAR '+str(dis))

def nebef(ediffg):
	ediffg = abs(ediffg)
	clean_files(['POSCAR'])
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