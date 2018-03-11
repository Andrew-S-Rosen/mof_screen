import numpy as np
from settings import omsdata_path

def get_CN(oms_path):
#read .oms file to get number of OMS

	f = open(oms_path,'r')
	oms_file = f.read()
	n_OMS = int(oms_file.split('OMS=')[1].split('\n')[0])
	f.close()

	return n_OMS

def get_omsex_line(line):
#read line in .omsex file for OMS details

	cus_symbol = line.split(' |')[0]
	cnum = int(line.split('CNUM: ')[1].split('|')[0])
	cus_coord = np.asarray(np.matrix(line.split('COORD: ')[1].split('|')[0][0:-1]))
	coords = np.asarray(np.matrix(line.split('NN: ')[1][0:-4]))

	return cus_symbol, cnum, cus_coord, coords

def get_omsex_data(refcode,n_OMS,mof):
#Get all .omsex data

	cus_coords_all = np.zeros((n_OMS,3))
	ase_cus_idx_all = []
	cus_sym_all = []
	cnums_all = []
	zeo_tol = 0.1
	with open(omsdata_path+refcode+'.omsex','r') as rf:
		for i, line in enumerate(rf):
			cus_sym_temp, cnum_temp, cus_coords_all[i,:], NN_coords_temp = get_omsex_line(line)
			cus_sym_all.append(cus_sym_temp)
			cnums_all.append(cnum_temp)
			if i == 0:
				NN_coords_all = NN_coords_temp
			else:
				NN_coords_all = np.vstack((NN_coords_all,NN_coords_temp))
			for j, element in enumerate(mof):
				if (sum(cus_coords_all[i,:] >= element.position-zeo_tol) == 3
					and sum(cus_coords_all[i,:] <= element.position+zeo_tol) == 3):
					ase_cus_idx_all.append(j)
					break
			if len(ase_cus_idx_all) < i+1:
				raise ValueError('ERROR with '+refcode+': a zeo++ OMS (#'+
					str(i)+') is not in same spot as in ASE')
				
	return cnums_all, cus_coords_all, ase_cus_idx_all, cus_sym_all, NN_coords_all
