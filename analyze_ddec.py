import re
import os
import numpy as np

results_path = '/projects/p30148/vasp_jobs/MOFs/reoptimized_core1/results/'

def get_ddec(filepath):
	cat_metals = ['Li', 'Be', 'Na', 'Mg', 'Al', 'K', 'Ca', 'Sc', 'Ti', 
		'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Rb', 'Sr', 
		'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 
		'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 
		'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu', 'Hf', 'Ta', 'W', 'Re', 'Os', 
		'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Fr', 'Ra', 'Ac', 'Th', 
		'Pa', 'U', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm', 'Md', 
		'No', 'Lr', 'Rf', 'Db', 'Sg', 'Bh', 'Hs', 'Mt', 'Ds', 'Rg', 'Cn']
	check = False
	charges = []
	elements = []
	with open(filepath+'DDEC6_even_tempered_net_atomic_charges.xyz','r') as rf:
		for line in rf:
			line = line.strip()
			if 'quadrupole moment tensor' in line:
				check = True
				continue
			if check == True:
				if not line:
					break
				else:
					line = line.strip()
					element_str = re.findall('[a-z]',line,re.I)
					element = ''
					for letter in element_str:
						element += letter
					elements.append(element)
					line = line.split(element_str[-1])[-1]
					charges.append(np.fromstring(line,dtype=float,sep=' ')[3])
	bad_charges = []
	bad_elements = []
	for i, element in enumerate(elements):
		if element in cat_metals:
			charge = charges[i]
			if charge < 0:
				bad_charges.append(charge)
				bad_elements.append(element)
	return bad_charges, bad_elements

def analyze_ddecs():
	bad_refcodes = []
	for refcode in os.listdir(results_path):
		if os.path.isdir(results_path+refcode):
			subdirs = os.listdir(results_path+refcode)
			for subdir in subdirs:
				path = results_path+refcode+'/'+subdir+'/ddec/'
				if os.path.exists(path):
					bad_charges, bad_elements = get_ddec(path)
					if len(bad_elements) != 0:
						if refcode not in bad_refcodes:
							bad_refcodes.append(refcode)
	return bad_refcodes

bad_refcodes = analyze_ddecs()
print(bad_refcodes)
