from __future__ import division
import os
from ase.io import read
import pandas as pd
import numpy as np
from ase.geometry import cell_to_cellpar
from ase.thermochemistry import IdealGasThermo
import matplotlib.pyplot as plt
from ase.vibrations import Vibrations
from pymatgen.analysis.local_env import MinimumVIRENN
import re
from pymatgen.io import ase as pm_ase

phase1_results = '/projects/p30470/phase1_results/'
phase2_results = '/projects/p30148/vasp_jobs/MOFs/oxidized_oms/results/'
phase3_results = '/projects/p30148/vasp_jobs/MOFs/phase3/results/'
phase4_results = '/projects/p30148/vasp_jobs/MOFs/phase4/results/'
gas_path = '/projects/p30148/vasp_jobs/MOFs/gasphase/'
T_gas = 150+273.15
P_gas = 1e5

def get_ddec(filepath):
	check = False
	charges = []
	i = 0
	with open(filepath+'/DDEC6_even_tempered_net_atomic_charges.xyz','r') as rf:
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
					line = line.split(element_str[-1])[-1]
					charges.append(np.fromstring(line,dtype=float,sep=' ')[3])
					i += 1
	return charges


def get_kpts(file_path):
	with open(file_path+'/KPOINTS','r') as rf:
		for i, line in enumerate(rf):
			if i < 2:
				continue
			elif i == 2:
				if 'gamma' in line.lower():
					gamma = True
				else:
					gamma = False
			elif i == 3:
				kpts = [int(num) for num in line.strip().split(' ')]
			else:
				break
	return kpts, gamma

def get_gas_df(gas_path,T,P):
	basename = 'gas_df'
	pckl_name = basename+'.pkl'
	csv_name = basename+'.csv'
	if os.path.exists(pckl_name):
		gas_df = pd.read_pickle(pckl_name)
	else:
		E_H2O, G_H2O = get_gas_data(gas_path,'h2o',T,P,'nonlinear',2,0)
		E_O2, G_O2 = get_gas_data(gas_path,'o2',T,P,'linear',2,1)
		E_CH4, G_CH4 = get_gas_data(gas_path,'ch4',T,P,'nonlinear',12,0)
		E_N2O, G_N2O = get_gas_data(gas_path,'n2o',T,P,'linear',1,0)
		E_N2, G_N2 = get_gas_data(gas_path,'n2',T,P,'linear',2,0)
		E_H2, G_H2 = get_gas_data(gas_path,'h2',T,P,'linear',2,0)
		gas_names = ['H2O','O2','CH4','N2O','N2','H2']
		E_list = [E_H2O,E_O2,E_CH4,E_N2O,E_N2,E_H2]
		G_list = [G_H2O,G_O2,G_CH4,G_N2O,G_N2,G_H2]
		gas_df = pd.DataFrame({'Name':gas_names,'E':E_list,'G':G_list})
		gas_df.set_index('Name',inplace=True)
		gas_df.to_pickle(pckl_name)
		gas_df.to_csv(csv_name)
	return gas_df

def get_gas_data(gas_path,gas_name,T,P,geom,sym,spin):
	pwd = os.getcwd()
	gas_E = read(gas_path+gas_name+'/OUTCAR')
	E = gas_E.get_potential_energy()
	gas = read(gas_path+gas_name+'/'+gas_name+'.xyz')
	os.chdir(gas_path+gas_name+'/vib')
	vib = Vibrations(gas)
	vib_energies = vib.get_energies()
	thermo = IdealGasThermo(
		vib_energies=vib_energies,
		potentialenergy=E,
		atoms=gas,
		geometry=geom,
		symmetrynumber=sym,
		spin=spin
		)
	G = thermo.get_gibbs_energy(temperature=T,pressure=P,verbose=False)
	os.chdir(pwd)
	return E, G

def get_NNs(mof,idx):
	bridge = pm_ase.AseAtomsAdaptor()
	struct = bridge.get_structure(mof)
	nn_object = MinimumVIRENN()
	neighbors = nn_object.get_nn_info(struct,idx)
	NN_idx = []
	if neighbors:
		for neighbor in neighbors:
			NN_idx.append(neighbor['site_index'])
	return NN_idx

def check_bad(screen_path,mof_name):
	isbad = False
	try:
		with open(screen_path+'bad_ads_addition.txt','r') as rf:
			for line in rf:
				if line.strip() == mof_name:
					isbad = True
	except:
		pass
	return isbad

def get_mof_data(screen_path,phase):
	refcode_list = []
	E_list = []
	net_magmom_list = []
	kpts_list = []
	gamma_list = []
	V_list = []
	lattice_list = []
	formula_list = []
	OMS_list = []
	n_atoms = []
	C_dist_list = []
	H_dist_list = []
	O_mag_list = []
	M_mag_list = []
	oxo_cnum_list = []
	MO_dist_list = []
	O_ddec6_list = []
	M_ddec6_list = []
	folder_list = os.listdir(screen_path)
	folder_list.sort()
	for folder in folder_list:
		if not os.path.exists(screen_path+'/'+folder+'/final_spe/'):
			continue
		refcode = folder
		E = np.inf
		for spin_result in os.listdir(screen_path+'/'+folder+'/final_spe/'):
			mof_path_temp = screen_path+folder+'/final_spe/'+spin_result
			mof_temp = read(mof_path_temp+'/OUTCAR')
			E_temp = mof_temp.get_potential_energy()
			if E_temp < E:
				E = E_temp
				mof_path = mof_path_temp
		mof_name = folder+'_'+spin_result
		isbad = check_bad(screen_path,mof_name)
		if isbad == True:
			continue
		mof = read(mof_path+'/OUTCAR')
		try:
			net_magmom = np.round(mof.get_magnetic_moment(),1)
		except:
			net_magmom = 0.0
		if phase == 3:
			O_idx = [atom.index for atom in mof if atom.symbol == 'O'][-1]
			H_idx = [atom.index for atom in mof if atom.symbol == 'H'][-1]
			H_dist = np.round(mof.get_distance(O_idx,H_idx,mic=True,vector=False),2)
			H_dist_list.append(H_dist)
			del mof[H_idx]
		elif phase == 4:
			O_idx = [atom.index for atom in mof if atom.symbol == 'O'][-1]
			C_idx = [atom.index for atom in mof if atom.symbol == 'C'][-1]
			H_idx = [atom.index for atom in mof if atom.symbol == 'H'][-4:]
			H_dist = np.round(np.min(mof.get_distances(O_idx,H_idx,mic=True,vector=False)),2)
			C_dist = np.round(mof.get_distance(O_idx,C_idx,mic=True,vector=False),2)
			H_dist_list.append(H_dist)
			C_dist_list.append(C_dist)
			del mof[[C_idx]+H_idx]
		if phase >= 2:
			O_idx = [atom.index for atom in mof if atom.symbol == 'O'][-1]
			ddec_charges = np.round(get_ddec(mof_path+'/ddec'),2)
			OMS_idx = get_NNs(mof,O_idx)
			OMS = mof[OMS_idx].get_chemical_symbols()
			OMS_list.append(OMS)
			oxo_cnum_list.append(len(OMS_idx))
			O_ddec6_list.append(ddec_charges[O_idx])
			for M_idx in OMS_idx:
				MO_dist = np.round(mof.get_distance(M_idx,O_idx,mic=True,vector=False),2)
				M_ddec6 = ddec_charges[M_idx]
			M_ddec6_list.append(M_ddec6)
			MO_dist_list.append(MO_dist)
			try:
				mof_magmoms = mof.get_magnetic_moments()
				O_mag = mof_magmoms[O_idx]
				M_mag = mof_magmoms[OMS_idx]
				O_mag_list.append(O_mag)
				M_mag_list.append(M_mag)
			except:
				O_mag = 0.0
				M_mag = 0.0
				O_mag_list.append(O_mag)
				M_mag_list.append(M_mag)
		V = np.round(mof.get_volume(),1)
		kpts, gamma = get_kpts(mof_path_temp)
		lattice = np.round(cell_to_cellpar(mof.cell),2).tolist()
		formula = mof.get_chemical_formula()
		n_atom = len(mof)

		refcode_list.append(refcode)
		kpts_list.append(kpts)
		gamma_list.append(gamma)
		E_list.append(E)
		net_magmom_list.append(net_magmom)
		V_list.append(V)
		lattice_list.append(lattice)
		formula_list.append(formula)
		n_atoms.append(n_atom)
	data_dict = {'Name':refcode_list,'E':E_list,'Formula':formula_list,'n_atoms':n_atoms,'Mag':net_magmom_list,'Cell':lattice_list,'V':V_list,'kpts':kpts_list,'gamma':gamma_list}
	if phase >= 2:
		data_dict['O_mag'] = O_mag_list
		data_dict['M'] = OMS_list
		data_dict['M_mag'] = M_mag_list
		data_dict['Oxo_cnum'] = oxo_cnum_list
		data_dict['MO_dist'] = MO_dist_list
		data_dict['O_ddec6'] = O_ddec6_list
		data_dict['M_ddec6'] = M_ddec6_list
	if phase == 3:
		data_dict['H_dist'] = H_dist_list
	if phase == 4:
		data_dict['H_dist'] = H_dist_list
		data_dict['C_dist'] = C_dist_list
	mof_df = pd.DataFrame(data_dict)
	mof_df.set_index('Name',inplace=True)
	return mof_df

def plot_descriptor(TS_df):
	TS_df.to_csv('E_TS.csv')
	plt.figure()
	plt.plot(TS_df['E_H'],TS_df['E_TS'],'ro')
	plt.xlabel('E_H (eV)')
	plt.ylabel('E_TS (eV)')
	plt.savefig('E_TS.png',bbox_inches='tight')
	plt.close()

def get_mof_df(phase,df_path):
	basename = 'phase'+str(phase)+'_df'
	pckl_name = basename+'.pkl'
	csv_name = basename+'.csv'
	if os.path.exists(pckl_name):
		mof_df = pd.read_pickle(pckl_name)
	else:
		mof_df = get_mof_data(df_path,phase)
		mof_df.to_pickle(pckl_name)
		mof_df.to_csv(csv_name)
	return mof_df

def clean_df(phase2_df,phase3_df):
	for phase3_name, row in phase3_df.iterrows():
		phase3_NNs = row.M
		phase2_name = phase3_name.rsplit('_spin',1)[0]
		phase2_NNs = phase2_df[phase2_df.index.str.match(phase2_name)].M.tolist()[0]
		if phase2_NNs != phase3_NNs:
			phase3_df = phase3_df[phase3_df.index != phase3_name]
	return phase3_df

def get_E_TS(phase2_df,phase3_df,phase4_df,gas_df):
	E_TS_list = []
	E_H_list = []
	Ea_list = []
	phase2_names = []
	H_ref = 0.5*gas_df.loc['H2'].E
	for phase3_name, row in phase3_df.iterrows():
		try:
			phase3_E = row.E
			phase2_name = phase3_name.rsplit('_spin',1)[0]
			phase2_E = float(phase2_df[phase2_df.index.str.match(phase2_name)].E)
			E_H = phase3_E-phase2_E-H_ref
			E_TS = 0.75*E_H+1.96
			phase4_E = float(phase4_df[phase4_df.index.str.match(phase2_name)].E)
			phase4_E_ref = phase4_E - (phase2_E + gas_df.loc['CH4'].E)
			Ea = E_TS - phase4_E_ref
			if Ea < 0:
				Ea = 0
			E_H_list.append(E_H)
			E_TS_list.append(E_TS)
			Ea_list.append(Ea)
			phase2_names.append(phase2_name)
		except:
			pass
	TS_df = pd.DataFrame({'Name':phase2_names,'E_H':E_H_list,'E_TS':E_TS_list,'E_a':Ea_list})
	TS_df = TS_df[['Name','E_H','E_TS','E_a']]
	TS_df.set_index('Name',inplace=True)
	return TS_df

def get_E_ads(phase2_df,phase4_df,gas_df):
	E_ads_list = []
	phase2_names = []
	for phase4_name, row in phase4_df.iterrows():
		try:
			phase4_E = row.E
			phase2_name = phase4_name.rsplit('_spin',1)[0]
			phase2_E = float(phase2_df[phase2_df.index.str.match(phase2_name)].E)
			E_ads = phase4_E - (phase2_E + gas_df.loc['CH4'].E)
			E_ads_list.append(E_ads)
			phase2_names.append(phase2_name)
		except:
			pass
	ads_df = pd.DataFrame({'Name':phase2_names,'E_ads':E_ads_list})
	ads_df = ads_df[['Name','E_ads']]
	ads_df.set_index('Name',inplace=True)
	return ads_df

def get_E_ox(phase1_df,phase2_df,gas_df):
	E_ox_list = []
	phase1_names = []
	for phase2_name, row in phase2_df.iterrows():
		phase2_E = row.E
		phase1_name = phase2_name.split('_spin')[0]
		phase1_E = float(phase1_df[phase1_df.index.str.match(phase1_name)].E)
		E_ox = (phase2_E+gas_df.loc['N2'].E) - (phase1_E+gas_df.loc['N2O'].E)
		E_ox_list.append(E_ox)
		phase1_names.append(phase1_name)
	ox_df = pd.DataFrame({'Name':phase1_names,'E_ox':E_ox_list})
	ox_df = ox_df[['Name','E_ox']]
	ox_df.set_index('Name',inplace=True)
	return ox_df

def plot_histograms(Ea_list,E_ox_list,E_ads_list):
	plt.figure()
	make_hist(E_ads_list,'r','auto')
	make_hist(Ea_list,'g',30)
	make_hist(E_ox_list,'b','auto')
	plt.minorticks_on()
	plt.xlabel(r'$\Delta E$ (eV)')
	plt.ylabel('Frequency')
	plt.legend([r'$\Delta E_{ads}$',r'$\Delta E_{a}$',r'$\Delta E_{ox}$'])
	plt.savefig('hist.png',bbox_inches='tight')
	plt.close()

def make_hist(data,color,bins):
	yhist, xhist, patches = plt.hist(data,bins=bins,color=color,alpha=0.5,edgecolor='k')

gas_df = get_gas_df(gas_path,T_gas,P_gas)
phase1_df = get_mof_df(1,phase1_results)
phase2_df = get_mof_df(2,phase2_results)
phase3_df = get_mof_df(3,phase3_results)
phase4_df = get_mof_df(4,phase4_results)
phase3_df = clean_df(phase2_df,phase3_df)
phase4_df = clean_df(phase2_df,phase4_df)
TS_df = get_E_TS(phase2_df,phase3_df,phase4_df,gas_df)
plot_descriptor(TS_df)
ads_df = get_E_ads(phase2_df,phase4_df,gas_df)
ox_df = get_E_ox(phase1_df,phase2_df,gas_df)
plot_histograms(TS_df['E_a'],ox_df['E_ox'],ads_df['E_ads'])