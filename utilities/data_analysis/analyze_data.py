import os
from ase.io import read
import pandas as pd
import numpy as np
from ase.geometry import cell_to_cellpar
import matplotlib.pyplot as plt

phase1_results = '/projects/p30470/phase1_results/'
phase2_results = '/projects/p30148/vasp_jobs/MOFs/oxidized_oms/results/'
phase3_results = '/projects/p30148/vasp_jobs/MOFs/phase3/results/'
phase4_results = '/projects/p30148/vasp_jobs/MOFs/phase4/results/'
gas_path = '/projects/p30148/vasp_jobs/MOFs/gasphase/'
oms_path = '/projects/p30148/vasp_jobs/structures/CoRE1-DFT-Clean/oms_data/'

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
				kpts = np.array(line.strip())
			else:
				break
	return kpts, gamma

def get_gas_dict(gas_path):
	h2o = read(gas_path+'h2o/OUTCAR')
	o2 = read(gas_path+'o2/OUTCAR')
	ch4 = read(gas_path+'ch4/OUTCAR')
	E_H2O = h2o.get_potential_energy()
	E_O2 = o2.get_potential_energy()
	E_CH4 = ch4.get_potential_energy()
	gas_dict = {'H2O':E_H2O,'O2':E_O2,'CH4':E_CH4}
	return gas_dict

def get_NNs(screen_path,mof_name):
	try:
		with open(screen_path+'NN_list.txt','r') as rf:
			for line in rf:
				if line.split('|')[0] == mof_name:
					NNs = line.split('|')[-1].strip()
					break
	except:
		NNs = np.nan
	return NNs

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

def get_mof_df(screen_path):
	refcode_list = []
	E_list = []
	net_magmom_list = []
	kpts_list = []
	gamma_list = []
	V_list = []
	lattice_list = []
	formula_list = []
	NNs_list = []
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
		NNs = get_NNs(screen_path,mof_name)
		try:
			net_magmom = mof.get_magnetic_moment()
		except:
			net_magmom = 0.0
		V = mof.get_volume()
		kpts, gamma = get_kpts(mof_path_temp)
		lattice = cell_to_cellpar(mof.cell)
		formula = mof.get_chemical_formula()

		refcode_list.append(refcode)
		kpts_list.append(kpts)
		gamma_list.append(gamma)
		E_list.append(E)
		net_magmom_list.append(net_magmom)
		V_list.append(V)
		lattice_list.append(lattice)
		formula_list.append(formula)
		NNs_list.append(NNs)
	mof_df = pd.DataFrame({'Name':refcode_list,'E':E_list,'Formula':formula_list,'NNs':NNs_list,'Mag':net_magmom_list,'Cell':lattice_list,'V':V_list,'kpts':kpts_list,'gamma':gamma_list})
	return mof_df

def get_NN_list(screen_path):
	names = []
	NNs = []
	with open(screen_path+'NN_list.txt','r') as rf:
		for line in rf:
			split_line = line.split('|').strip()
			names.append(split_line[0])
			NNs.append(split_line[1])
	return names, NNs

def plot_descriptor(E_H_list,E_TS_list):
	plt.figure()
	plt.plot(E_H_list,E_TS_list,'ro')
	plt.xlabel('E_H (eV)')
	plt.ylabel('E_TS (eV)')
	plt.savefig('E_TS.png')
	plt.close()

def write_descriptor(refcodes,E_H_list,E_TS_list,Ea_list):
	with open('E_TS.csv','w') as wf:
		for i,refcode in enumerate(refcodes):
			wf.write(refcode+','+str(E_H_list[i])+','+str(E_TS_list[i])+','+str(Ea_list[i])+'\n')

def get_df(phase,df_path):
	pckl_name = 'phase'+str(phase)+'_df.pkl'
	if os.path.exists(pckl_name):
		mof_df = pd.read_pickle(pckl_name)
	else:
		mof_df = get_mof_df(df_path)
		mof_df.to_pickle(pckl_name)
	return mof_df

def clean_df(phase2_df,phase3_df):
	for i, row in phase3_df.iterrows():
		phase3_NNs = row.NNs
		phase3_name = row.Name
		phase2_name = phase3_name.rsplit('_spin',1)[0]
		phase2_NNs = phase2_df[phase2_df['Name'].str.match(phase2_name)].NNs.tolist()[0]
		if phase2_NNs != phase3_NNs:
			phase3_df = phase3_df[phase3_df.Name != phase3_name]
			phase2_df = phase2_df[phase2_df.Name != phase2_name]
	return phase2_df, phase3_df

def get_E_TS(phase2_df,phase3_df,gas_dict):
	E_TS_list = []
	E_H_list = []
	Ea_list = []
	phase2_names = []
	H_ref = 0.5*gas_dict['H2O']-0.25*gas_dict['O2']
	for i, row in phase3_df.iterrows():
		try:
			phase3_E = row.E
			phase3_name = row.Name
			phase2_name = phase3_name.rsplit('_spin',1)[0]
			phase2_E = float(phase2_df[phase2_df['Name'].str.match(phase2_name)].E)
			E_H = phase3_E-phase2_E-H_ref
			E_TS = 0.75*E_H+1.09
			E_H_list.append(E_H)
			E_TS_list.append(E_TS)
			phase4_E = float(phase4_df[phase4_df['Name'].str.match(phase2_name)].E)
			Ea = E_TS + phase2_E + gas_dict['CH4'] - phase4_E
			Ea_list.append(Ea)
			phase2_names.append(phase2_name)
		except:
			pass
	return phase2_names, E_H_list, E_TS_list, Ea_list

def get_E_ads(phase2_df,phase4_df,gas_dict):
	E_ads_list = []
	phase2_names = []
	for i, row in phase4_df.iterrows():
		try:
			phase4_E = row.E
			phase4_name = row.Name
			phase2_name = phase4_name.rsplit('_spin',1)[0]
			phase2_E = float(phase2_df[phase2_df['Name'].str.match(phase2_name)].E)
			E_ads = phase4_E - (phase2_E + gas_dict['CH4'])
			E_ads_list.append(E_ads)
			phase2_names.append(phase2_name)
		except:
			pass
	return phase2_names, E_ads_list

def get_E_ox(phase1_df,phase2_df,gas_dict):
	E_ox_list = []
	phase1_names = []
	for i, row in phase2_df.iterrows():
		phase2_E = row.E
		phase2_name = row.Name
		phase1_name = phase2_name.split('_spin')[0]
		phase1_E = float(phase1_df[phase1_df['Name'].str.match(phase1_name)].E)
		E_ox = phase2_E - (phase1_E + 0.5*gas_dict['O2'])
		E_ox_list.append(E_ox)
		phase1_names.append(phase1_name)
	return phase1_names, E_ox_list

gas_dict = get_gas_dict(gas_path)
phase1_df = get_df(1,phase1_results)
phase2_df = get_df(2,phase2_results)
phase3_df = get_df(3,phase3_results)
phase4_df = get_df(4,phase4_results)
phase2_df, phase3_df = clean_df(phase2_df,phase3_df)
phase2_TS_names, E_H_list, E_TS_list, Ea_list = get_E_TS(phase2_df,phase3_df,gas_dict)
plot_descriptor(E_H_list,E_TS_list)
phase2_ads_names, E_ads_list = get_E_ads(phase2_df,phase4_df,gas_dict)
phase1_ox_names, E_ox_list = get_E_ox(phase1_df,phase2_df,gas_dict)
write_descriptor(phase2_TS_names,E_H_list,E_TS_list,Ea_list)
