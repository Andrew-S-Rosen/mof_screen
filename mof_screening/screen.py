import os
import numpy as np
from compute_environ import get_nprocs, choose_vasp_version
from writers import pprint, write_success, write_errors, write_energy
from settings import basepath
from kpts_handler import get_kpts, get_gpt_version
from janitor import clean_files, manage_restart_files
from runner import mof_run, prep_next_run, mof_bfgs_run
from ase.io import read
from cif_handler import cif_to_mof
from magmom_handler import set_initial_magmoms, continue_magmoms, get_mag_indices
from calculators import calcs_ads, calcs_vol, defaults
from error_handler import get_warning_msgs

def run_ads_screen(cif_files):
#Run high-throughput screening

	#Files, spin levels, and accuracy levels to iterate over
	vasp_files = ['INCAR','POSCAR','KPOINTS','POTCAR','OUTCAR',
	'CONTCAR','CHGCAR','AECCAR0','AECCAR2','WAVECAR']
	spin_levels = ['spin1','spin2']
	acc_levels = ['scf_test','isif2_lowacc','isif2_medacc','final','final_spe']
	f_tol = np.abs(defaults['ediffg'])

	#for each CIF file, optimize the structure
	for cif_file in cif_files:

		nprocs, ppn = get_nprocs()

		refcode = cif_file.split('.cif')[0]
		pprint('***STARTING '+refcode+'***')

		#Make sure MOF isn't running on other process
		working_cif_path = basepath+'working/'+cif_file
		if os.path.isfile(working_cif_path) == True:
			pprint('SKIPPED: Running on another process')
			continue

		#Partial paths to write the OUTCARs
		results_partial_paths = []
		error_outcar_partial_paths = []
		for acc_level in acc_levels:
			results_partial_paths.append(basepath+'results/'+refcode+'/'+acc_level)
			error_outcar_partial_paths.append(basepath+'errors/'+refcode+'/'+acc_level)
		spin1_final_mof_path = results_partial_paths[-1]+'/'+spin_levels[0]+'/OUTCAR'

		#Get the kpoints
		kpts_lo, gamma = get_kpts(cif_file,defaults['kppa_lo'])
		kpts_hi, gamma = get_kpts(cif_file,defaults['kppa_hi'])
		defaults['gamma'] = gamma
		defaults['kpts_lo'] = kpts_lo
		defaults['kpts_hi'] = kpts_hi

		#for each spin level, optimize the structure
		for spin_level in spin_levels:

			#***********PREP FOR RUN***********
			calc_swaps = []
			outcar_paths = []
			error_outcar_paths = []
			run_i = 0
			clean_files(vasp_files)
			for results_partial_path in results_partial_paths:
				outcar_paths.append(results_partial_path+'/'+spin_level+'/OUTCAR')
			for error_outcar_partial_path in error_outcar_partial_paths:
				error_outcar_paths.append(error_outcar_partial_path+'/'+spin_level+'/OUTCAR')

			#***********SCF TEST************
			acc_level = acc_levels[run_i]
			if os.path.isfile(outcar_paths[run_i]) != True and os.path.isfile(error_outcar_paths[run_i]) != True:
				if os.path.isfile(spin1_final_mof_path):
					mof = read(spin1_final_mof_path)
				else:
					mof = cif_to_mof(cif_file)
				mof = set_initial_magmoms(mof,spin_level)
				gpt_version, nprocs = get_gpt_version(kpts_lo,len(mof),nprocs,ppn)
				choose_vasp_version(gpt_version,nprocs,calc_swaps)
				pprint('Running '+spin_level+', '+acc_level)
				mof, calc_swaps = mof_run(mof,calcs_ads(run_i),cif_file,gpt_version,nprocs,calc_swaps)
				if mof != None:
					write_success(refcode,spin_level,acc_level,vasp_files,cif_file)
				else:
					pprint('^ VASP crashed')
					write_errors(refcode,spin_level,acc_level,vasp_files,cif_file)
			elif os.path.isfile(outcar_paths[run_i]) == True:
				pprint('COMPLETED: '+spin_level+', '+acc_level)
			mof, run_i, skip_spin2 = prep_next_run(acc_level,run_i,refcode,spin_level)
			if mof == None:
				pprint('Skipping rest because of errors')
				break
			warnings = get_warning_msgs(outcar_paths[run_i-1])
			calc_swaps.extend(warnings)

			#***********ISIF 2 (lowacc)************
			acc_level = acc_levels[run_i]
			if os.path.isfile(outcar_paths[run_i-1]) == True and os.path.isfile(outcar_paths[run_i]) != True and os.path.isfile(error_outcar_paths[run_i]) != True:
				if os.path.isfile(spin1_final_mof_path):
					mof = read(spin1_final_mof_path)
				else:
					mof = cif_to_mof(cif_file)
				manage_restart_files(results_partial_paths[run_i-1]+'/'+spin_level)
				mof = set_initial_magmoms(mof,spin_level)
				gpt_version, nprocs = get_gpt_version(kpts_lo,len(mof),nprocs,ppn)
				choose_vasp_version(gpt_version,nprocs,calc_swaps)
				pprint('Running '+spin_level+', '+acc_level)
				steps = 100
				fmax = 5.0
				mof, dyn, calc_swaps = mof_bfgs_run(mof,calcs_ads(run_i),cif_file,calc_swaps,steps,fmax)
				if mof != None and dyn and mof.calc.scf_converged == True:
					loop_i = 0
					converged = False
					clean_files(['opt.traj'])
					while mof != None and loop_i < 5 and converged == False and mof.calc.scf_converged == True:
						mof = read('OUTCAR')
						mof, abs_magmoms = continue_magmoms(mof,'INCAR')
						choose_vasp_version(gpt_version,nprocs,calc_swaps)
						mof, calc_swaps = mof_run(mof,calcs_ads(1.5),cif_file,gpt_version,nprocs,calc_swaps)
						if mof == None:
							break
						converged = mof.calc.converged
						loop_i += 1
				if mof != None and mof.calc.scf_converged == True and mof.calc.converged == True:
					write_success(refcode,spin_level,acc_level,vasp_files,cif_file)
				else:
					write_errors(refcode,spin_level,acc_level,vasp_files,cif_file)
					if mof == None:
						pprint('^ VASP crashed')
					elif mof.calc.scf_converged == False:
						pprint('^ SCF did not converge')
					elif mof.calc.converged == False:
						pprint('^ Convergence not reached')
			elif os.path.isfile(outcar_paths[run_i]) == True:
				pprint('COMPLETED: '+spin_level+', '+acc_level)
			mof, run_i, skip_spin2 = prep_next_run(acc_level,run_i,refcode,spin_level)
			if mof == None:
				pprint('Skipping rest because of errors')
				break
			if spin_level == 'spin2':
				mag_indices = get_mag_indices(mof)
				old_mof = read(spin1_final_mof_path)
				if np.sum(np.abs(mof.get_initial_magnetic_moments()[mag_indices] - old_mof.get_magnetic_moments()[mag_indices]) >= 0.05) == 0:
					pprint('Skipping rest because SPIN2 converged to SPIN1')
					continue

			#***********ISIF 2 (medacc)************
			acc_level = acc_levels[run_i]
			if os.path.isfile(outcar_paths[run_i-1]) == True and os.path.isfile(outcar_paths[run_i]) != True and os.path.isfile(error_outcar_paths[run_i]) != True:
				gpt_version, nprocs = get_gpt_version(kpts_hi,len(mof),nprocs,ppn)
				pprint('Running '+spin_level+', '+acc_level)
				if sum(kpts_lo) == 3 and sum(kpts_hi) > 3:
					clean_files(['CHGCAR','WAVECAR'])
				else:
					manage_restart_files(results_partial_paths[run_i-1]+'/'+spin_level)
				choose_vasp_version(gpt_version,nprocs,calc_swaps)
				mof,calc_swaps = mof_run(mof,calcs_ads(run_i),cif_file,gpt_version,nprocs,calc_swaps)
				if mof != None and mof.calc.scf_converged == True and mof.calc.converged == True:
					write_success(refcode,spin_level,acc_level,vasp_files,cif_file)
				else:
					write_errors(refcode,spin_level,acc_level,vasp_files,cif_file)
					if mof == None:
						pprint('^ VASP crashed')
					elif mof.calc.scf_converged == False:
						pprint('^ SCF did not converge')
					elif mof.calc.converged == False:
						pprint('^ Convergence not reached')
			elif os.path.isfile(outcar_paths[run_i]) == True:
				pprint('COMPLETED: '+spin_level+', '+acc_level)
			mof, run_i, skip_spin2 = prep_next_run(acc_level,run_i,refcode,spin_level)
			if mof == None:
				pprint('Skipping rest because of errors')
				break

			#***********ISIF 2 (final)************
			acc_level = acc_levels[run_i]
			if os.path.isfile(outcar_paths[run_i-1]) == True and os.path.isfile(outcar_paths[run_i]) != True and os.path.isfile(error_outcar_paths[run_i]) != True:
				gpt_version, nprocs = get_gpt_version(kpts_hi,len(mof),nprocs,ppn)
				choose_vasp_version(gpt_version,nprocs,calc_swaps)
				manage_restart_files(results_partial_paths[run_i-1]+'/'+spin_level)
				pprint('Running '+spin_level+', '+acc_level)
				mof,calc_swaps = mof_run(mof,calcs_ads(run_i),cif_file,gpt_version,nprocs,calc_swaps)
				if mof != None and mof.calc.scf_converged == True and mof.calc.converged == True:
					if 'large_supercell' in calc_swaps:
						pprint('Running '+spin_level+', '+acc_level+' (LREAL=False)')
						calc_swaps.remove('large_supercell')
						mof = read('OUTCAR')
						mof, abs_magmoms = continue_magmoms(mof,'INCAR')
						mof, calc_swaps = mof_run(mof,calcs_ads(run_i),cif_file,gpt_version,nprocs,calc_swaps)
						if mof != None and mof.calc.scf_converged == True and mof.calc.converged == True:
							write_success(refcode,spin_level,acc_level,vasp_files,cif_file)
						else:
							write_errors(refcode,spin_level,acc_level,vasp_files,cif_file)
							if mof == None:
								pprint('^ VASP crashed')
							elif mof.calc.scf_converged == False:
								pprint('^ SCF did not converge')
							elif mof.calc.converged == False:
								pprint('^ Convergence not reached')
					else:
						write_success(refcode,spin_level,acc_level,vasp_files,cif_file)
				else:
					write_errors(refcode,spin_level,acc_level,vasp_files,cif_file)
					if mof == None:
						pprint('^ VASP crashed')
					elif mof.calc.scf_converged == False:
						pprint('^ SCF did not converge')
					elif mof.calc.converged == False:
						pprint('^ Convergence not reached')
			elif os.path.isfile(outcar_paths[run_i]) == True:
				pprint('COMPLETED: '+spin_level+', '+acc_level)
			mof, run_i, skip_spin2 = prep_next_run(acc_level,run_i,refcode,spin_level)
			if mof == None:
				pprint('Skipping rest because of errors')
				break

			#***********FINAL SPE***********
			acc_level = acc_levels[run_i]
			calc_swaps = []
			if os.path.isfile(outcar_paths[run_i-1]) == True and os.path.isfile(outcar_paths[run_i]) != True and os.path.isfile(error_outcar_paths[run_i]) != True:
				gpt_version, nprocs = get_gpt_version(kpts_hi,len(mof),nprocs,ppn)
				choose_vasp_version(gpt_version,nprocs,calc_swaps)
				manage_restart_files(results_partial_paths[run_i-1]+'/'+spin_level)
				pprint('Running '+spin_level+', '+acc_level)
				mof,calc_swaps = mof_run(mof,calcs_ads(run_i),cif_file,gpt_version,nprocs,calc_swaps)
				if mof != None and mof.calc.scf_converged == True:
					write_success(refcode,spin_level,acc_level,vasp_files,cif_file)
				else:
					write_errors(refcode,spin_level,acc_level,vasp_files,cif_file)
					if mof == None:
						pprint('^ VASP crashed')
					elif mof.calc.scf_converged == False:
						pprint('^ SCF did not converge')
			elif os.path.isfile(outcar_paths[run_i]) == True:
				pprint('COMPLETED: '+spin_level+', '+acc_level)
			mof, run_i, skip_spin2 = prep_next_run(acc_level,run_i,refcode,spin_level)
			if mof == None:
				pprint('Skipping rest because of errors')
				break
			success_path = basepath+'results/'+refcode+'/'+acc_level+'/'+spin_level
			f_final = np.max(np.linalg.norm(read(success_path+'/OUTCAR').get_forces(),axis=1))
			if f_final > f_tol:
				pprint('WARNING: |F_max| = '+str(f_final))
				
			#***********SAVE and CONTINUE***********
			if os.path.isfile(outcar_paths[-1]) == True:
				write_energy(refcode,acc_level,spin_level)
			if skip_spin2 == True:
				pprint('Skipping '+spin_levels[1]+' run')
				break

def run_vol_screen(cif_files):
#Run high-throughput screening

	#Files, spin levels, and accuracy levels to iterate over
	vasp_files = ['INCAR','POSCAR','KPOINTS','POTCAR','OUTCAR',
	'CONTCAR','CHGCAR','AECCAR0','AECCAR2','WAVECAR','opt.traj']
	spin_levels = ['spin1','spin2']
	acc_levels = ['scf_test','isif2','isif3_lowacc','isif3_highacc','final','final_spe']
	f_tol = np.abs(defaults['ediffg'])

	#for each CIF file, optimize the structure
	for cif_file in cif_files:

		nprocs, ppn = get_nprocs()

		refcode = cif_file.split('.cif')[0]
		pprint('***STARTING '+refcode+'***')

		#Make sure MOF isn't running on other process
		working_cif_path = basepath+'working/'+cif_file
		if os.path.isfile(working_cif_path) == True:
			pprint('SKIPPED: Running on another process')
			continue

		#Partial paths to write the OUTCARs
		results_partial_paths = []
		error_outcar_partial_paths = []
		for acc_level in acc_levels:
			results_partial_paths.append(basepath+'results/'+refcode+'/'+acc_level)
			error_outcar_partial_paths.append(basepath+'errors/'+refcode+'/'+acc_level)
		spin1_final_mof_path = results_partial_paths[-1]+'/'+spin_levels[0]+'/OUTCAR'

		#Get the kpoints
		kpts_lo, gamma = get_kpts(cif_file,defaults['kppa_lo'])
		kpts_hi, gamma = get_kpts(cif_file,defaults['kppa_hi'])
		defaults['gamma'] = gamma
		defaults['kpts_lo'] = kpts_lo
		defaults['kpts_hi'] = kpts_hi

		#for each spin level, optimize the structure
		for spin_level in spin_levels:

			#***********PREP FOR RUN***********
			calc_swaps = []
			outcar_paths = []
			error_outcar_paths = []
			run_i = 0
			clean_files(vasp_files)
			for results_partial_path in results_partial_paths:
				outcar_paths.append(results_partial_path+'/'+spin_level+'/OUTCAR')
			for error_outcar_partial_path in error_outcar_partial_paths:
				error_outcar_paths.append(error_outcar_partial_path+'/'+spin_level+'/OUTCAR')

			#***********SCF TEST************
			acc_level = acc_levels[run_i]
			if os.path.isfile(outcar_paths[run_i]) != True and os.path.isfile(error_outcar_paths[run_i]) != True:
				if os.path.isfile(spin1_final_mof_path):
					mof = read(spin1_final_mof_path)
				else:
					mof = cif_to_mof(cif_file)
				mof = set_initial_magmoms(mof,spin_level)
				choose_vasp_version(kpts_lo,len(mof),nprocs,ppn)
				pprint('Running '+spin_level+', '+acc_level)
				mof, calc_swaps = mof_run(mof,calcs_vol(run_i),cif_file,calc_swaps)
				if mof != None:
					write_success(refcode,spin_level,acc_level,vasp_files,cif_file)
				else:
					pprint('^ VASP crashed')
					write_errors(refcode,spin_level,acc_level,vasp_files,cif_file)
			elif os.path.isfile(outcar_paths[run_i]) == True:
				pprint('COMPLETED: '+spin_level+', '+acc_level)
			mof, run_i, skip_spin2 = prep_next_run(acc_level,run_i,refcode,spin_level)
			if mof == None:
				pprint('Skipping rest because of errors')
				break
			warnings = get_warning_msgs(outcar_paths[run_i-1])
			calc_swaps.extend(warnings)

			#***********ISIF 2 (initial)************
			acc_level = acc_levels[run_i]
			if os.path.isfile(outcar_paths[run_i-1]) == True and os.path.isfile(outcar_paths[run_i]) != True and os.path.isfile(error_outcar_paths[run_i]) != True:
				if os.path.isfile(spin1_final_mof_path):
					mof = read(spin1_final_mof_path)
				else:
					mof = cif_to_mof(cif_file)
				manage_restart_files(results_partial_paths[run_i-1]+'/'+spin_level)
				mof = set_initial_magmoms(mof,spin_level)
				choose_vasp_version(kpts_lo,len(mof),nprocs,ppn)
				pprint('Running '+spin_level+', '+acc_level)
				steps = 100
				fmax = 5.0
				mof, dyn, calc_swaps = mof_bfgs_run(mof,calcs_vol(run_i),cif_file,calc_swaps,steps,fmax)
				if mof != None and dyn and mof.calc.scf_converged == True:
					mof = read('OUTCAR')
					mof, abs_magmoms = continue_magmoms(mof,'INCAR')
					mof, calc_swaps = mof_run(mof,calcs_vol(1.5),cif_file,calc_swaps)
					if mof != None and mof.calc.converged  == False and mof.calc.scf_converged == True:
						calc_swaps.append('nsw=90')
						if 'fire' not in calc_swaps:
							calc_swaps.append('fire')
						mof = read('OUTCAR')
						mof, abs_magmoms = continue_magmoms(mof,'INCAR')
						mof, calc_swaps = mof_run(mof,calcs_vol(1.5),cif_file,calc_swaps)
						calc_swaps.remove('nsw=90')
						if 'fire' in calc_swaps:
							calc_swaps.remove('fire')
				if mof != None and mof.calc.scf_converged == True and mof.calc.converged == True:
					write_success(refcode,spin_level,acc_level,vasp_files,cif_file)
				else:
					write_errors(refcode,spin_level,acc_level,vasp_files,cif_file)
					if mof == None:
						pprint('^ VASP crashed')
					elif mof.calc.scf_converged == False:
						pprint('^ SCF did not converge')
					elif mof.calc.converged == False:
						pprint('^ Convergence not reached')
			elif os.path.isfile(outcar_paths[run_i]) == True:
				pprint('COMPLETED: '+spin_level+', '+acc_level)
			mof, run_i, skip_spin2 = prep_next_run(acc_level,run_i,refcode,spin_level)
			if mof == None:
				pprint('Skipping rest because of errors')
				break
			if spin_level == 'spin2':
				mag_indices = get_mag_indices(mof)
				old_mof = read(spin1_final_mof_path)
				if np.sum(np.abs(mof.get_initial_magnetic_moments()[mag_indices] - old_mof.get_magnetic_moments()[mag_indices]) >= 0.05) == 0:
					pprint('Skipping rest because SPIN2 converged to SPIN1')
					continue

			#***********LOW ACCURACY ISIF3***********
			acc_level = acc_levels[run_i]
			if os.path.isfile(outcar_paths[run_i-1]) == True and os.path.isfile(outcar_paths[run_i]) != True and os.path.isfile(error_outcar_paths[run_i]) != True:
				converged = False
				loop_i = 0
				n_runs = 11
				choose_vasp_version(kpts_lo,len(mof),nprocs,ppn)
				manage_restart_files(results_partial_paths[run_i-1]+'/'+spin_level)
				while converged == False and loop_i < n_runs:
					pprint('Running '+spin_level+', '+acc_level+': iteration '+str(loop_i)+'/'+str(n_runs-1))
					if loop_i == n_runs - 1 and 'fire' not in calc_swaps:
						calc_swaps.append('fire')
					mof,calc_swaps = mof_run(mof,calcs_vol(run_i),cif_file,calc_swaps)
					if mof == None:
						break
					converged = mof.calc.converged
					mof = read('OUTCAR')
					mof, abs_magmoms = continue_magmoms(mof,'INCAR')
					loop_i += 1
				if 'fire' in calc_swaps:
					calc_swaps.remove('fire')
				if mof != None and converged == True:
					write_success(refcode,spin_level,acc_level,vasp_files,cif_file)
				else:
					write_errors(refcode,spin_level,acc_level,vasp_files,cif_file)
					if mof == None:
						pprint('^ VASP crashed')
					elif converged == False:
						pprint('^ Convergence not reached')
			elif os.path.isfile(outcar_paths[run_i]) == True:
				pprint('COMPLETED: '+spin_level+', '+acc_level)
			mof, run_i, skip_spin2 = prep_next_run(acc_level,run_i,refcode,spin_level)
			if mof == None:
				pprint('Skipping rest because of errors')
				break

			#***********HIGH ACCURACY ISIF3***********
			acc_level = acc_levels[run_i]
			if os.path.isfile(outcar_paths[run_i-1]) == True and os.path.isfile(outcar_paths[run_i]) != True and os.path.isfile(error_outcar_paths[run_i]) != True:
				converged = False
				loop_i = 0
				n_runs = 11
				V_diff = np.inf
				V_cut = 0.01
				V0 = mof.get_volume()
				choose_vasp_version(kpts_hi,len(mof),nprocs,ppn)
				if sum(kpts_lo) == 3 and sum(kpts_hi) > 3:
					clean_files(['CHGCAR','WAVECAR'])
				else:
					manage_restart_files(results_partial_paths[run_i-1]+'/'+spin_level)
				while (converged == False or V_diff > V_cut) and loop_i < n_runs:
					pprint('Running '+spin_level+', '+acc_level+': iteration '+str(loop_i)+'/'+str(n_runs-1))
					if loop_i == n_runs - 1 and 'fire' not in calc_swaps:
						calc_swaps.append('fire')
					mof,calc_swaps = mof_run(mof,calcs_vol(run_i),cif_file,calc_swaps)
					if mof == None:
						break
					if loop_i > 0:
						converged = mof.calc.converged
					mof = read('OUTCAR')
					V = mof.get_volume()
					mof, abs_magmoms = continue_magmoms(mof,'INCAR')					
					if loop_i > 0:
						V_diff = np.abs((V-V0))/V0
					V0 = V
					loop_i += 1
				if mof != None and converged == True and V_diff <= V_cut and 'large_supercell' in calc_swaps:
					pprint('Running '+spin_level+', '+acc_level+' (LREAL=False)')
					calc_swaps.append('nsw=100')
					calc_swaps.remove('large_supercell')
					mof,calc_swaps = mof_run(mof,calcs_vol(run_i),cif_file,calc_swaps)
					calc_swaps.remove('nsw=100')
					if mof != None and mof.calc.converged == True:
						if 'fire' in calc_swaps:
							calc_swaps.remove('fire')
						write_success(refcode,spin_level,acc_level,vasp_files,cif_file)
					else:
						write_errors(refcode,spin_level,acc_level,vasp_files,cif_file)
						if mof == None:
							pprint('^ VASP crashed')
						elif mof.calc.converged == False:
							pprint('^ Convergence not reached')					
				else:
					write_errors(refcode,spin_level,acc_level,vasp_files,cif_file)
					if mof == None:
						pprint('^ VASP crashed')
					elif converged == False:
						pprint('^ Convergence not reached')
					elif V_diff > V_cut:
						pprint('^ Change in V of '+str(V_diff)+' percent')
			elif os.path.isfile(outcar_paths[run_i]) == True:
				pprint('COMPLETED: '+spin_level+', '+acc_level)
			mof, run_i, skip_spin2 = prep_next_run(acc_level,run_i,refcode,spin_level)
			if mof == None:
				pprint('Skipping rest because of errors')
				break

			#***********FINAL OPT***********
			acc_level = acc_levels[run_i]
			if os.path.isfile(outcar_paths[run_i-1]) == True and os.path.isfile(outcar_paths[run_i]) != True and os.path.isfile(error_outcar_paths[run_i]) != True:
				choose_vasp_version(kpts_hi,len(mof),nprocs,ppn)
				manage_restart_files(results_partial_paths[run_i-1]+'/'+spin_level)
				pprint('Running '+spin_level+', '+acc_level)
				if 'large_supercell' in calc_swaps:
					calc_swaps.remove('large_supercell')
				mof,calc_swaps = mof_run(mof,calcs_vol(run_i),cif_file,calc_swaps)
				if mof == None:
					break
				if mof != None and mof.calc.converged == True and mof.calc.scf_converged == True:
					write_success(refcode,spin_level,acc_level,vasp_files,cif_file)
				else:
					write_errors(refcode,spin_level,acc_level,vasp_files,cif_file)
					if mof == None:
						pprint('^ VASP crashed')
					elif mof.calc.scf_converged == False:
						pprint('^ SCF did not converge')
					elif mof.calc.converged == False:
						pprint('^ Convergence not reached')
			elif os.path.isfile(outcar_paths[run_i]) == True:
				pprint('COMPLETED: '+spin_level+', '+acc_level)
			mof, run_i, skip_spin2 = prep_next_run(acc_level,run_i,refcode,spin_level)
			if mof == None:
				pprint('Skipping rest because of errors')
				break

			#***********FINAL SPE***********
			acc_level = acc_levels[run_i]
			if os.path.isfile(outcar_paths[run_i-1]) == True and os.path.isfile(outcar_paths[run_i]) != True and os.path.isfile(error_outcar_paths[run_i]) != True:
				choose_vasp_version(kpts_hi,len(mof),nprocs,ppn)
				manage_restart_files(results_partial_paths[run_i-1]+'/'+spin_level)
				pprint('Running '+spin_level+', '+acc_level)
				if 'large_supercell' in calc_swaps:
					calc_swaps.remove('large_supercell')
				mof,calc_swaps = mof_run(mof,calcs_vol(run_i),cif_file,calc_swaps)
				if mof == None:
					break
				if mof != None and mof.calc.scf_converged == True:
					write_success(refcode,spin_level,acc_level,vasp_files,cif_file)
				else:
					write_errors(refcode,spin_level,acc_level,vasp_files,cif_file)
					if mof == None:
						pprint('^ VASP crashed')
					elif mof.calc.scf_converged == False:
						pprint('^ SCF did not converge')
			elif os.path.isfile(outcar_paths[run_i]) == True:
				pprint('COMPLETED: '+spin_level+', '+acc_level)
			mof, run_i, skip_spin2 = prep_next_run(acc_level,run_i,refcode,spin_level)
			if mof == None:
				pprint('Skipping rest because of errors')
				break
			success_path = basepath+'results/'+refcode+'/'+acc_level+'/'+spin_level
			f_final = np.max(np.linalg.norm(read(success_path+'/OUTCAR').get_forces(),axis=1))
			if f_final > f_tol:
				pprint('WARNING: |F_max| = '+str(f_final))

			#***********SAVE and CONTINUE***********
			if os.path.isfile(outcar_paths[-1]) == True:
				write_energy(refcode,acc_level,spin_level)
			if skip_spin2 == True:
				pprint('Skipping '+spin_levels[1]+' run')
				break

