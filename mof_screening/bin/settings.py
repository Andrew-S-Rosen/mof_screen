#PATHS
kpts_path = '/projects/p30148/vasp_jobs/MOFs/reoptimized_core1/results_cleaned/data_for_phase2/kpts.txt'
mofpath = '/projects/p30148/vasp_jobs/MOFs/reoptimized_core1/results_cleaned/data_for_phase2/cifs/'
basepath = '/projects/p30148/vasp_jobs/MOFs/oxidized_oms/'
module_paths = ''

#CALCULATION DETAILS
ads_species = 'O' #adsorbate species
skip_mofs = [] #MOFs to skip in analysis
phase = 2 #phase of screening

#EXECUTABLES
vasp_module = 'vasp/5.4.1'
vtst_module = 'vasp/5.4.1_vtst'
vasp_ex = ['vasp_std','vasp_gam']
vtst_ex = ['vasp_std_vtst','vasp_gam_vtst']
parallel_cmd = 'mpirun -n'

#FILE NAMES
submit_script = 'sub_screen2_temp.job'
stdout_file = 'driver.out'