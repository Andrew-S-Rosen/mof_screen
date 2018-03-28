#PATHS
kpts_path = '/projects/p30148/vasp_jobs/MOFs/reoptimized_core1/results_cleaned/data_for_phase2/kpts.txt'
mofpath = '/projects/p30148/vasp_jobs/MOFs/oxidized_oms/results_cleaned/data_for_phase3/'
basepath = '/projects/p30148/vasp_jobs/MOFs/phase3/runner/runner1/'

#CALCULATION DETAILS
ads_species = 'O' #adsorbate species
skip_mofs = [] #MOFs to skip in analysis
phase = 3 #phase of screening
f_tol = 0.03 #force tolerance

#FILE NAMES
submit_script = 'sub_screen3_temp.job'
stdout_file = 'driver.out'