#PATHS
kpts_path = '/projects/p30148/vasp_jobs/MOFs/reoptimized_core1/results_cleaned/data_for_phase2/kpts.txt'
mofpath = '/projects/p30148/vasp_jobs/MOFs/reoptimized_core1/results_cleaned/data_for_phase2/cifs/'
basepath = '/projects/p30148/vasp_jobs/MOFs/oxidized_oms/'

#CALCULATION DETAILS
ads_species = 'O' #adsorbate species
skip_mofs = [] #MOFs to skip in analysis
phase = 2 #phase of screening

#FILE NAMES
submit_script = 'sub_screen2_temp.job'
stdout_file = 'driver.out'