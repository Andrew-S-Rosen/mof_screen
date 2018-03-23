#PATHS
kpts_path = '/projects/p30148/vasp_jobs/MOFs/reoptimized_core1/results_cleaned/data_for_phase2/kpts.txt'
#mofpath = '/projects/p30148/vasp_jobs/MOFs/oxidized_oms/results_cleaned/data_for_phase3/'
mofpath = '/projects/p30148/vasp_jobs/MOFs/testing/runners/here/'
basepath = '/projects/p30148/vasp_jobs/MOFs/testing/'

#CALCULATION DETAILS
ads_species = 'O' #adsorbate species
skip_mofs = [] #MOFs to skip in analysis
phase = 3 #phase of screening

#FILE NAMES
submit_script = 'sub_screen3_temp.job'
stdout_file = 'driver.out'