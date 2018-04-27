from pymofscreen.cif_handler import get_cif_files
from pymofscreen.screen import screener

#Set up paths
mofpath = '../example_structures/'
basepath = '../'
submit_script = 'sub_screen.job'

#Get CIF files
cif_files = get_cif_files(mofpath)

#Construct screener object
s = screener(mofpath,basepath,submit_script=submit_script)

#Run screening
for cif_file in cif_files:
	mof = s.run_screen(cif_file,'volume')