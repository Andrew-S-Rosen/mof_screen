from pymofscreen.cif_handler import get_cif_files
from pymofscreen.screen import screener

#Set up paths
mofpath = 'folder_of_CIF_files'
basepath = 'folder_to_store_data'
submit_script = 'submit_screen.sh'

#Get CIF files
cif_files = get_cif_files(mofpath)

#Construct screener object
s = screener(mofpath,basepath,submit_script=submit_script)

#Run screening
for cif_file in cif_files:
	mof = s.run_screen(cif_file,'volume')