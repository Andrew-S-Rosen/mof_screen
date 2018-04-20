from cif_handler import get_cif_files
from janitor import prep_paths
from screen import run_ads_screen, run_vol_screen
from settings import phase
import time
from random import randint

#Run screening analysis
time.sleep(randint(1,60))
prep_paths()
cif_files = get_cif_files()
if phase == 1:
	run_vol_screen(cif_files)
else:
	run_ads_screen(cif_files)