from cif_handler import get_cif_files
from janitor import prep_paths
from screen import run_screen

#Run screening analysis
prep_paths()
cif_files = get_cif_files()
run_screen(cif_files)