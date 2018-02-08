from cif_handler import get_cif_files
from screen import run_screen
from janitor import prep_paths
prep_paths()
cif_files = get_cif_files()
run_screen(cif_files)