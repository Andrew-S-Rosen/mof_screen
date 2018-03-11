import os
import numpy as np
from ase.io import read
from settings import coremof_path, OMS_ads_species
from add_adsorbate import write_ads_file
from ads_sites import get_opt_ads_site
from geom import get_vire_NNs
from path_prep import prep_paths, get_refcode

prep_paths()
for struct_file in os.listdir(coremof_path):
	if not os.path.isfile(coremof_path+struct_file):
		continue
	refcode = get_refcode(struct_file)
	mof = read(coremof_path+struct_file)
	center_idx = [atom.index for atom in mof if atom.symbol == OMS_ads_species][-1]
	neighbors_idx = get_vire_NNs(struct_file,center_idx,OMS_ads_species)
	cnum = len(neighbors_idx)
	center_coord = mof[center_idx].position
	mic_coords = np.squeeze(mof.get_distances(center_idx,neighbors_idx,mic=True,vector=True))
	ads_site = get_opt_ads_site(struct_file,cnum,mic_coords,center_idx,center_coord)
	write_ads_file(refcode,struct_file,ads_site)