#-------------Paths for files-------------
coremof_path = 'C:/Users/asros/Desktop/1/' #path to MOF CIFs to oxygenate
newmofs_path = 'C:/Users/asros/Desktop/2/' #path to store oxygenated CIFs
omsdata_path = 'C:/Users/asros/Desktop/1/OMS_data/' #path to Zeo++ .omsex and .oms files
error_path = newmofs_path+'errors/' #path to store failed CIFs

#-------------Parameters-------------
guess_length = 2.0 #M-adsorbate bond distance
ads_species = 'O' #adsorbate species
rcut = 2.5 #cutoff for calculating nearest neighbors
sum_cutoff = 0.5 #sum(r_i) cutoff for planar structures
rmse_tol = 0.25 #RMSE cutoff for planar structures
overlap_tol = 1.3 #distance tolerance for overlapping atoms