#-------------Paths for files-------------
coremof_path = 'C:/Users/asros/Desktop/1/' #path to MOF CIFs to oxygenate
newmofs_path = 'C:/Users/asros/Desktop/2/' #path to store oxygenated CIFs
error_path = newmofs_path+'errors/' #path to store failed CIFs

#-------------Parameters-------------
OMS_ads_species = 'O'
ads_species = 'H'
guess_length = 1.0
rcut = 2.5 #cutoff for calculating nearest neighbors
sum_cutoff = 0.5 #sum(r_i) cutoff for planar structures
rmse_tol = 0.25 #RMSE cutoff for planar structures
overlap_tol = 0.75 #distance tolerance for overlapping atoms
