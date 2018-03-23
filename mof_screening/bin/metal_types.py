import numpy as np
spblock_metals = [3,4,11,12,19,20,37,38,55,56,87,88]
dblock_metals = np.concatenate((np.arange(21,30,1),np.arange(39,48,1),np.arange(71,80,1),np.arange(103,112,1)),axis=0).tolist()
fblock_metals = np.concatenate((np.arange(57,71,1),np.arange(89,103,1)),axis=0).tolist()
nonmetal_list = [1,2,6,7,8,9,10,15,16,17,18,34,35,36,53,54,86]
metal_list = [val for val in np.arange(1,119,1) if val not in nonmetal_list]
mag_list = [metal for metal in metal_list if metal not in spblock_metals]
nomag_list = [val for val in np.arange(1,119,1) if val not in mag_list]
poor_metals = [metal for metal in metal_list if metal not in dblock_metals+fblock_metals+spblock_metals]
