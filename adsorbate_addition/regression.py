import numpy as np

def OLS_fit(xyz):
#ordinary LS regression to z=a+bx+cy and calculation of normal vector

	x = xyz[:,0][np.newaxis].T
	y = xyz[:,1][np.newaxis].T
	z = xyz[:,2][np.newaxis].T
	onevec = np.ones((len(x),1))
	A = np.hstack((x,y,onevec))
	B = z
	fit = np.squeeze(np.dot(np.linalg.pinv(np.dot(A.T,A)),np.dot(A.T,B)))
	z_fit = fit[0]*x+fit[1]*y+fit[2]
	ss_res = sum((z_fit-z)**2)[0]
	ss_tot = sum((z-np.mean(z))**2)[0]
	r2 = 1-ss_res/ss_tot
	normal_vec = np.array([fit[0],fit[1],-1])
	if r2 < 1 and len(x) == 2:
		raise ValueError('Poor linear fit to two points')
		
	return normal_vec

def TLS_fit(xyz):
#total LS regression to ax+by+cz+d=0 and calculation of normal vector

	xyz_mean = np.mean(xyz,axis=0)
	xyz_sub = xyz-xyz_mean
	[u,s,v] = np.linalg.svd(xyz_sub,full_matrices=False)
	v = v.T
	normal_vec = v[:,-1]
	a = normal_vec[0]
	b = normal_vec[1]
	c = normal_vec[2]
	d = -np.dot(xyz_mean,normal_vec)
	fit = a*xyz[:,0]+b*xyz[:,1]+c*xyz[:,2]+d
	ss_res = sum(fit**2)
	rmse = (ss_res/np.shape(xyz)[0])**0.5

	return rmse, normal_vec
