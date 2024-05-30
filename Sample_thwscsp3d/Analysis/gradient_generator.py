import numpy as np
from tqdm import tqdm 
from findiff import FinDiff, coefficients, Coefficient
#--------------------------------------------------------------------

with tqdm(desc="Loading data", unit=" rows") as pbar:
    data = np.loadtxt('ac74.dat') # input the file that is generated from fortran code
    pbar.update(data.shape[0])

print('Data loaded')

with tqdm(desc="Processing data", total=4) as pbar:
    data2 = data[:,3]
    F = data2.reshape(160,288,511)
    r = data[:511,1]
    pbar.update(1)

    phi = data[:,2]
    phi2 = phi.reshape(160,288,511)
    p = phi2[0,:,0]
    pbar.update(1)

    theta = data[:,0]
    theta2 = theta.reshape(160,288,511)
    t = theta2[:,0,0]
    pbar.update(1)

    phi_,theta_,r_ = np.meshgrid(p,t,r)
    X = r*np.sin(theta_)*np.cos(phi_)
    Y = r*np.sin(theta_)*np.sin(phi_)
    Z = r*np.cos(theta_)
    pbar.update(1)

print('Data processed')
np.save('/home/ijas/Fortran_IIA/testing_code/res/Components2/phi.npy', phi_)
np.save('/home/ijas/Fortran_IIA/testing_code/res/Components2/theta.npy', theta_)
np.save('/home/ijas/Fortran_IIA/testing_code/res/Components2/r.npy', r_)

#--------------------------------------------------------------------

with tqdm(desc="Calculating Gradients",total=3) as pbar:
    dr = r[2]-r[1]
    d_dr = FinDiff(2, dr, 1, acc=4)
    diff_r = d_dr(F)
    pbar.update(1)
    #print('grad_r_max',np.max(diff_r))
    #print('grad_r_min',np.min(diff_r))

    dt = t[2]-t[1]
    d_dt = FinDiff(0, dt, 1, acc=4)

    diff_t = d_dt(F)
    pbar.update(1)
    #print('grad_t_max',np.max(diff_t))
    #print('grad_t_min',np.min(diff_t))

    dp = p[2]-p[1]
    d_dp = FinDiff(1, dp, 1, acc=4)

    diff_p = d_dp(F)
    pbar.update(1)
    #print('grad_p_max',np.max(diff_p))
    #print('grad_p_min',np.min(diff_p))
#--------------------------------------------------------------------

Br = -diff_r # give the location where it is needed to be saved 
np.save('/home/ijas/Fortran_IIA/testing_code/res/Components2/74Br.npy',Br)
Bt = -(1/r_)*diff_t
np.save('/home/ijas/Fortran_IIA/testing_code/res/Components2/74Bt.npy',Bt)
Bp = -(1/(r_*np.sin(theta_)))*diff_p
np.save('/home/ijas/Fortran_IIA/testing_code/res/Components2/74Bp.npy',Bp)
print('dimension of Grid :', phi_.shape)
print('dimension of B :', Br.shape)
print('max of Br: ', np.max(Br), 'min of Br:', np.min(Br))
#print('Total Magnetic Field :', np.sqrt((Br)**2 + (Bt)**2 + (Bp)**2))
