import numpy as np
from tqdm import tqdm 
from findiff import FinDiff, coefficients, Coefficient
from scipy.integrate import trapezoid 
from scipy.integrate import cumulative_trapezoid
#--------------------------------------------------------------------

with tqdm(desc="Loading data", unit=" rows") as pbar:
    data = np.loadtxt('ac143.dat')
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
print('shape of r',r.shape)
print('shape of t',t.shape)
print('shape of p',p.shape)

#--------------------------------------------------------------------

with tqdm(desc="Calculating Gradients",total=3) as pbar:
    dr = r[1]-r[0]
    d_dr = FinDiff(2, dr, 1, acc=4)
    diff_r = d_dr(F)
    pbar.update(1)
    #print('grad_r_max',np.max(diff_r))
    #print('grad_r_min',np.min(diff_r))

    dt = t[1]-t[0]
    d_dt = FinDiff(0, dt, 1, acc=4)

    diff_t = d_dt(F)
    pbar.update(1)
    #print('grad_t_max',np.max(diff_t))
    #print('grad_t_min',np.min(diff_t))

    dp = p[1]-p[0]
    d_dp = FinDiff(1, dp, 1, acc=4)

    diff_p = d_dp(F)
    pbar.update(1)
    #print('grad_p_max',np.max(diff_p))
    #print('grad_p_min',np.min(diff_p))
#--------------------------------------------------------------------

Br = -diff_r
Bt = -(1/r_)*diff_t
Bp = -(1/(r_*np.sin(theta_)))*diff_p
B = np.sqrt((Br)**2 + (Bt)**2 + (Bp)**2)
dv = dt*dr*dp

print('dimension of Grid :', phi_.shape)
print('dimension of Br :', Br.shape)
print('max of Br: ', np.max(Br), 'min of Br:', np.min(Br))
print('Total Magnetic Field :', np.sum(B), 'Dim of B:', B.shape)

integrand = (B**2) * (r_**2) * (np.sin(theta_)) 
print('dimension of integrand :', integrand.shape)
integrand_r = cumulative_trapezoid(integrand,r,axis = 2)
print('dimension of integrand_r :', integrand_r.shape)
integrand_phi = cumulative_trapezoid(integrand_r[:,:,-1],p,axis = 1)
print('dimension of integrand_phi :', integrand_phi.shape)
integrand_theta = cumulative_trapezoid(integrand_phi[:,-1],t,axis = 0)
print('dimension of integrand_theta :', integrand_theta.shape)


E = integrand_theta[-1]/20
print('Magnetic Energy :', E, 'Dim of E:', E.shape)
#print(np.sum(E)/20)
