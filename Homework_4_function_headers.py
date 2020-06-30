'''Homework 4, Computational Photonics, SS 2020:  Fourier modal method.
'''

import numpy as np
from numpy.linalg import eig, solve
from scipy.linalg import toeplitz
from scipy.fftpack import fft


def fmm1d_te_layer_modes(perm, period, k0, kx, N):
    '''Calculates the TE eigenmodes of a one-dimensional grating layer.

    Arguments
    ---------
        perm: 1d-array
            permittivity distribution
        period: float
            grating period
        k0: float
            vacuum wavenumber
        kx: float
            transverse wave vector
        N: int
            number of positive Fourier orders

    Returns
    -------
        beta: 1d-array
            propagation constants of the eigenmodes
        phie: 2d-array
            Fourier coefficients of the eigenmodes (each column
            corresponds to one mode)
    '''
    Nx = perm.size

    m = np.linspace(-N,N,2*N+1) # we don't want all the fourier components, but only 2N+1 of them

    perm_ft = fft(perm)/Nx

    if (perm_ft.size/2 >= N):

        print("ok, there are enough frequencies to fill the Toeplitz matrix")

    else: print("not enough frequencies from fft")

    perm_ft_pos = perm_ft[0:2*N+1]

    perm_ft_neg=np.zeros_like(perm_ft_pos)

    perm_ft_neg[0] = perm_ft_pos[0]

    perm_ft_neg[1:] = perm_ft[-1:-2*N -1:-1]

    perm_toepl = toeplitz(perm_ft_pos, perm_ft_neg)

    if (perm_ft[0] == np.mean(perm)):

        print("ok, scaling is correct")

    else: print("scaling is incorrect")



    ### creating the K^2 matrix

    G = m*2*np.pi / period

    K_sq = (np.diag(G + kx))**2



    ### solving the eigenvalue problem

    beta_sq, phi_e = eig(perm_toepl * k0**2 - K_sq)



    ### homogeneous layer

    if (np.unique(perm).size ==1) or (perm.size == 1):

        print("the layer is homogeneous with eps = ", perm[0])

        phi_e = np.eye(2*N +1)

        beta_sq = np.eye(2*N+1) * perm[0] * k0**2 - K_sq +0j

    else: print("ok, the layer is inhomogeneous")



    ### rooting beta

    beta = np.sqrt(beta_sq)

    beta[(beta.real + beta.imag < 0)] = -beta[(beta.real + beta.imag < 0)]



    return beta, phi_e


def fmm1d_te(lam, theta, period, perm_in, perm_out,
             layer_perm, layer_ticknesses, N):
    '''Calculates the TE diffraction efficiencies for a one-dimensional
    layered grating structure using the T-matrix method.

    Arguments
    ---------
        lam: float
            vacuum wavelength
        theta: float
            angle of incidence in rad
        period: float
            grating period
        perm_in: float
            permittivity on the incidence side
        perm_out: float
            permittivity on the exit side
        layer_perm: 2d-array
            permittivity distribution within the grating
            layers (matrix, each row corresponds to one layer)
        layer_thicknesses: 1d-array
            thicknesses of the grating layers
        N: int
            number of positive Fourier orders

    Returns
    -------
        eta_r: 1d-array
            diffraction efficiencies of the reflected diffraction orders
        eta_t: 1d-array
            diffraction efficiencies of the transmitted diffraction orders
        r: 1d-array
            amplitude reflection coefficients of the reflected
            diffraction orders
        t: 1d-array
            amplitude transmission coefficients of the transmitted
            diffraction orders
    '''
    #parameters
    k0=2*np.pi/lam # wave vector in vacuum
    kx=k0*perm_in*np.sin(theta) # x component of incident wave vector
    K=np.diag(np.linspace(-N,N,2*N+1)*2*np.pi/period+kx+0j) # K matrix
    
    #before the grating(homogeneous medium)
    phie1=np.eye(2*N+1) # eigenvector
    beta_diag1=np.sqrt(k0**2*perm_in*np.eye(2*N+1)-K**2) #eigenvalue
    beta1=beta_0=np.diag(beta_diag1)
    t1=np.block([[phie1,phie1],[phie1@beta_diag1,-phie1@beta_diag1]]) # first matrix for interface
    T0=np.eye(2*(2*N+1)) #stack T-matrix
    
    # iterations for each layer
    for i in range(0,np.size(layer_perm,0)):
        #update thickness
        if i==0:
            z=0
        else:
            z=layer_thicknesses[i-1]
        perm2=layer_perm[i,:] #permittivity of the layer
        beta2, phie2 = fmm1d_te_layer_modes2(perm2, period, k0, kx, N) #calculate eigenvalues and eigenvectors
        beta_diag2=np.diag(beta2)
        t2=np.block([[phie2,phie2],[phie2@beta_diag2,-phie2@beta_diag2]]) # second matrix for interface
        p_plus=np.diag(np.exp(1j*beta1*z))
        p_minus=np.diag(np.exp(-1j*beta1*z))
        p=np.block([[p_plus,np.zeros_like(p_plus)],[np.zeros_like(p_plus),p_minus]])
        
        #update layer T-matrix
        T=np.linalg.solve(t2,t1) #interface T-matrix
        T=T@p
        
        # update stack transfer matrix
        T0=T@T0 
        
        #update necessary parameters for next iteration
        t1=t2
        beta1=beta2
    
    #after the grating(homogeneous medium)
    phie2=np.eye(2*N+1)
    beta_diag2=np.sqrt(k0**2*perm_out*np.eye(2*N+1)-K**2)
    beta2=beta_out=np.diag(beta_diag2)
    
    t2=np.block([[phie2,phie2],[phie2@beta_diag2,-phie2@beta_diag2]])
    z=layer_thicknesses[-1]
    p_plus=np.diag(np.exp(1j*beta1*z))
    p_minus=np.diag(np.exp(-1j*beta1*z))
    p=np.block([[p_plus,np.zeros_like(p_plus)],[np.zeros_like(p_plus),p_minus]])
    T=np.linalg.solve(t2,t1)
    T=T@p #update layer T-matrix
    T0=T@T0 # update stack transfer matrix
    
    #different elements of stack transfer matrix
    t11=T0[0:2*N+1,0:2*N+1]
    t12=T0[0:2*N+1,2*N+1:4*N+2]
    t21=T0[2*N+1:4*N+2,0:2*N+1]
    t22=T0[2*N+1:4*N+2,2*N+1:4*N+2]
    
    #incident field
    ain=np.zeros([2*N+1,1])
    ain[N]=1 

    #calculate transmission and reflection coefficients
    r=-np.linalg.solve(t22,t21)
    t=t11+t12@r
    r=r@ain
    t=t@ain
    
    #calculate efficiencies
    r=r.flatten()
    t=t.flatten()
    beta_in=k0*perm_in*np.cos(theta) # z-component of incident wave vector
    eta_r=beta_0.real/beta_in.real*np.abs(r)**2 # efficiency of transmitted light
    eta_t=beta_out.real/beta_in.real*np.abs(t)**2 # efficiency of reflected light 
    
    return eta_r,eta_t,r,t
