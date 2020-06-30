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

    ### creating the toeplitz matrix from the fft-permittivity
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

    pass


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
    pass