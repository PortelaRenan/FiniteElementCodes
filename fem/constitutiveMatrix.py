import numpy as np

def planeStress(E, nu):
    return E/(1 - nu**2) * np.array([[1,  nu, 0],
                                     [nu, 1,  0],
                                     [0,  0, (1 - nu)/2]])

def planeStrain(E, nu):
    return E/((1 + nu)*(1 - 2*nu)) * np.array([[1 - nu,  nu, 0],
                                               [nu, 1 - nu,  0],
                                               [0,  0, (1 - 2*nu)/2]])

def isotropicMaterial(E, thickness, poisson):
    D11 = E*thickness**3/12/(1-poisson**2)
    D22 = D11; 
    D12 = poisson*D11
    D66 = (1-poisson)*D11/2

    return np.array([[D11, D12, 0],
                     [D12, D22, 0],
                     [0, 0, D66]])

def orthotropicMaterial(E1, E2, G12, poisson12, thickness):
    poisson21 = poisson12*E2/E1
    I = thickness**3/12
    D11 = E1*I/(1-poisson12*poisson21)
    D22 = E2*I/(1-poisson12*poisson21)
    D12 = poisson12*D22
    D66 = G12*I

    return np.array([[D11, D12, 0],
                     [D12, D22, 0],
                     [0, 0, D66]])