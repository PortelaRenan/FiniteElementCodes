import sys
from pathlib import Path

sys.path.append(str(Path(__file__).resolve().parent.parent))

import numpy as np
from fem.elements.cst import GradientMatrix, elementArea, principalStresses
from fem.constitutiveMatrix import planeStress
from fem.plotting import plotCST

# material properties
E = 105e9
nu = 0.30
weightDensity = 154.2 # kN/m³

# geometry properties
t = 5e-3 # thickness

# mesh
nodeCoordinates = np.array([
    [   0,   0],
    [ 400,   0],
    [ 400, 400],
    [   0, 400],
    [ 200, 200]
])*1e-3

elementNodes = np.array([
    [0, 1, 4],
    [1, 2, 4],
    [4, 2, 3],
    [0, 4, 3]
])

GDOF = 2*len(nodeCoordinates)
K = np.zeros((GDOF, GDOF))
displacements = np.zeros(GDOF)

# load 
P = -50e3
load = [9]
F = np.zeros(GDOF)
F[load] = P
g = 9.81 # gravity
bodyForce = np.array([0, -weightDensity*g])

# stiffness matrix assembly
for element in elementNodes:
    nodes = nodeCoordinates[element]

    A2 = elementArea(nodes)
    B = GradientMatrix(nodes)    
    
    D = planeStress(E, nu)

    DOF = [2*element[0], 2*element[0] + 1, 2*element[1], 2*element[1] + 1, 2*element[2], 2*element[2] + 1]
    
    K[np.ix_(DOF, DOF)] += t*A2/2*B.T @ D @ B
    F[DOF] += abs(A2)/2*t/3*np.tile(bodyForce, 3)

# boundary conditions
bc = [4, 5, 6, 7]
mask = np.arange(GDOF)
mask = [False if x in bc else True for x in mask]

# solution
displacements[mask] = np.linalg.solve(K[np.ix_(mask, mask)], 
                                       F[mask])

np.set_printoptions(precision=4)
print(f'Displacement = \n{displacements} mm')

# plot
plotCST(nodeCoordinates, elementNodes, displacements, scale=1e2)

# reaction
reaction = np.matmul(K, displacements) - F
print(f'Reaction = \n{reaction/1e3} kN')

for element in elementNodes:
    nodes = nodeCoordinates[element]
    
    B = GradientMatrix(nodes)    
    
    DOF = [2*element[0], 2*element[0] + 1, 2*element[1], 2*element[1] + 1, 2*element[2], 2*element[2] + 1]

    elementDisplacement = displacements[DOF]

    sigma = D @ B @ elementDisplacement

    print(f'\nPrincipal stresses\n')

    principalStresses(sigma)