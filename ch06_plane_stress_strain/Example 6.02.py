import sys
from pathlib import Path

sys.path.append(str(Path(__file__).resolve().parent.parent))

import numpy as np
from fem.elements.cst import GradientMatrix, elementArea, principalStresses
from fem.constitutiveMatrix import planeStress

# material properties
E = 210e9
nu = 0.30

# geometry properties
t = 20e-3 # thickness

# mesh
nodeCoordinates = np.array([
    [  0,   0],
    [  0, 200],
    [400, 200],
    [400,   0]
])*1e-3

elementNodes = np.array([
    [0, 2, 1],
    [0, 3, 2]
])

GDOF = 2*len(nodeCoordinates)
K = np.zeros((GDOF, GDOF))
displacements = np.zeros(GDOF)

# stiffness matrix assembly
for element in elementNodes:
    nodes = nodeCoordinates[element]

    A2 = elementArea(nodes)
    B = GradientMatrix(nodes)    
    
    D = planeStress(E, nu)

    DOF = [2*element[0], 2*element[0] + 1, 2*element[1], 2*element[1] + 1, 2*element[2], 2*element[2] + 1]
    
    K[np.ix_(DOF, DOF)] += t*A2/2*B.T @ D @ B

# boundary conditions
bc = [0, 1, 2, 3]
mask = np.arange(GDOF)
mask = [False if x in bc else True for x in mask]

# load 
P = 14e3
load = [4, 6]
F = np.zeros(GDOF)
F[load] = P

# solution
displacements[mask] = np.linalg.solve(K[np.ix_(mask, mask)], 
                                       F[mask])
# stress calculation
n = 0
for element in elementNodes:
    n += 1
    nodes = nodeCoordinates[element]

    B = GradientMatrix(nodes)    
    
    D = planeStress(E, nu)

    DOF = [2*element[0], 2*element[0] + 1, 2*element[1], 2*element[1] + 1, 2*element[2], 2*element[2] + 1]
    elementDisplacement = displacements[DOF]

    sigma = D @ B @ elementDisplacement
    
    print(f'\nPrincipal stresses at Element {n}\n')
    principalStresses(sigma)