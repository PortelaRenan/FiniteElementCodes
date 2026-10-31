import sys
from pathlib import Path

sys.path.append(str(Path(__file__).resolve().parent.parent))

import numpy as np
from fem.elements.cst import GradientMatrix, elementArea, principalStresses
from fem.constitutiveMatrix import planeStress

# material properties
E = 210e9
nu = 0.25

# geometry properties
t = 1e-2 # thickness

# mesh
nodeCoordinates = np.array([
    [  0, -1],
    [  2,  0],
    [  0,  1]
])*0.0254

elementNodes = np.array([
    [0, 1, 2]
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

displacements = np.array([
    0, 0.0625, 0.03, 0, 0, 0.0625
])*1e-3

# stress calculation
B = GradientMatrix(nodeCoordinates)    

D = planeStress(E, nu)

DOF = [2*element[0], 2*element[0] + 1, 2*element[1], 2*element[1] + 1, 2*element[2], 2*element[2] + 1]
elementDisplacement = displacements[DOF]

sigma = D @ B @ elementDisplacement

print(f'\nPrincipal stresses 6.4(a) \n')

principalStresses(sigma)

# mesh
nodeCoordinates = np.array([
    [  1.2, 0],
    [  2.4, 0],
    [  1.2, 1]
])*0.0254

# stress calculation
B = GradientMatrix(nodeCoordinates)    

D = planeStress(E, nu)

DOF = [2*element[0], 2*element[0] + 1, 2*element[1], 2*element[1] + 1, 2*element[2], 2*element[2] + 1]
elementDisplacement = displacements[DOF]

sigma = D @ B @ elementDisplacement

print(f'\nPrincipal stresses 6.4(b) \n')

principalStresses(sigma)

# mesh
nodeCoordinates = np.array([
    [  0, 0],
    [  2, 0],
    [  0, 1]
])*0.0254

# stress calculation
B = GradientMatrix(nodeCoordinates)    

D = planeStress(E, nu)

DOF = [2*element[0], 2*element[0] + 1, 2*element[1], 2*element[1] + 1, 2*element[2], 2*element[2] + 1]
elementDisplacement = displacements[DOF]

sigma = D @ B @ elementDisplacement

print(f'\nPrincipal stresses 6.4(c) \n')

principalStresses(sigma)