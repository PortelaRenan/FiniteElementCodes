import sys
from pathlib import Path

sys.path.append(str(Path(__file__).resolve().parent.parent))

import numpy as np
from fem.elements.cst import GradientMatrix, elementArea, principalStresses
from fem.constitutiveMatrix import planeStress
from fem.plotting import plotCST
from fem.gauss import gaussPoints

# material properties
E = 105e9
nu = 0.30

# geometry properties
t = 5e-3 # thickness

# ---------- 6.11 (a) ----------

# mesh
nodeCoordinates = np.array([
    [   0,   0],
    [ 100,   0],
    [ 100, 150]
])*1e-3

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

# boundary conditions
bc = [0, 1, 2, 3]
mask = np.arange(GDOF)
mask = [False if x in bc else True for x in mask]

# load 
Px = 10e6
L = np.sqrt( (nodeCoordinates[0][0] - nodeCoordinates[1][0])**2 + (nodeCoordinates[1][1] - nodeCoordinates[0][1])**2 )
F = np.array([0, 0, 1, 0, 2, 0])*t*Px*L/6


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

# ---------- 6.14 (b) ----------

# mesh
nodeCoordinates = np.array([
    [   0,   0],
    [ 100,   0],
    [ 100, 150]
])*1e-3

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

# boundary conditions
bc = [0, 1, 2, 3]
mask = np.arange(GDOF)
mask = [False if x in bc else True for x in mask]

# load 
P = 10e6
L = np.linalg.norm(
    np.array(nodeCoordinates[1]) -
    np.array(nodeCoordinates[0])
)
load = [2, 4]
F = np.zeros(GDOF)
p, w = gaussPoints(n = 2) 
J = L/2

for pi, wi in zip(p, w):
    
    # Shape functions
    N2 = (1 - pi)/2
    N3 = (1 + pi)/2

    # Quadratic pressure distribution
    px = P*((1 + pi)/2)**2

    # Equivalent nodal forces
    P2 = N2 * px * wi * J * t
    P3 = N3 * px * wi * J * t

    F[load] += [P2, P3]

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