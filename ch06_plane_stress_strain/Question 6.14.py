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

# geometry properties
t = 5e-3 # thickness

# ---------- 6.14 (a) ----------

# mesh
nodeCoordinates = np.array([
    [   0,   0],
    [ 500,   0],
    [ 500, 250],
    [   0, 250],
    [ 250, 125]
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

# stiffness matrix assembly
for element in elementNodes:
    nodes = nodeCoordinates[element]

    A2 = elementArea(nodes)
    B = GradientMatrix(nodes)    
    
    D = planeStress(E, nu)

    DOF = [2*element[0], 2*element[0] + 1, 2*element[1], 2*element[1] + 1, 2*element[2], 2*element[2] + 1]
    
    K[np.ix_(DOF, DOF)] += t*A2/2*B.T @ D @ B

# boundary conditions
bc = [0, 1, 6, 7, 4, 5]
mask = np.arange(GDOF)
mask = [False if x in bc else True for x in mask]

# load 
Px, Py = 40e3*np.cos(30*np.pi/180), -40e3*np.sin(30*np.pi/180)
load = [2, 3]
F = np.zeros(GDOF)
F[load] = Px, Py

print(F)

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
    [ 100, 100],
    [   0, 100],
    [  50,  50]
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

# stiffness matrix assembly
for element in elementNodes:
    nodes = nodeCoordinates[element]

    A2 = elementArea(nodes)
    B = GradientMatrix(nodes)    
    
    D = planeStress(E, nu)

    DOF = [2*element[0], 2*element[0] + 1, 2*element[1], 2*element[1] + 1, 2*element[2], 2*element[2] + 1]
    
    K[np.ix_(DOF, DOF)] += t*A2/2*B.T @ D @ B

# boundary conditions
bc = [0, 1, 6, 7]
mask = np.arange(GDOF)
mask = [False if x in bc else True for x in mask]

# load 
P = 20e3
load = [2, 4]
F = np.zeros(GDOF)
F[load] = P

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


# ---------- 6.14 (c) ----------

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

# stiffness matrix assembly
for element in elementNodes:
    nodes = nodeCoordinates[element]

    A2 = elementArea(nodes)
    B = GradientMatrix(nodes)    
    
    D = planeStress(E, nu)

    DOF = [2*element[0], 2*element[0] + 1, 2*element[1], 2*element[1] + 1, 2*element[2], 2*element[2] + 1]
    
    K[np.ix_(DOF, DOF)] += t*A2/2*B.T @ D @ B

# boundary conditions
bc = [4, 5, 6, 7]
mask = np.arange(GDOF)
mask = [False if x in bc else True for x in mask]

# load 
P = -50e3
load = [9]
F = np.zeros(GDOF)
F[load] = P

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

# ---------- 6.14 (d) ----------

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

# stiffness matrix assembly
for element in elementNodes:
    nodes = nodeCoordinates[element]

    A2 = elementArea(nodes)
    B = GradientMatrix(nodes)    
    
    D = planeStress(E, nu)

    DOF = [2*element[0], 2*element[0] + 1, 2*element[1], 2*element[1] + 1, 2*element[2], 2*element[2] + 1]
    
    K[np.ix_(DOF, DOF)] += t*A2/2*B.T @ D @ B

# boundary conditions
bc = [0, 1, 6, 7]
mask = np.arange(GDOF)
mask = [False if x in bc else True for x in mask]

# load 
P = 10e6
load = [2, 4]
F = np.zeros(GDOF)
F[load] = P*t*0.4/2

print(F)

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