import sys
from pathlib import Path

sys.path.append(str(Path(__file__).resolve().parent.parent))

import numpy as np
from fem.gauss import gaussPoints

E = 210e9   # Young's modulus
A = 12.5e-4 # Area

# coordinates and connectivities
L = 1.5
numberElements = 2
nodeCoordinates = np.linspace(0,L,5)
numberNodes = np.size(nodeCoordinates)
elementNodes = np.array([
    [0, 1, 2],
    [2, 3, 4]
])

# structure
GDOF = len(nodeCoordinates)
displacements = np.zeros(GDOF)
K = np.zeros((GDOF, GDOF))

# stiffness matrix assembly
for element in elementNodes:
    nodes = nodeCoordinates[element]

    DOF = [element[0], element[2], element[1]]

    K[np.ix_(DOF, DOF)] += np.array([   [ 14/3,  2/3, -16/3],
                                        [  2/3, 14/3, -16/3],
                                        [-16/3, -16/3, 32/3]
                                    ])*E*A/2/(L/2)

# load 
F = np.zeros(GDOF)
q = 2e3 # kN/m
p, w = gaussPoints(n = 2)

# force vector assembly
for element in elementNodes:
    xa, xb = nodeCoordinates[element[0]], nodeCoordinates[element[2]]
    Le = xb - xa

    force = np.zeros(3)
    for xi, wi in zip(p, w):
        x = xa + Le/2*(1 + xi)                # physical position of Gauss point
        distributedLoad = q * x / L           # 0 at x = 0, q0 at x = L_total
        N = np.array([xi*(xi - 1)/2, xi*(xi + 1)/2, 1 - xi**2])
        force += N * distributedLoad * (Le/2) * wi
   
    
    nodes = nodeCoordinates[element]

    DOF = [element[0], element[2], element[1]]

    F[DOF] += force

# boundary conditions
bc = [4]
mask = np.arange(GDOF)
mask = [False if x in bc else True for x in mask]

# solution
displacements[mask] = np.linalg.solve(K[np.ix_(mask, mask)], 
                                       F[mask])

np.set_printoptions(precision=4)
print(f'Displacement = \n{displacements} mm')

# reaction
reaction = K @ displacements - F
print(f'Reaction = \n{reaction/1e3} kN')


# stress at x = 0, L/2, L  ->  xi = -1, 0, 1
xi_pts = np.array([-1.0, 0.0, 1.0])          # first end, midside, last end
stress = np.zeros(GDOF)
count  = np.zeros(GDOF)

for element in elementNodes:
    xa, xb = nodeCoordinates[element[0]], nodeCoordinates[element[2]]
    Le = xb - xa

    DOF = [element[0], element[2], element[1]]
    ue  = displacements[DOF]

    B = np.column_stack([2*xi_pts - 1, 2*xi_pts + 1, -4*xi_pts]) / Le   # (3, 3)
    stressElement = E * (B @ ue)                                         # (3,)

    nodesE = [element[0], element[1], element[2]]
    stress[nodesE] += stressElement
    count[nodesE]  += 1

stress /= count
print(f'Stress (MPa) = {stress / 1e6}')