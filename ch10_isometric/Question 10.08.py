import sys
from pathlib import Path

sys.path.append(str(Path(__file__).resolve().parent.parent))

import numpy as np
from fem.gauss import gaussPoints

# material properties & constitutive matrix
E = 205e9 # Young's modulus
A = 2e-4

nodeCoordinates = np.array([[0,0], 
                            [2,0], 
                            [4,0]], dtype=float)


L = np.linalg.norm(
    np.array(nodeCoordinates[2]) -
    np.array(nodeCoordinates[0])
)

# stiffness
K = np.array([
    [ 14/3,  2/3, -16/3],
    [  2/3, 14/3, -16/3],
    [-16/3, -16/3, 32/3]
])*E*A/2/L

# structure
GDOF = len(nodeCoordinates)
displacements = np.zeros(GDOF)

# load 
F = np.zeros(GDOF)
q = 2e3 # kN/m
p, w = gaussPoints(n = 2)

for xi, wi in zip(p, w):

    F += np.array([xi*(xi - 1)/2, xi*(xi + 1)/2, 1 - xi**2])*q*L/2*wi

# boundary conditions
bc = [0]
mask = np.arange(GDOF)
mask = [False if x in bc else True for x in mask]

# solution
displacements[mask] = np.linalg.solve(K[np.ix_(mask, mask)], 
                                       F[mask])

np.set_printoptions(precision=4)
print(f'Displacement = \n{displacements} mm')

# stress at x = 0, L/2, L  ->  xi = -1, 0, 1
xi_pts = np.array([-1.0, 0.0, 1.0])
stress = np.zeros(GDOF)

for i, xi in enumerate(xi_pts):
    B = np.array([(2*xi - 1), (2*xi + 1), -4*xi]) / L
    stress[i] = E * (B @ displacements)

np.set_printoptions(precision=4)
print(f'Stress = {stress / 1e6}')   # MPa