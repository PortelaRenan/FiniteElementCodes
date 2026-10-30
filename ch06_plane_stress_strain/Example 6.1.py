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
t = 20e-3 # thickness

# mesh
nodeCoordinates = np.array([
    [0, -20],
    [40, 0],
    [0, 20]
])*1e-3

# stress calculation
A2 = elementArea(nodeCoordinates)
B = GradientMatrix(nodeCoordinates)

D = planeStress(E, nu)

k = t*A2/2*B.T @ D @ B

u = np.array([0, 0.05, 0.025, 0, 0, 0.05])*1e-3

sigma = D @ B @ u

principalStresses(sigma)