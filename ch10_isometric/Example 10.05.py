import sys
from pathlib import Path

sys.path.append(str(Path(__file__).resolve().parent.parent))

import numpy as np
from fem.elements.Q4 import shapeFunctionQ4, Bmatrix
from fem.constitutiveMatrix import planeStress
from fem.gauss import gaussPoints


E = 200e9 # Young's modulus
nu = 0.3 # Poisson's ration
h = 1e-2  # element thickness

nodeCoordinates = np.array([[3,2], [5,2], [5,4], [3,4]], dtype=float)

D = planeStress(E, nu)

p, w = gaussPoints(n = 2)

K = np.zeros((8, 8))

gaussPoints2D = [
    (p[0], p[0]),  # (-1/sqrt(3), -1/sqrt(3))
    (p[1], p[0]),  # (+1/sqrt(3), -1/sqrt(3))
    (p[1], p[1]),  # (+1/sqrt(3), +1/sqrt(3))
    (p[0], p[1])   # (-1/sqrt(3), +1/sqrt(3))
]

B, detJ = Bmatrix(0, 0, nodeCoordinates)

d = np.array([0, 0, 0.02, 0.03, 0.06, 0.032, 0, 0])*1e-3

sigma = D @ B @ d
print(sigma/1e6)