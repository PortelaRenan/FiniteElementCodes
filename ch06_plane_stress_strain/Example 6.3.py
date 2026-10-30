import sys
from pathlib import Path

sys.path.append(str(Path(__file__).resolve().parent.parent))

import numpy as np
from fem.elements.Q4 import GradientMatrixQ4
from fem.constitutiveMatrix import planeStress
from fem.gauss import gaussPoints

b = 100e-3
h =  50e-3
t =  25e-3

# material properties
E = 210e9
nu = 0.3
D = planeStress(E, nu)

p, w = gaussPoints(n = 2)

gaussPoints2D = [
    (p[0], p[0]),  # (-1/sqrt(3), -1/sqrt(3))
    (p[1], p[0]),  # (+1/sqrt(3), -1/sqrt(3))
    (p[1], p[1]),  # (+1/sqrt(3), +1/sqrt(3))
    (p[0], p[1])   # (-1/sqrt(3), +1/sqrt(3))
]

K = np.zeros((8, 8))

for xi, eta in gaussPoints2D:
    B = GradientMatrixQ4(b*xi, h*eta, b, h)
    
    K += B.T @ D @ B * t * h * b

print(K)