import sys
from pathlib import Path

sys.path.append(str(Path(__file__).resolve().parent.parent))

import numpy as np

from fem.gauss import gaussPoints

# geometry properties
t = 5e-3 # thickness

# ---------- 6.12 (b) ----------

# mesh
nodeCoordinates = np.array([
    [   0,   0],
    [ 100,   0]
])*1e-3

GDOF = 2*len(nodeCoordinates)

# load 
P = 10e6
L = np.linalg.norm(
    np.array(nodeCoordinates[1]) -
    np.array(nodeCoordinates[0])
)
load = [0, 2]
F = np.zeros(GDOF)
p, w = gaussPoints(n = 10) 
J = L/2

for pi, wi in zip(p, w):
    
    # Shape functions
    N2 = (1 - pi)/2
    N3 = (1 + pi)/2

    # Quadratic pressure distribution
    px = P*np.sin(np.pi/2*(1 + pi))

    # Equivalent nodal forces
    P2 = N2 * px * wi * J * t
    P3 = N3 * px * wi * J * t

    F[load] += [P2, P3]

print(F)