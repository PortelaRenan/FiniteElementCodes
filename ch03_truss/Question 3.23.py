import sys
from pathlib import Path

sys.path.append(str(Path(__file__).resolve().parent.parent))

import numpy as np
from fem.elements.bar import Truss
from fem.plotting import plotTruss

E = 70e9   # modulus od elasticity
A = 6e-4    # cross sectional area


# coordinates and connectivities
nodeCoordinates = np.array([
    [1.25, 1.25*np.sqrt(3)],
    [0, 0],
    [2.5, 0]
])
elementNodes = np.array([
    [0, 1],
    [0, 2]
])

numberNodes = len(nodeCoordinates)
numberElements = len(elementNodes)

# structure
GDOF = 2*numberNodes 
displacements = np.zeros(GDOF)

# load 
P = 60e3
load = [0]
F = np.zeros(GDOF)
F[load] = P

# stiffness
stiffness = Truss(GDOF, elementNodes, nodeCoordinates, E, A)

# boundary conditions
bc = [2, 3, 4, 5]
mask = np.arange(GDOF)
mask = [False if x in bc else True for x in mask]

# solution
displacements[mask] = np.linalg.solve(stiffness[np.ix_(mask, mask)], 
                                       F[mask])

np.set_printoptions(precision=4)
print(f'Displacement = \n{displacements} m')

# reaction
reaction = np.matmul(stiffness, displacements)
print(f'Reaction = \n{reaction/1e3} kN') 

plotTruss(nodeCoordinates, elementNodes, displacements, scale=1e1)