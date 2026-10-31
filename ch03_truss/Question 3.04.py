import sys
from pathlib import Path

sys.path.append(str(Path(__file__).resolve().parent.parent))

import numpy as np
from fem.elements.bar import Bar

E = 200e9   # modulus od elasticity
A = 30e-4   # cross sectional area


# coordinates and connectivities
L = 0.75*3
numberElements = 3
nodeCoordinates = np.linspace(0,L,numberElements+1)
numberNodes = np.size(nodeCoordinates)
elementNodes = [[i, i + 1] for i in range(numberElements)]

# structure
GDOF = numberNodes 
displacements = np.zeros(GDOF)

# load 
P = [20e3, -50e3]
load = [1, 2]
F = np.zeros(GDOF)
F[load] = P

# stiffness
stiffness = Bar(GDOF, elementNodes, nodeCoordinates, E, A)

# boundary conditions
bc = [0, 3]
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
