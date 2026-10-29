import sys
from pathlib import Path

sys.path.append(str(Path(__file__).resolve().parent.parent))

import numpy as np
import matplotlib.pyplot as plt
from fem.elements.beam import BernoulliBeam
from fem.plotting import beamPlot

E = 70e9   # modulus od elasticity
I = 4e-4    # second moment of area

# coordinates and connectivities
L = 5
numberElements = 2
nodeCoordinates = np.array((0, 3, 5))
numberNodes = np.size(nodeCoordinates)
elementNodes = [[i, i + 1] for i in range(numberElements)]

# structure
GDOF = 2*numberNodes 
displacements = np.zeros(GDOF)

# load 
P = -80e3
load = [4]
F = np.zeros(GDOF)
F[load] = P/2

# stiffness
stiffness = BernoulliBeam(GDOF, elementNodes, nodeCoordinates, E, I)

# boundary conditions
bc = [0, 2, 5]
mask = np.arange(GDOF)
mask = [False if x in bc else True for x in mask]

# solution
displacements[mask] = np.linalg.solve(stiffness[np.ix_(mask, mask)], 
                                       F[mask])

np.set_printoptions(precision=4)
print(f'Displacement = \n{displacements.reshape((-1, 2))} mm')

# reaction
reaction = np.matmul(stiffness, displacements)
print(f'Reaction = \n{reaction.reshape((-1, 2))/1e3} kN')

# Post-processing: nodal displacement, shear force and bending moment diagrams
beamPlot(displacements, reaction, nodeCoordinates)

