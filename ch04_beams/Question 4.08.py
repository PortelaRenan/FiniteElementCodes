import sys
from pathlib import Path

sys.path.append(str(Path(__file__).resolve().parent.parent))

import numpy as np
import matplotlib.pyplot as plt
from fem.elements.beam import BernoulliBeam
from fem.plotting import beamPlot

E = 210e9   # modulus od elasticity
I = 4e-5    # second moment of area

# coordinates and connectivities
L = 12
numberElements = 2
nodeCoordinates = np.linspace(0,L,numberElements+1)
numberNodes = np.size(nodeCoordinates)
elementNodes = [[i, i + 1] for i in range(numberElements)]

# structure
GDOF = 2*numberNodes 
displacements = np.zeros(GDOF)

# load 
P = -5e3
load = [4]
F = np.zeros(GDOF)
F[load] = P

delta = -10e-3
setDisplacement = [2]
vectorDisplacement = np.zeros(GDOF)
vectorDisplacement[setDisplacement] = delta

# stiffness
stiffness = BernoulliBeam(GDOF, elementNodes, nodeCoordinates, E, I)
F0 = stiffness @ vectorDisplacement

# boundary conditions
bc = [0, 1, 2]
mask = np.arange(GDOF)
mask = [False if x in bc else True for x in mask]

# solution
FR = F - F0
displacements[mask] = np.linalg.solve(stiffness[np.ix_(mask, mask)], 
                                       FR[mask])

displacements += vectorDisplacement

np.set_printoptions(precision=4)
print(f'Displacement = \n{displacements.reshape((-1, 2))} mm')

# reaction
reaction = np.matmul(stiffness, displacements)
print(f'Reaction = \n{reaction.reshape((-1, 2))/1e3} kN')

# Post-processing: nodal displacement, shear force and bending moment diagrams
beamPlot(displacements, reaction, nodeCoordinates)