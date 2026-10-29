import sys
from pathlib import Path

sys.path.append(str(Path(__file__).resolve().parent.parent))

import numpy as np
import matplotlib.pyplot as plt
from fem.elements.beam import BernoulliBeam, distributedLoad
from fem.plotting import beamPlot

E = 200e9   # modulus od elasticity
I = 8e-5    # second moment of area

# coordinates and connectivities
L = 10
numberElements = 2
nodeCoordinates = np.linspace(0,L,numberElements+1)
numberNodes = np.size(nodeCoordinates)
elementNodes = [[i, i + 1] for i in range(numberElements)]

# structure
GDOF = 2*numberNodes
displacements = np.zeros(GDOF)

# load 
w = -30e3
F = distributedLoad(GDOF, elementNodes, nodeCoordinates, w) 

# stiffness
stiffness = BernoulliBeam(GDOF, elementNodes, nodeCoordinates, E, I)

# boundary conditions
bc = [0, 1, 4]
mask = np.arange(GDOF)
mask = [False if x in bc else True for x in mask]

# solution
displacements[mask] = np.linalg.solve(stiffness[np.ix_(mask, mask)], 
                                       F[mask])

np.set_printoptions(precision=4)
print(f'Displacement = \n{displacements} mm')

# reaction
reaction = np.matmul(stiffness, displacements) - F
print(f'Reaction = \n{reaction/1e3} kN')

# Post-processing: nodal displacement, shear force and bending moment diagrams
beamPlot(displacements, reaction, nodeCoordinates) 