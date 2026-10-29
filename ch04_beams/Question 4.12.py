import sys
from pathlib import Path

sys.path.append(str(Path(__file__).resolve().parent.parent))

import numpy as np
import matplotlib.pyplot as plt
from fem.elements.beam import BernoulliBeam
from fem.plotting import beamPlot

E = 200e9   # modulus od elasticity
I = 8e-5    # second moment of area
k = 200e3   # spring stiffness

# coordinates and connectivities
L = 6
numberElements = 1
nodeCoordinates = np.linspace(0,L,numberElements+1)
numberNodes = np.size(nodeCoordinates)
elementNodes = [[i, i + 1] for i in range(numberElements)]

# structure
GDOF = 2*numberNodes + 1 
displacements = np.zeros(GDOF)

# load 
P = -10e3
load = [2]
F = np.zeros(GDOF)
F[load] = P

# stiffness
stiffness = BernoulliBeam(GDOF, elementNodes, nodeCoordinates, E, I)

DOF = [2, 4]    
stiffness[np.ix_(DOF, DOF)] += k*np.array([
    [1, -1],
    [-1, 1]
])

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
reaction = np.matmul(stiffness, displacements)
print(f'Reaction = \n{reaction/1e3} kN')

# Post-processing: nodal displacement, shear force and bending moment diagrams
beamPlot(displacements[0:GDOF - 1], reaction[0:GDOF - 1], nodeCoordinates[0:GDOF - 1]) 