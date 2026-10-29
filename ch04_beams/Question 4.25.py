import sys
from pathlib import Path

sys.path.append(str(Path(__file__).resolve().parent.parent))

import numpy as np
import matplotlib.pyplot as plt
from fem.elements.beam import BernoulliBeam, distributedLoad
from fem.plotting import beamPlot

E = 10e9    # modulus od elasticity
I = 4e-5    # second moment of area
k = 800e3   # spring stiffness

# coordinates and connectivities
L = 12
numberElements = 2
nodeCoordinates = np.linspace(0,L,numberElements+1)
numberNodes = np.size(nodeCoordinates)
elementNodes = [[i, i + 1] for i in range(numberElements)]

# structure
GDOF = 2*numberNodes + 1
displacements = np.zeros(GDOF)

 # load 
w = -7.5e3
F = distributedLoad(GDOF, elementNodes, nodeCoordinates, w) 

# stiffness
stiffness = BernoulliBeam(GDOF, elementNodes, nodeCoordinates, E, I)

DOF = [2, 6]    
stiffness[np.ix_(DOF, DOF)] += k*np.array([
    [1, -1],
    [-1, 1]
])

# boundary conditions
bc = [0, 1, 4, 5, 6]
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
beamPlot(displacements[0:GDOF - 1], reaction[0:GDOF - 1], nodeCoordinates[0:GDOF - 1]) 