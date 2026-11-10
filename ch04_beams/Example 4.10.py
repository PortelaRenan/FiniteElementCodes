import sys
from pathlib import Path

sys.path.append(str(Path(__file__).resolve().parent.parent))

import numpy as np
import matplotlib.pyplot as plt
from fem.elements.beam import localBernoulliBeam, rightHingedBeam, calculateBendingMoment
from fem.plotting import beamResults


E = 1    # modulus od elasticity
I = 1    # moment of Inertia

# coordinates and connectivities
L = 3               # Length of beam in meters
numberElements = 2
nodeCoordinates = np.array([0, 2, 3])
numberNodes = np.size(nodeCoordinates)
elementNodes = [[i, i + 1] for i in range(numberElements)]

# structure
GDOF = 2*numberNodes 
displacements = np.zeros(GDOF)

# load 
P = -1         
load = [2]
F = np.zeros(GDOF)
F[load] = P

# boundary conditions
bc = [0, 1, 4, 5]
mask = np.arange(GDOF)
mask = [False if x in bc else True for x in mask]   

# stiffness assembly
stiffness = np.zeros((GDOF, GDOF))

for element in range(np.size(elementNodes,0)):

    if element == 0:
        DOF, localStiffness = rightHingedBeam(E, I, elementNodes[element], nodeCoordinates)
    else:
        DOF, localStiffness = localBernoulliBeam(E, I, elementNodes[element], nodeCoordinates)

    stiffness[np.ix_(DOF, DOF)] += localStiffness

# solution
displacements[mask] = np.linalg.solve(stiffness[np.ix_(mask, mask)], 
                                       F[mask])
print(f'Max Displacement = \n{displacements[2]*1e3:.2f} mm')

# reaction
reaction = stiffness @ displacements 
reaction = reaction.reshape((-1, 2))
maxReaction = max(abs(reaction[:, 0]))
print(f'Max Reaction = \n{maxReaction/1e3:.2f} kN')

xMoment, M = calculateBendingMoment(displacements, elementNodes, nodeCoordinates, E, I)

maxBending = max(abs(M))
print(f'Max Bending Moment = \n{maxBending/1e3:.2f} kN.m')

beamResults(displacements, reaction, nodeCoordinates, xMoment, M)