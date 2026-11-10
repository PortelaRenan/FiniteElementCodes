import sys
from pathlib import Path

sys.path.append(str(Path(__file__).resolve().parent.parent))

import numpy as np
import matplotlib.pyplot as plt
from fem.elements.beam import localBernoulliBeam, rightHingedBeam, localDistributedLoad
from fem.displacementResults import displacementResults


E = 210e9    # modulus od elasticity
I = 2e-4     # moment of Inertia

# coordinates and connectivities
L = 4               # Length of beam in meters
numberElements = 3
nodeCoordinates = np.array([0, 2, 3, 4])
numberNodes = np.size(nodeCoordinates)
elementNodes = [[i, i + 1] for i in range(numberElements)]

# structure
GDOF = 2*numberNodes 
displacements = np.zeros(GDOF)

# load 
q = -10e3
F = np.zeros(GDOF)

# boundary conditions
bc = [0, 1, 2, 6, 7]
mask = np.arange(GDOF)
mask = [False if x in bc else True for x in mask]   

# stiffness assembly
stiffness = np.zeros((GDOF, GDOF))

for element in range(np.size(elementNodes,0)):

    # local stiffness and DOFs
    if element == 1:
        DOF, localStiffness = rightHingedBeam(E, I, elementNodes[element], nodeCoordinates)
    else:
        DOF, localStiffness = localBernoulliBeam(E, I, elementNodes[element], nodeCoordinates)

    stiffness[np.ix_(DOF, DOF)] += localStiffness

    # local force
    if element == 2:
        DOF, localLoad = localDistributedLoad(q, elementNodes[element], nodeCoordinates)

        F[DOF] = localLoad

# solution
displacements[mask] = np.linalg.solve(stiffness[np.ix_(mask, mask)], 
                                       F[mask])
displacementResults(displacements)