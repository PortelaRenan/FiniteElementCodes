import sys
from pathlib import Path

sys.path.append(str(Path(__file__).resolve().parent.parent))

import numpy as np
import matplotlib.pyplot as plt
from FEA_functions import TimoshenkoBeam

E = 207e9   # modulus of elasticity
G =  80e9   # shear modulus
I = 25e-3*(50e-3)**3    # second moment of area
A = 25e-3*50e-3    # cross section area

# coordinates and connectivities
L = 200e-3
numberElements = 1
nodeCoordinates = np.linspace(0,L,numberElements+1)
numberNodes = np.size(nodeCoordinates)
elementNodes = [[i, i + 1] for i in range(numberElements)]

# structure
GDOF = 2*numberNodes 
displacements = np.zeros(GDOF)

# load 
P = -10e3
load = [2]
F = np.zeros(GDOF)
F[load] = P/2 # symmetry

# stiffness
stiffness = TimoshenkoBeam(GDOF, elementNodes, nodeCoordinates, E, G, I, A)

# boundary conditions
bc = [0, 3] # support @ node 1 + rotation @ node 2
mask = np.arange(GDOF)
mask = [False if x in bc else True for x in mask]

# solution
displacements[mask] = np.linalg.solve(stiffness[np.ix_(mask, mask)], 
                                       F[mask])

np.set_printoptions(precision=4)
print(f'Displacement = \n{displacements.reshape((-1, 2))} mm')

# reaction
reaction = np.matmul(stiffness, displacements)
print(f'Reaction = \n{reaction.reshape((-1, 2))} kN')

# Post-processing: nodal displacement, shear force and bending moment diagrams
nodal_forces = reaction.reshape((-1, 2))

fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(10, 8))

# displacement
xDis = displacements.reshape((-1, 2))[:,0]
ax1.plot(nodeCoordinates, xDis, linestyle = '-')
ax1.axhline(0, color='black', linewidth=0.8)
ax1.set_xlabel('Beam coordinate, x (m)')
ax1.set_ylabel('Displacement (m)')
ax1.set_title('Nodal Shear Force')
ax1.grid(True)

# Shear force: first column of reaction vector
ax2.plot(nodeCoordinates, nodal_forces[:, 0], 'o-', linewidth=2)
ax2.fill_between(nodeCoordinates, 0, nodal_forces[:, 0], alpha=0.2)
ax2.axhline(0, color='black', linewidth=0.8)
ax2.set_xlabel('Beam coordinate, x (m)')
ax2.set_ylabel('Shear force (N)')
ax2.set_title('Nodal Shear Force')
ax2.grid(True)

# Bending moment: second column of reaction vector
ax3.plot(nodeCoordinates, nodal_forces[:, 1], 'o-', linewidth=2)
ax3.fill_between(nodeCoordinates, 0, nodal_forces[:, 1], alpha=0.2)
ax3.axhline(0, color='black', linewidth=0.8)
ax3.set_xlabel('Beam coordinate, x (m)')
ax3.set_ylabel('Bending moment (N·m)')
ax3.set_title('Nodal Bending Moment')
ax3.grid(True)

plt.tight_layout()
plt.show()

displacement = displacements[2]
print(f'Difference = \n{(-(-2.597e-6 - displacement)/displacement)*100}%')