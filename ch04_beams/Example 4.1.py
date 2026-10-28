import sys
from pathlib import Path

sys.path.append(str(Path(__file__).resolve().parent.parent))

import numpy as np
import matplotlib.pyplot as plt
from fem.elements.beam import BernoulliBeam

E = 1 # modulus od elasticity
I = 1    # second moment of area

# coordinates and connectivities
L = 2
numberElements = 2
nodeCoordinates = np.linspace(0,L,numberElements+1)
numberNodes = np.size(nodeCoordinates)
elementNodes = [[i, i + 1] for i in range(numberElements)]

# structure
GDOF = 2*numberNodes 
displacements = np.zeros(GDOF)

# load 
P = -1
load = [0]
F = np.zeros(GDOF)
F[load] = P

# stiffness
stiffness = BernoulliBeam(GDOF, elementNodes, nodeCoordinates, E, I)

# boundary conditions
bc = [2, 4, 5]
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

# Post-processing
xDis = displacements.reshape((-1, 2))[:,0]
plt.plot(nodeCoordinates, xDis, linestyle = '-')
plt.show()

# Plot shear force and bending moment diagrams
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(9, 7))

for e, element in enumerate(elementNodes):

    i, j = element
    x1 = nodeCoordinates[i]
    x2 = nodeCoordinates[j]
    Le = x2 - x1

    # Recalculate element end forces
    dofs = [2*i, 2*i + 1, 2*j, 2*j + 1]
    de = displacements[dofs]

    ke = (E * I / Le**3) * np.array([
        [12,       6*Le,  -12,       6*Le],
        [6*Le,  4*Le**2,  -6*Le,  2*Le**2],
        [-12,     -6*Le,   12,      -6*Le],
        [6*Le,  2*Le**2,  -6*Le,  4*Le**2]
    ])

    fe = ke @ de

    # Local coordinate along the element
    x_local = np.linspace(0, Le, 50)
    x_global = x1 + x_local

    # Internal shear and bending moment (sign convention)
    V = np.full_like(x_local, fe[0])
    M = -fe[1] + fe[0] * x_local

    # Shear force diagram
    ax1.plot(x_global, V, linewidth=2)
    ax1.fill_between(x_global, 0, V, alpha=0.2)

    # Bending moment diagram
    ax2.plot(x_global, M, linewidth=2)
    ax2.fill_between(x_global, 0, M, alpha=0.2)

# Formatting
ax1.axhline(0, color='black', linewidth=0.8)
ax1.set_ylabel('Shear force')
ax1.set_title('Shear Force Diagram')
ax1.grid(True)

ax2.axhline(0, color='black', linewidth=0.8)
ax2.set_xlabel('Beam coordinate, x')
ax2.set_ylabel('Bending moment')
ax2.set_title('Bending Moment Diagram')
ax2.grid(True)

plt.tight_layout()
plt.show()