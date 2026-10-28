import numpy as np
import matplotlib.pyplot as plt

E1 = 200e9   # modulus od elasticity
A1 =  6e-4    # cross sectional area
E2 = 100e9   # modulus od elasticity
A2 = 12e-4    # cross sectional area


# coordinates and connectivities
L = 1.8
numberElements = 3
nodeCoordinates = np.linspace(0,L,numberElements+1)
numberNodes = np.size(nodeCoordinates)
elementNodes = [[i, i + 1] for i in range(numberElements)]

# structure
GDOF = numberNodes 
displacements = np.zeros(GDOF)

# load 
P = 15000
load = [1]
F = np.zeros(GDOF)
F[load] = P

# stiffness
Le = L/numberElements
A = [A1, A1, A2]
E = [E1, E1, E2]

stiffness = np.zeros((GDOF, GDOF))
for i in range(numberElements):
    localStiffness = E[i]*A[i]/Le*np.array([[ 1, -1],
                                             [-1, 1]])

    indices = elementNodes[i]
    DOF = [indices[0], indices[1]]  

    stiffness[np.ix_(DOF, DOF)] += localStiffness
    print(stiffness)

# boundary conditions
bc = [0, 3]
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
