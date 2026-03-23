# Libraries
import numpy as np
from scipy.linalg import eig
from time import perf_counter as pc

numberElements    = 2
numberNodes   = numberElements + 1
L      = 1
GDOF = 2*numberNodes # all degrees of freedom
kf = np.zeros((2*numberNodes,2*numberNodes))    # global stiffness matrix pre-allocation
kg = np.zeros((2*numberNodes,2*numberNodes))    # global stiffness matrix pre-allocation
coord = [[i + 1, i*L/(numberNodes-1)] for i in range(numberNodes)]
inci = [[i + 1, i + 1, i + 2, coord[int(i)][1], coord[int(i) + 1][1]] for i in range(numberElements)]

# Material properties
E = 200e9

# Geometry properties
# first element
D1 = 20e-3
d1 = 18e-3
# second element
D2 = 20e-3
d2 = 16e-3

l = L/numberElements

geo = np.array([[D1, d1],[D2, d2]])

# Boundary conditions
fixedDOF = np.array([0, 1])
mask = np.arange(GDOF)
mask = np.bool([0 if x in fixedDOF else 1 for x in mask])

for i in range(numberElements):
    node1 = inci[i][1] # first node element
    node2 = inci[i][2] # second node element
    D = geo[i][0]
    d = geo[i][1]
    inertia = np.pi/64*(D**4 - d**4)
    
    # local stiffness matrix
    kf_e = E*inertia/l**3*np.array([[  12, 6*l, -12, 6*l],
                                    [ 6*l,  4*l**2, -6*l,  2*l**2],
                                    [-12,-6*l,  12,-6*l],
                                    [ 6*l,  2*l**2, -6*l,  4*l**2]])
        
    # local geometric matrix
    kg_e = 1/30/l*np.array([[ 36, 3*l, -36, 3*l],
                            [ 3*l,4*l**2, -3*l, -l**2],
                            [-36,-3*l,  36,-3*l],
                            [ 3*l,-l**2, -3*l, 4*l**2]])
        
    # localization vector
    DOF = [2*node1-2,2*node1-1,2*node2-2,2*node2-1]

    kf[np.ix_(DOF, DOF)] += kf_e
    kg[np.ix_(DOF, DOF)] += kg_e

w, v = eig(kf[np.ix_(mask, mask)], 
                     kg[np.ix_(mask, mask)])

Pcr = np.min(np.real(w))

print(f'Pcr = {Pcr: .2f} N')