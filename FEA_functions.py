import os
import numpy as np
from datetime import datetime
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection

def gaussPoints(n: int, a: float =-1, b: float = 1) -> tuple[list[float], list[float]]:
    """
    Compute Gauss-Legendre quadrature points and weights.
    
    Parameters:
        n : int   - number of Gauss points
        a : float - lower bound of integration (default -1)
        b : float - upper bound of integration (default  1)
    
    Returns:
        points  : array of quadrature points
        weights : array of quadrature weights
    """
    # Get points and weights on [-1, 1]
    points, weights = np.polynomial.legendre.leggauss(n)
    
    # Scale to [a, b]
    points  = 0.5 * (b - a) * points + 0.5 * (b + a)
    weights = 0.5 * (b - a) * weights
    
    return points, weights

# Beam Problems

def BernoulliBeam(GDOF: int, elementNodes: list[tuple[int, int]], nodeCoordinates: list[tuple[float, float]], E: float, I: float) -> np.ndarray:
    """
    Assemble the global stiffness matrix for an Euler-Bernoulli beam structure.

    Parameters:
        GDOF            : Total number of global degrees of freedom
        elementNodes    : Connectivity table — each entry is a (node_i, node_j) pair
                          defining the two nodes of an element
        nodeCoordinates : Nodal coordinates — each entry is an x coordinate value
        E               : Young's modulus of the material
        I               : Second moment of area of the beam cross-section

    Returns:
        stiffness : np.ndarray, shape (GDOF, GDOF)
            Global stiffness matrix assembled from all element contributions
    """


    stiffness = np.zeros((GDOF, GDOF))

    for element in range(np.size(elementNodes,0)):
        indices = elementNodes[element]
        DOF = [2*indices[0], 2*indices[0] + 1, 2*indices[1], 2*indices[1] + 1]    

        length = nodeCoordinates[indices[1]] - nodeCoordinates[indices[0]]
        localStiffness = E*I/length**3*np.array([
            [12, 6*length, -12, 6*length],
            [6*length, 4*length**2, -6*length, 2*length**2],
            [-12, -6*length, 12, -6*length],
            [6*length, 2*length**2, -6*length, 4*length**2]
        ])

        stiffness[np.ix_(DOF, DOF)] += localStiffness

    return stiffness

def TimoshenkoBeam(GDOF: int, elementNodes: list[tuple[int, int]], nodeCoordinates: list[tuple[float, float]], E: float, G: float, I: float, A: float) -> np.ndarray:
    """
    Assemble the global stiffness matrix for an Timoshenko beam structure.

    Parameters:
        GDOF            : Total number of global degrees of freedom
        elementNodes    : Connectivity table — each entry is a (node_i, node_j) pair
                          defining the two nodes of an element
        nodeCoordinates : Nodal coordinates — each entry is an x coordinate value
        E               : Young's modulus of the material
        G               : Shear modulus of the material
        I               : Second moment of area of the beam cross-section
        A               : Cross section area

    Returns:
        stiffness : np.ndarray, shape (GDOF, GDOF)
            Global stiffness matrix assembled from all element contributions
    """


    stiffness = np.zeros((GDOF, GDOF))

    for element in range(np.size(elementNodes,0)):
        indices = elementNodes[element]
        DOF = [2*indices[0], 2*indices[0] + 1, 2*indices[1], 2*indices[1] + 1]    

        ks = 5/6
        length = nodeCoordinates[indices[1]] - nodeCoordinates[indices[0]]
        phi = 12*E*I/(ks*A*G*length**2)
        localStiffness = E*I/length**3*(1 + phi)*np.array([
            [12, 6*length, -12, 6*length],
            [6*length, (4 + phi)*length**2, -6*length, (2 - phi)*length**2],
            [-12, -6*length, 12, -6*length],
            [6*length, (2 - phi)*length**2, -6*length, (4 + phi)*length**2]
        ])

        stiffness[np.ix_(DOF, DOF)] += localStiffness

    return stiffness