import os
import numpy as np
from datetime import datetime
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection

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

def distributedLoad(GDOF: int, elementNodes: list[tuple[int, int]], nodeCoordinates: list[tuple[float, float]], P: float) -> np.ndarray:
    """
    Assemble the global load vector for a uniformly distributed load on a beam.

    Parameters:
        GDOF            : Total number of global degrees of freedom
        elementNodes    : Connectivity table — each entry is a (node_i, node_j) pair
                          defining the two nodes of an element
        nodeCoordinates : Nodal coordinates — each entry is an x coordinate value
        P               : Magnitude of the uniformly distributed load

    Returns:
        load : np.ndarray, shape (GDOF,)
            Global equivalent nodal load vector assembled from all element contributions
    """

    load = np.zeros(GDOF)
    for element in range(np.size(elementNodes,0)):
        indices = elementNodes[element]
        DOF = [2*indices[0], 2*indices[0] + 1, 2*indices[1], 2*indices[1] + 1]    

        length = nodeCoordinates[indices[1]] - nodeCoordinates[indices[0]]
        
        localLoad = P*length*np.array([1/2, length/12, 1/2, -length/12])

        load[DOF] += localLoad

    return load


def triangularLoad(GDOF: int, elementNodes: list[tuple[int, int]], nodeCoordinates: list[tuple[float, float]], P: float) -> np.ndarray:
    """
    Assemble the global load vector for a linearly increasing distributed load on a beam.

    Parameters:
        GDOF            : Total number of global degrees of freedom
        elementNodes    : Connectivity table — each entry is a (node_i, node_j) pair
                          defining the two nodes of an element
        nodeCoordinates : Nodal coordinates — each entry is an x coordinate value
        P               : Magnitude of the uniformly distributed load

    Returns:
        load : np.ndarray, shape (GDOF,)
            Global equivalent nodal load vector assembled from all element contributions
    """

    load = np.zeros(GDOF)
    for element in range(np.size(elementNodes,0)):
        
        indices = elementNodes[element]
        DOF = [2*indices[0], 2*indices[0] + 1, 2*indices[1], 2*indices[1] + 1]    

        length = nodeCoordinates[indices[1]] - nodeCoordinates[indices[0]]
        
        localLoad = P*length*np.array([1/2, length/12, 1/2, -length/12])

        load[DOF] += localLoad

    return NotImplemented