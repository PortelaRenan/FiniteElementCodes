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


def calculateBendingMoment(displacements, elementNodes,
                           nodeCoordinates, E, I):

    xMoment = []
    bendingMoment = []

    for nodes in elementNodes:

        node_i, node_j = nodes

        Le = nodeCoordinates[node_j] - nodeCoordinates[node_i]

        # Element DOFs: [v_i, theta_i, v_j, theta_j]
        dofs = [
            2*node_i,
            2*node_i + 1,
            2*node_j,
            2*node_j + 1
        ]

        de = displacements[dofs]

        # Evaluate moment at 10 points per element
        xi_values = np.linspace(0, 1, 10)

        for xi in xi_values:

            # Second derivatives of Hermite shape functions
            d2N = np.array([
                (-6 + 12*xi) / Le**2,
                (-4 + 6*xi) / Le,
                ( 6 - 12*xi) / Le**2,
                (-2 + 6*xi) / Le
            ])

            curvature = d2N @ de

            M = E * I * curvature

            x = nodeCoordinates[node_i] + xi * Le

            xMoment.append(x)
            bendingMoment.append(M)

    return np.array(xMoment), np.array(bendingMoment)

def rightHingedBeam(E: float, I: float, elementNodes: tuple[int, int], nodeCoordinates: list[float]) -> tuple[list[int], np.ndarray]:
    """
    Calculate the local stiffness matrix of a two-node Euler-Bernoulli beam element
    with a hinge at its right end (node j), and the global DOF indices it maps to.

    Parameters:
        E               : Young's modulus of the material
        I               : Second moment of area of the beam cross-section
        elementNodes    : Node indices (node_i, node_j) defining the two ends of this
                          single element
        nodeCoordinates : x coordinate of every node in the mesh, indexed by node number

    Returns:
        DOF            : list[int], length 3
            Global DOF indices of the element, in the same order as the rows and
            columns of localStiffness: [v_i, theta_i, v_j]
        localStiffness : np.ndarray, shape (3, 3)
            Element stiffness matrix with the rotation at node j condensed out
            (hinge), for the DOFs [v_i, theta_i, v_j]
    """

    # Global DOF indices (node i: displacement and rotation; node j: displacement only, rotation released by the hinge)
    DOF = [2*elementNodes[0], 2*elementNodes[0] + 1, 2*elementNodes[1]]

    # Element length from the x coordinates of its end nodes
    length = nodeCoordinates[elementNodes[1]] - nodeCoordinates[elementNodes[0]]

    # Stiffness of a propped-cantilever style element: (3EI / L^3) * [...]
    localStiffness = 3*E*I/length**3*np.array([
        [1, length, -1],
        [length, length**2, -length],
        [-1, -length, 1]
    ])

    return DOF, localStiffness


def leftHingedBeam(E: float, I: float, elementNodes: tuple[int, int], nodeCoordinates: list[float]) -> tuple[list[int], np.ndarray]:
    """
    Calculate the local stiffness matrix of a two-node Euler-Bernoulli beam element
    with a hinge at its left end (node i), and the global DOF indices it maps to.

    Parameters:
        E               : Young's modulus of the material
        I               : Second moment of area of the beam cross-section
        elementNodes    : Node indices (node_i, node_j) defining the two ends of this
                          single element
        nodeCoordinates : x coordinate of every node in the mesh, indexed by node number

    Returns:
        DOF            : list[int], length 3
            Global DOF indices of the element, in the same order as the rows and
            columns of localStiffness: [v_i, v_j, theta_j]
        localStiffness : np.ndarray, shape (3, 3)
            Element stiffness matrix with the rotation at node i condensed out
            (hinge), for the DOFs [v_i, v_j, theta_j]
    """

    # Global DOF indices (node i: displacement only, rotation released by the hinge; node j: displacement and rotation)
    DOF = [2*elementNodes[0], 2*elementNodes[1], 2*elementNodes[1] + 1]

    # Element length from the x coordinates of its end nodes
    length = nodeCoordinates[elementNodes[1]] - nodeCoordinates[elementNodes[0]]

    # Stiffness of a propped-cantilever style element: (3EI / L^3) * [...]
    localStiffness = 3*E*I/length**3*np.array([
        [ 1, -1,  length],
        [-1,  1, -length],
        [length, -length, length**2],
    ])

    return DOF, localStiffness

def localBernoulliBeam(E: float, I: float, elementNodes: tuple[int, int], nodeCoordinates: list[float]) -> tuple[list[int], np.ndarray]:
    """
    Calculate the local stiffness matrix of a two-node Euler-Bernoulli beam element.

    Parameters:
            E               : Young's modulus of the material
            I               : Second moment of area of the beam cross-section
            elementNodes    : Node indices (node_i, node_j) defining the two ends of this
                              single element
            nodeCoordinates : x coordinate of every node in the mesh, indexed by node number

    Returns:
        DOF            : list[int], length 4
                    Global DOF indices of the element, in the same order as the rows and
                    columns of localStiffness (used to assemble into the global matrix)

        localStiffness : np.ndarray, shape (4, 4)
            Element stiffness matrix for the DOFs [v_i, theta_i, v_j, theta_j]
            (transverse displacement and rotation at each node)
    """

    # Global DOF indices of the element (2 DOFs per node: displacement and rotation)
    DOF = [2*elementNodes[0], 2*elementNodes[0] + 1, 2*elementNodes[1], 2*elementNodes[1] + 1]

    # Element length from the x coordinates of its end nodes
    length = nodeCoordinates[elementNodes[1]] - nodeCoordinates[elementNodes[0]]

    # Standard Euler-Bernoulli beam stiffness matrix: (EI / L^3) * [...]
    localStiffness = E*I/length**3*np.array([
                [12, 6*length, -12, 6*length],
                [6*length, 4*length**2, -6*length, 2*length**2],
                [-12, -6*length, 12, -6*length],
                [6*length, 2*length**2, -6*length, 4*length**2]
            ])

    return DOF, localStiffness

def localDistributedLoad(P: float, elementNodes: list[tuple[int, int]], nodeCoordinates: list[tuple[float, float]] ) -> tuple[list[int], np.ndarray]:
    """
    Calculate local load vector for a uniformly distributed load on a beam.

    Parameters:
        P               : Magnitude of the uniformly distributed load
        elementNodes    : Node indices (node_i, node_j) defining the two ends of this
                                      single element
        nodeCoordinates : Nodal coordinates — each entry is an x coordinate value
        

    Returns:
        DOF            : list[int], length 4
                            Global DOF indices of the element, in the same order as the rows and
                            columns of localStiffness (used to assemble into the global matrix)

        load : np.ndarray, shape (4,)
            Global equivalent nodal load vector assembled from all element contributions
    """

    # Global DOF indices of the element (2 DOFs per node: displacement and rotation)
    DOF = [2*elementNodes[0], 2*elementNodes[0] + 1, 2*elementNodes[1], 2*elementNodes[1] + 1]

    # Element length from the x coordinates of its end nodes
    length = nodeCoordinates[elementNodes[1]] - nodeCoordinates[elementNodes[0]]

    # Distributed load vector for a single element
    load = P*length*np.array([1/2, length/12, 1/2, -length/12])
    
    return DOF, load