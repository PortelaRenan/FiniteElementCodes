import numpy as np

def Bar(GDOF: int, elementNodes: list[tuple[int, int]], nodeCoordinates: list[tuple[float, float]], E: float, A: float) -> np.ndarray:
    """
    Assemble the global stiffness matrix for a bar structure.

    Parameters:
        GDOF            : Total number of global degrees of freedom
        elementNodes    : Connectivity table — each entry is a (node_i, node_j) pair
                          defining the two nodes of an element
        nodeCoordinates : Nodal coordinates — each entry is an x coordinate value
        E               : Young's modulus of the material
        A               : Cross-section area

    Returns:
        stiffness : np.ndarray, shape (GDOF, GDOF)
            Global stiffness matrix assembled from all element contributions
    """


    stiffness = np.zeros((GDOF, GDOF))

    for element in range(np.size(elementNodes,0)):
        indices = elementNodes[element]
        DOF = [indices[0], indices[1]]    

        length = nodeCoordinates[indices[1]] - nodeCoordinates[indices[0]]
        localStiffness = E*A/length*np.array([
            [ 1, -1],
            [-1,  1]
        ])

        stiffness[np.ix_(DOF, DOF)] += localStiffness

    return stiffness

def Truss(GDOF: int, elementNodes: list[tuple[int, int]], nodeCoordinates: list[tuple[float, float]], E: float, A: float) -> np.ndarray:
    """
    Assemble the global stiffness matrix for a truss structure.

    Parameters:
        GDOF            : Total number of global degrees of freedom
        elementNodes    : Connectivity table — each entry is a (node_i, node_j) pair
                          defining the two nodes of an element
        nodeCoordinates : Nodal coordinates — each entry is an x coordinate value
        E               : Young's modulus of the material
        A               : Cross-section area

    Returns:
        stiffness : np.ndarray, shape (GDOF, GDOF)
            Global stiffness matrix assembled from all element contributions
    """


    stiffness = np.zeros((GDOF, GDOF))

    for element in range(np.size(elementNodes,0)):
        indices = elementNodes[element]

        x1 = nodeCoordinates[indices[0], 0]
        x2 = nodeCoordinates[indices[1], 0]

        y1 = nodeCoordinates[indices[0], 1]
        y2 = nodeCoordinates[indices[1], 1]

        length = np.sqrt((x2 - x1)**2 + (y2 - y1)**2)

        if (x2 - x1) == 0:
            if y2 > y1:
                beta = 2*np.atan(1)
            else:
                beta = -2*np.atan(1)
        else:
            beta = np.atan((y2 - y1)/ (x2 - x1))

        c = np.cos(beta)
        s = np.sin(beta)

        localStiffness = E*A/length*np.array([
            [c*c, c*s, -c*c, -c*s],
            [c*s, s*s, -c*s, -s*s],
            [-c*c, -c*s, c*c, c*s],
            [-c*s, -s*s, c*s, s*s]
        ])

        DOF = [2*indices[0], 2*indices[0] + 1, 2*indices[1], 2*indices[1] + 1]    

        stiffness[np.ix_(DOF, DOF)] += localStiffness

    return stiffness