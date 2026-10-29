import os
import numpy as np

def shapeFunctionQ4(xi, eta):
    """
    Shape function matrix of a Q4 element.

    Parameters:
        xi     : gaussian Point float value     
        eta    : gaussian Point float value

    Returns:
        shapeFunction      : (4 x 2) array 
        naturalDerivatives : (4 x 2) array 
    """        
    shapeFunction = 1/4*np.array([[(1 - xi)*(1 - eta)],
                                  [(1 + xi)*(1 - eta)],
                                  [(1 + xi)*(1 + eta)],
                                  [(1 - xi)*(1 + eta)]
    ]) 
    
    # natural derivatives order:
    # [d/dx, d/dy, d^2/dx^2, d^2/dy^2, d^2/dxdy]
    naturalDerivatives = (1/4) * np.array([
        [-(1 - eta), -(1 - xi)],
        [ (1 - eta), -(1 + xi)],
        [ (1 + eta),  (1 + xi)],
        [-(1 + eta),  (1 - xi)]
    ])

    return shapeFunction, naturalDerivatives

def Bmatrix(xi: float, eta: float, nodeXY: list[tuple[float, float]]):
    """
    B matrix of a Q4 element.

    Parameters:
        xi     : gaussian Point float value     
        eta    : gaussian Point float value
        nodeXY : (4, 2) array of element node coordinates [[x1,y1], ... [x4,y4]]

    Returns:
        B matrix: (3 x 8) array B matrix of a Q4 element
        detJ    : determinant of the Jacobian 
    """        

    _, dN_nat = shapeFunctionQ4(xi, eta)      # (4, 2)
    
    # Jacobian: J = [[dx/dxi, dy/dxi], [dx/deta, dy/deta]]
    J = dN_nat.T @ nodeXY                     # (2, 2)
    detJ = np.linalg.det(J)
    invJ = np.linalg.inv(J)
    
    # derivatives with respect to x, y: row 0 = dN/dx, row 1 = dN/dy
    dN_xy = np.linalg.solve(J, dN_nat.T)                   # (2, 4)
    B = np.zeros((3, 8))
    B[0, 0::2] = dN_xy[0, :]  # eps_x = du/dx
    B[1, 1::2] = dN_xy[1, :]  # eps_y = dv/dy
    B[2, 0::2] = dN_xy[1, :]  # gamma_xy = du/dy
    B[2, 1::2] = dN_xy[0, :]  # gamma_xy = dv/dx
    
    return B, detJ 