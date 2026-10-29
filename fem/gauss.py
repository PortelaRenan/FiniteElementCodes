import os
import numpy as np

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