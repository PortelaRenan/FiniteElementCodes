import numpy as np

def elementArea(nodeCoordinates: list[tuple[float, float]])-> np.float64:
    """ 
    Function gives twice the area of a triangular element
     
    """
   
    x1 = nodeCoordinates[0, 0]
    x2 = nodeCoordinates[1, 0]
    x3 = nodeCoordinates[2, 0]

    y1 = nodeCoordinates[0, 1]
    y2 = nodeCoordinates[1, 1]
    y3 = nodeCoordinates[2, 1]
    return np.linalg.det([
        [1, x1, y1],
        [1, x2, y2],
        [1, x3, y3]
    ])

def GradientMatrix(nodeCoordinates: list[tuple[float, float]])-> np.ndarray:

    A2 = elementArea(nodeCoordinates)
    
    x1 = nodeCoordinates[0, 0]
    x2 = nodeCoordinates[1, 0]
    x3 = nodeCoordinates[2, 0]

    y1 = nodeCoordinates[0, 1]
    y2 = nodeCoordinates[1, 1]
    y3 = nodeCoordinates[2, 1]

    B = 1/A2*np.array([
        [y2 - y3, 0, y3 - y1, 0, y1 - y2, 0],
        [0, x3 - x2, 0, x1 - x3, 0, x2 - x1],
        [x3 - x2, y2 - y3, x1 - x3, y3 - y1, x2 - x1, y1 - y2]
    ])

    return B

def principalStresses(sigma: np.array) -> None:

    sigma1 = (sigma[0] + sigma[1])/2 + np.sqrt(((sigma[0] - sigma[1])/2)**2 + sigma[2]**2)
    sigma2 = (sigma[0] + sigma[1])/2 - np.sqrt(((sigma[0] - sigma[1])/2)**2 + sigma[2]**2)

    angle = 0.5*np.atan(2*sigma[2]/(sigma[0] - sigma[1]))*180/np.pi

    print(f'Sigma 1 = {sigma1/1e6:.2f} MPa, \nSigma 2 = {sigma2/1e6:.2f} MPa, \ntheta = {angle:.2f} degrees')

    vonMises = 1/np.sqrt(2)*( (sigma1 - sigma2)**2 + sigma1**2 + sigma2**2 )**(1/2)

    print(f'von Mises stress = {vonMises/1e6:.2f} MPa')