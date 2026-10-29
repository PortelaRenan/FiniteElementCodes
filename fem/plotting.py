import os
import numpy as np
from datetime import datetime
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection

def beamPlot(displacements, reaction, nodeCoordinates: list[tuple[float, float]]) -> None:
    """
    Plot the displacement of the beam, the shear force and bending moment diagrams.

    Parameters:
        displacements   : 
        reaction        :      
        nodeCoordinates : Nodal coordinates — each entry is an x coordinate value

    Returns:
        None
    """        
    reaction = reaction.reshape((-1, 2))
        
    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(10, 8))

    # displacement
    xDis = displacements.reshape((-1, 2))[:,0]
    ax1.plot(nodeCoordinates, xDis, linestyle = '-')
    ax1.axhline(0, color='black', linewidth=0.8)
    ax1.set_xlabel('Beam coordinate, x (m)')
    ax1.set_ylabel('Displacement (m)')
    ax1.set_title('Nodal Displacement')
    ax1.grid(True)

    # Shear force: first column of reaction vector
    ax2.plot(nodeCoordinates, reaction[:, 0], 'o-', linewidth=2)
    ax2.fill_between(nodeCoordinates, 0, reaction[:, 0], alpha=0.2)
    ax2.axhline(0, color='black', linewidth=0.8)
    ax2.set_xlabel('Beam coordinate, x (m)')
    ax2.set_ylabel('Shear force (N)')
    ax2.set_title('Nodal Shear Force')
    ax2.grid(True)

    # Bending moment: second column of reaction vector
    ax3.plot(nodeCoordinates, reaction[:, 1], 'o-', linewidth=2)
    ax3.fill_between(nodeCoordinates, 0, reaction[:, 1], alpha=0.2)
    ax3.axhline(0, color='black', linewidth=0.8)
    ax3.set_xlabel('Beam coordinate, x (m)')
    ax3.set_ylabel('Bending moment (N·m)')
    ax3.set_title('Nodal Bending Moment')
    ax3.grid(True)

    plt.tight_layout()
    plt.show()