import os
import numpy as np
from datetime import datetime
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
from matplotlib.lines import Line2D

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

def plotTruss(nodeCoordinates, elementNodes, displacements, scale=None):
    u = displacements.reshape((-1, 2))                # (numberNodes, 2): [ux, uy]
    deformed = nodeCoordinates + (scale or 1)*u

    fig, ax = plt.subplots(figsize=(6, 6))
    for n1, n2 in elementNodes:
        ax.plot(*zip(nodeCoordinates[n1], nodeCoordinates[n2]),
                'b--', lw=1.5, label='original' if n1 == elementNodes[0][0] and n2 == elementNodes[0][1] else None)
        ax.plot(*zip(deformed[n1], deformed[n2]),
                'r-', lw=2, label=f'deformed (x{scale:g})' if (n1, n2) == tuple(elementNodes[0]) else None)

    ax.plot(*nodeCoordinates.T, 'bo')
    ax.plot(*deformed.T, 'ro')
    for i, (x, y) in enumerate(nodeCoordinates):
        ax.annotate(str(i + 1), (x, y), textcoords='offset points', xytext=(6, 6))

    ax.set_aspect('equal')
    ax.grid(True)
    ax.set_xlabel('x (m)'); ax.set_ylabel('y (m)')
    ax.set_title('Truss: original and deformed')
    ax.legend()
    plt.show()

def plotCST(nodeCoordinates, elementNodes, displacements, scale=None,
            elementValues=None, valueLabel='σ (MPa)', cmap='viridis'):
    """Original (dashed) and deformed (solid) CST mesh.
    elementValues : optional (numberElements,) array that shades the deformed elements."""
    elementNodes = np.asarray(elementNodes)
    u = displacements.reshape((-1, 2))                       # (numberNodes, 2): [ux, uy]

    if scale is None:                                        # auto: largest displacement = 10 % of plate size
        size = np.ptp(nodeCoordinates, axis=0).max()
        scale = 0.1*size/np.abs(u).max()
    deformed = nodeCoordinates + scale*u

    fig, ax = plt.subplots(figsize=(7, 6))

    # shaded deformed elements (optional)
    if elementValues is not None:
        pc = PolyCollection([deformed[e] for e in elementNodes], array=np.asarray(elementValues),
                            cmap=cmap, edgecolors='none', alpha=0.85)
        ax.add_collection(pc)
        fig.colorbar(pc, ax=ax, label=valueLabel, shrink=0.8)

    # element outlines
    for e in elementNodes:
        tri = np.append(e, e[0])                             # close the triangle
        ax.plot(*nodeCoordinates[tri].T, 'b--', lw=1.5)
        ax.plot(*deformed[tri].T, 'r-', lw=2)

    ax.plot(*nodeCoordinates.T, 'bo')
    ax.plot(*deformed.T, 'ro')
    for i, (x, y) in enumerate(nodeCoordinates):
        ax.annotate(str(i + 1), (x, y), textcoords='offset points', xytext=(6, 6))

    # element numbers at the centroids
    for i, e in enumerate(elementNodes):
        cx, cy = nodeCoordinates[e].mean(axis=0)
        ax.text(cx, cy, f'({i + 1})', ha='center', va='center', fontsize=9,
                bbox=dict(boxstyle='round,pad=0.15', fc='white', ec='none', alpha=0.7))

    ax.legend(handles=[Line2D([], [], color='b', ls='--', marker='o', label='original'),
                       Line2D([], [], color='r', ls='-', marker='o', label=f'deformed (x{scale:g})')],
              loc='upper center', bbox_to_anchor=(0.5, -0.12), ncol=2)
    ax.autoscale()
    ax.set_aspect('equal')
    ax.grid(True)
    ax.set_xlabel('x (m)'); ax.set_ylabel('y (m)')
    ax.set_title('CST mesh: original and deformed')
    plt.tight_layout()
    plt.show()
    return fig, ax