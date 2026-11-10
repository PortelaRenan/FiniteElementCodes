import numpy as np
from math import log10, floor

def eng(x: float, digits: int = 4):
    """
    Format a number as a string in engineering notation (exponent a multiple of 3).

    Parameters:
        x      : Number to format
        digits : Number of decimal places shown in the mantissa (default 4)

    Returns:
        formatted : str
            Number as a string, e.g. "1.2500e-03" or "-3.4200e+00"
    """
    if x == 0:
        return f"{0:.{digits}f}e+00"
    # Exponent rounded down to the nearest multiple of 3
    exp = int(floor(log10(abs(x)) / 3) * 3)
    # Mantissa is x scaled by the exponent; exponent printed with sign and at least 2 digits
    return f"{x / 10**exp:.{digits}f}e{exp:+03d}"

def displacementResults(displacements: np.ndarray):
    """
    Print the nodal displacements and rotations as a table, one row per node,
    in engineering notation.

    Parameters:
        displacements : Global displacement vector ordered [v0, theta0, v1, theta1, ...],
                        where v is the transverse displacement and theta is the rotation
                        of each node. Must be a NumPy array (uses .reshape) with an even
                        number of entries. Assumes every node has both DOFs in the vector.

    Returns:
        None
            Nothing is returned; the table is printed to standard output.
    """
    print(f"{'Node':>4} {'Displacement':>16} {'Rotation':>16}")
    # Each row of the reshaped array holds (displacement, rotation) of one node
    for node, (v, theta) in enumerate(displacements.reshape(-1, 2)):
        # node + 1 makes the printed node numbering start at 1
        print(f"{node + 1:>4} {eng(v):>16} {eng(theta):>16}")