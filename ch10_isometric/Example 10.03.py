import sys
from pathlib import Path

sys.path.append(str(Path(__file__).resolve().parent.parent))

import numpy as np
from fem.gauss import gaussPoints

I = 0

points, weights = gaussPoints(n = 3)
for i in range(3):
    I += (points[i]**2 + np.cos(points[i]/2))*weights[i]

np.set_printoptions(precision=4)
print(f'I = {I}')

I = 0

for i in range(3):
    I += (3**(points[i]) - points[i])*weights[i]

print(f'I = {I}')